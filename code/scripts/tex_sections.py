#!/usr/bin/env python3
"""Manifest-based LaTeX split/preview/export/apply tool.

The split directories (for example, ``main/`` and ``supplement/``) plus their
``.texsplit.json`` manifests are the source of truth. Root-level ``*.tex``
files are drivers whose managed body is generated from the manifest.
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


MANIFEST_VERSION = 1
MANIFEST_NAME = ".texsplit.json"
BACKUP_DIR = ".texsplit-backup"
ARCHIVE_DIR = ".texsplit-archive"
PREVIEW_DIR = "preview"
CONTENT_BEGIN_RE = re.compile(
    r'^% TEXSPLIT-CONTENT-BEGIN root="(?P<root>[^"]+)"\s*$', re.MULTILINE
)
CONTENT_END_RE = re.compile(r"^% TEXSPLIT-CONTENT-END\s*$", re.MULTILINE)
SOURCE_BEGIN_RE = re.compile(
    r'^% TEXSPLIT:BEGIN file="(?P<file>[^"]+)" id="(?P<id>[^"]+)" '
    r'hash="(?P<hash>sha256:[0-9a-f]{64})"\s*$'
)
SOURCE_END_RE = re.compile(
    r'^% TEXSPLIT:END file="(?P<file>[^"]+)" id="(?P<id>[^"]+)"\s*$'
)
NEW_SECTION_RE = re.compile(
    r'^% TEXSPLIT:NEW-SECTION title="(?P<title>[^"]+)"(?: id="(?P<id>[^"]+)")?\s*$'
)
NEW_SUBSECTION_RE = re.compile(
    r'^% TEXSPLIT:NEW-SUBSECTION section="(?P<section>[^"]+)" '
    r'title="(?P<title>[^"]+)"(?: id="(?P<id>[^"]+)")?\s*$'
)
NEW_END_RE = re.compile(r"^% TEXSPLIT:END-NEW\s*$")
INPUT_RE = re.compile(r"\\input\s*\{(?P<target>[^}]+)\}")
HEADING_START_RE = re.compile(
    r"(?m)^[ \t]*\\(?P<command>section|subsection)\*?"
    r"(?:\s*\[[^\]\n]*\])?\s*\{"
)
TAIL_MARKERS = (
    r"\bibliographystyle",
    r"\bibliography",
    r"\printbibliography",
    r"\onecolumn",
    r"\appendix",
    r"\end{document}",
)
COMMON_ENVIRONMENTS = (
    "document",
    "itemize",
    "enumerate",
    "table",
    "table*",
    "tabular",
    "tabular*",
    "figure",
    "figure*",
    "equation",
    "align",
)


class TexSectionsError(Exception):
    """User-facing command failure."""


@dataclass(frozen=True)
class Heading:
    command: str
    title: str
    start: int
    end: int


@dataclass
class Block:
    file: str
    marker_id: str
    marker_hash: str
    content: str
    start_line: int


@dataclass
class NewBlock:
    kind: str
    title: str
    content: str
    start_line: int
    item_id: str | None = None
    section_id: str | None = None


@dataclass
class ParsedPreview:
    blocks: list[Block]
    new_blocks: list[NewBlock]
    prefix: str
    suffix: str


@dataclass
class CheckResult:
    errors: list[str]
    warnings: list[str]

    def error(self, message: str) -> None:
        self.errors.append(message)

    def warn(self, message: str) -> None:
        self.warnings.append(message)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    project_root = args.project_root.resolve()

    try:
        dispatch(args, project_root)
    except TexSectionsError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Manage manifest-based split LaTeX documents."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("."),
        help="Project directory containing drivers and split roots (default: cwd).",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    split = subparsers.add_parser("split", help="Split a full source into manifest layout.")
    split.add_argument("source", help="Full source .tex file to split.")
    split.add_argument("--root", required=True, help="Split root directory, e.g. main.")
    split.add_argument("--force", action="store_true", help="Replace existing root after backup.")

    for name, help_text in (
        ("preview", "Render editable and clean previews."),
        ("export", "Render clean merged export without markers."),
        ("check", "Validate manifest, files, driver, previews, and rendered LaTeX."),
        ("sync-driver", "Rewrite the driver managed block from manifest."),
        ("status", "Print document structure summary."),
    ):
        sub = subparsers.add_parser(name, help=help_text)
        sub.add_argument("root", help="Document root, e.g. main or supplement.")

    apply_preview = subparsers.add_parser(
        "apply-preview", help="Apply editable preview marker block changes to source files."
    )
    apply_preview.add_argument("root")
    apply_preview.add_argument("--dry-run", action="store_true")
    apply_preview.add_argument("--force", action="store_true")
    apply_preview.add_argument("--only", help="Only apply one project-relative source path.")
    apply_preview.add_argument("--allow-driver", action="store_true")
    apply_preview.add_argument("--adopt-order", action="store_true")
    apply_preview.add_argument("--adopt-new", action="store_true")
    apply_preview.add_argument("--delete-missing", action="store_true")

    add_section = subparsers.add_parser("add-section", help="Add a section.")
    add_section.add_argument("root")
    add_section.add_argument("title")
    add_position_args(add_section)

    add_subsection = subparsers.add_parser("add-subsection", help="Add a subsection.")
    add_subsection.add_argument("root")
    add_subsection.add_argument("section")
    add_subsection.add_argument("title")
    add_position_args(add_subsection)

    remove_section = subparsers.add_parser("remove-section", help="Remove a section.")
    remove_section.add_argument("root")
    remove_section.add_argument("section")
    remove_section.add_argument("--delete", action="store_true")

    remove_subsection = subparsers.add_parser("remove-subsection", help="Remove a subsection.")
    remove_subsection.add_argument("root")
    remove_subsection.add_argument("section")
    remove_subsection.add_argument("subsection")
    remove_subsection.add_argument("--delete", action="store_true")

    rename_section = subparsers.add_parser("rename-section", help="Rename a section heading.")
    rename_section.add_argument("root")
    rename_section.add_argument("section")
    rename_section.add_argument("title")
    rename_section.add_argument("--rename-file", action="store_true")

    rename_subsection = subparsers.add_parser(
        "rename-subsection", help="Rename a subsection heading."
    )
    rename_subsection.add_argument("root")
    rename_subsection.add_argument("section")
    rename_subsection.add_argument("subsection")
    rename_subsection.add_argument("title")
    rename_subsection.add_argument("--rename-file", action="store_true")

    move_section = subparsers.add_parser("move-section", help="Move a section in manifest order.")
    move_section.add_argument("root")
    move_section.add_argument("section")
    add_position_args(move_section, required=True)

    move_subsection = subparsers.add_parser(
        "move-subsection", help="Move a subsection in manifest order."
    )
    move_subsection.add_argument("root")
    move_subsection.add_argument("section")
    move_subsection.add_argument("subsection")
    add_position_args(move_subsection, required=True)

    adopt_section = subparsers.add_parser("adopt-section", help="Adopt an existing section file.")
    adopt_section.add_argument("root")
    adopt_section.add_argument("file")
    adopt_section.add_argument("--title")
    add_position_args(adopt_section)

    adopt_subsection = subparsers.add_parser(
        "adopt-subsection", help="Adopt an existing subsection file."
    )
    adopt_subsection.add_argument("root")
    adopt_subsection.add_argument("section")
    adopt_subsection.add_argument("file")
    adopt_subsection.add_argument("--title")
    add_position_args(adopt_subsection)

    reconcile = subparsers.add_parser("reconcile", help="Report or fix manifest/file drift.")
    reconcile.add_argument("root")
    reconcile.add_argument("--apply-title-from-files", action="store_true")
    reconcile.add_argument("--apply-title-from-manifest", action="store_true")
    reconcile.add_argument("--adopt-orphans", action="store_true")

    return parser


def add_position_args(parser: argparse.ArgumentParser, *, required: bool = False) -> None:
    group = parser.add_mutually_exclusive_group(required=required)
    group.add_argument("--before", help="Insert before this id.")
    group.add_argument("--after", help="Insert after this id.")


def dispatch(args: argparse.Namespace, project_root: Path) -> None:
    command = args.command.replace("-", "_")
    handlers = {
        "split": command_split,
        "preview": command_preview,
        "export": command_export,
        "apply_preview": command_apply_preview,
        "check": command_check,
        "sync_driver": command_sync_driver,
        "add_section": command_add_section,
        "add_subsection": command_add_subsection,
        "remove_section": command_remove_section,
        "remove_subsection": command_remove_subsection,
        "rename_section": command_rename_section,
        "rename_subsection": command_rename_subsection,
        "move_section": command_move_section,
        "move_subsection": command_move_subsection,
        "adopt_section": command_adopt_section,
        "adopt_subsection": command_adopt_subsection,
        "reconcile": command_reconcile,
        "status": command_status,
    }
    handlers[command](project_root, args)


def command_split(project_root: Path, args: argparse.Namespace) -> None:
    source = (project_root / args.source).resolve()
    root_name = args.root
    root_dir = project_root / root_name
    if not source.exists():
        raise TexSectionsError(f"source not found: {source.relative_to(project_root)}")
    if root_dir.exists() and not args.force:
        raise TexSectionsError(f"{root_name}/ already exists; use --force")

    full_text = expand_inputs(project_root, source)
    manifest, file_contents, driver = build_split_manifest(
        project_root, source, root_name, full_text
    )

    backup_dir = create_backup(project_root, [source])
    print(f"backup: {backup_dir.relative_to(project_root)}")
    if root_dir.exists():
        move_path(root_dir, backup_dir / root_name)

    root_dir.mkdir(parents=True, exist_ok=True)
    for rel_file, content in file_contents.items():
        atomic_write_text(root_dir / rel_file, content)
    save_manifest(project_root, manifest)
    atomic_write_text(project_root / manifest["driver"], driver)
    print(f"split: {source.relative_to(project_root)} -> {root_name}/")
    print(f"manifest: {root_name}/{MANIFEST_NAME}")


def command_preview(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    editable_output = project_root / manifest["editable_preview"]
    editable_text = make_preview_paths_portable(
        render_document(project_root, manifest, markers=True)
    )
    clean_output = project_root / manifest["clean_preview"]
    clean_text = make_preview_paths_portable(
        render_document(project_root, manifest, markers=False)
    )
    atomic_write_text(editable_output, editable_text)
    atomic_write_text(clean_output, clean_text)
    print(f"preview: {editable_output.relative_to(project_root)}")
    print(f"preview: {clean_output.relative_to(project_root)}")


def command_export(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    output = project_root / manifest["clean_preview"]
    atomic_write_text(
        output,
        make_preview_paths_portable(render_document(project_root, manifest, markers=False)),
    )
    print(f"export: {output.relative_to(project_root)}")


def command_sync_driver(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    backup_dir = create_backup(project_root, [project_root / manifest["driver"]])
    sync_driver(project_root, manifest)
    print(f"backup: {backup_dir.relative_to(project_root)}")
    print(f"sync-driver: {manifest['driver']}")


def command_apply_preview(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    preview_path = project_root / manifest["editable_preview"]
    if not preview_path.exists():
        raise TexSectionsError(f"editable preview not found: {manifest['editable_preview']}")
    parsed = parse_preview(preview_path.read_text(encoding="utf-8"))
    apply_preview(project_root, manifest, parsed, args)


def command_check(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    result = check_manifest_project(project_root, manifest)
    for warning in result.warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    for error in result.errors:
        print(f"ERROR: {error}", file=sys.stderr)
    if result.errors:
        raise TexSectionsError(f"check failed with {len(result.errors)} error(s)")
    print(f"check: {args.root}")
    print("check complete: no errors")


def command_status(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    print(f"Document: {manifest['root']}")
    print(f"Driver: {manifest['driver']}")
    print(f"Manifest: {manifest['root']}/{MANIFEST_NAME}")
    print(f"Editable preview: {manifest['editable_preview']}")
    print(f"Clean export: {manifest['clean_preview']}")
    print("")
    if manifest.get("backbone"):
        print("Backbone:")
        for fragment in manifest["backbone"]:
            print(f"  {fragment['id']}: {fragment['title']} ({fragment['position']})")
        print("")
    print("Sections:")
    for section in manifest["sections"]:
        print(f"  {section['id']}: {section['title']}")
        for subsection in section.get("subsections", []):
            print(f"    - {subsection['id']}: {subsection['title']}")
    result = check_manifest_project(project_root, manifest)
    if result.warnings:
        print("")
        print("Warnings:")
        for warning in result.warnings:
            print(f"  {warning}")


def command_add_section(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    section = new_section(manifest, args.title)
    insert_item(manifest["sections"], section, before=args.before, after=args.after)
    backup_dir = backup_for_manifest_write(project_root, manifest)
    atomic_write_text(project_root / manifest["root"] / section["file"], section_text(args.title))
    save_manifest(project_root, manifest)
    sync_driver(project_root, manifest)
    print(f"backup: {backup_dir.relative_to(project_root)}")
    print(f"added section: {section['id']}")


def command_add_subsection(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    section = find_section(manifest, args.section)
    subsection = new_subsection(manifest, section, args.title)
    insert_item(section["subsections"], subsection, before=args.before, after=args.after)
    backup_dir = backup_for_manifest_write(project_root, manifest)
    atomic_write_text(
        project_root / manifest["root"] / subsection["file"], subsection_text(args.title)
    )
    save_manifest(project_root, manifest)
    sync_driver(project_root, manifest)
    print(f"backup: {backup_dir.relative_to(project_root)}")
    print(f"added subsection: {subsection['id']}")


def command_remove_section(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    index, section = find_section_with_index(manifest, args.section)
    backup_dir = backup_for_manifest_write(project_root, manifest)
    manifest["sections"].pop(index)
    section_dir = project_root / manifest["root"] / Path(section["file"]).parent
    if args.delete:
        remove_path(section_dir)
    else:
        archive_path(project_root, manifest["root"], section_dir)
    save_manifest(project_root, manifest)
    sync_driver(project_root, manifest)
    print(f"backup: {backup_dir.relative_to(project_root)}")
    print(f"removed section: {section['id']}")


def command_remove_subsection(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    section = find_section(manifest, args.section)
    index, subsection = find_subsection_with_index(section, args.subsection)
    backup_dir = backup_for_manifest_write(project_root, manifest)
    section["subsections"].pop(index)
    path = project_root / manifest["root"] / subsection["file"]
    if args.delete:
        remove_path(path)
    else:
        archive_path(project_root, manifest["root"], path)
    save_manifest(project_root, manifest)
    sync_driver(project_root, manifest)
    print(f"backup: {backup_dir.relative_to(project_root)}")
    print(f"removed subsection: {subsection['id']}")


def command_rename_section(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    section = find_section(manifest, args.section)
    backup_dir = backup_for_manifest_write(project_root, manifest)
    section_path = project_root / manifest["root"] / section["file"]
    replace_heading_in_file(section_path, "section", args.title)
    section["title"] = args.title
    if args.rename_file:
        rename_section_path(project_root, manifest, section)
    save_manifest(project_root, manifest)
    sync_driver(project_root, manifest)
    print(f"backup: {backup_dir.relative_to(project_root)}")
    print(f"renamed section: {section['id']}")


def command_rename_subsection(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    section = find_section(manifest, args.section)
    subsection = find_subsection(section, args.subsection)
    backup_dir = backup_for_manifest_write(project_root, manifest)
    subsection_path = project_root / manifest["root"] / subsection["file"]
    replace_heading_in_file(subsection_path, "subsection", args.title)
    subsection["title"] = args.title
    if args.rename_file:
        rename_subsection_path(project_root, manifest, section, subsection)
    save_manifest(project_root, manifest)
    sync_driver(project_root, manifest)
    print(f"backup: {backup_dir.relative_to(project_root)}")
    print(f"renamed subsection: {subsection['id']}")


def command_move_section(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    backup_dir = backup_for_manifest_write(project_root, manifest)
    move_item(manifest["sections"], args.section, before=args.before, after=args.after)
    save_manifest(project_root, manifest)
    sync_driver(project_root, manifest)
    print(f"backup: {backup_dir.relative_to(project_root)}")
    print(f"moved section: {args.section}")


def command_move_subsection(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    section = find_section(manifest, args.section)
    backup_dir = backup_for_manifest_write(project_root, manifest)
    move_item(section["subsections"], args.subsection, before=args.before, after=args.after)
    save_manifest(project_root, manifest)
    sync_driver(project_root, manifest)
    print(f"backup: {backup_dir.relative_to(project_root)}")
    print(f"moved subsection: {args.subsection}")


def command_adopt_section(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    rel = normalize_doc_relative_path(manifest, args.file)
    path = project_root / manifest["root"] / rel
    if not path.exists():
        raise TexSectionsError(f"section file not found: {manifest['root']}/{rel}")
    title = args.title or read_required_heading(path, "section")
    section_id = unique_id(slugify(title, "section"), existing_section_ids(manifest))
    section = {"id": section_id, "title": title, "file": rel.as_posix(), "subsections": []}
    insert_item(manifest["sections"], section, before=args.before, after=args.after)
    backup_dir = backup_for_manifest_write(project_root, manifest)
    save_manifest(project_root, manifest)
    sync_driver(project_root, manifest)
    print(f"backup: {backup_dir.relative_to(project_root)}")
    print(f"adopted section: {section_id}")


def command_adopt_subsection(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    section = find_section(manifest, args.section)
    rel = normalize_doc_relative_path(manifest, args.file)
    path = project_root / manifest["root"] / rel
    if not path.exists():
        raise TexSectionsError(f"subsection file not found: {manifest['root']}/{rel}")
    title = args.title or read_required_heading(path, "subsection")
    subsection_id = unique_id(slugify(title, "subsection"), existing_subsection_ids(section))
    subsection = {"id": subsection_id, "title": title, "file": rel.as_posix()}
    insert_item(section["subsections"], subsection, before=args.before, after=args.after)
    backup_dir = backup_for_manifest_write(project_root, manifest)
    save_manifest(project_root, manifest)
    sync_driver(project_root, manifest)
    print(f"backup: {backup_dir.relative_to(project_root)}")
    print(f"adopted subsection: {subsection_id}")


def command_reconcile(project_root: Path, args: argparse.Namespace) -> None:
    manifest = load_manifest(project_root, args.root)
    changed = False
    backup_dir: Path | None = None
    report = reconcile_report(project_root, manifest)
    if args.apply_title_from_files or args.apply_title_from_manifest or args.adopt_orphans:
        backup_dir = backup_for_manifest_write(project_root, manifest)
    if args.apply_title_from_files:
        changed |= update_manifest_titles_from_files(project_root, manifest)
    if args.apply_title_from_manifest:
        changed |= update_file_titles_from_manifest(project_root, manifest)
    if args.adopt_orphans:
        changed |= adopt_orphans(project_root, manifest)
    if changed:
        save_manifest(project_root, manifest)
        sync_driver(project_root, manifest)
    if backup_dir:
        print(f"backup: {backup_dir.relative_to(project_root)}")
    for line in report:
        print(line)
    if changed:
        print("reconcile: applied changes")
    else:
        print("reconcile: report only")


def build_split_manifest(
    project_root: Path, source: Path, root_name: str, text: str
) -> tuple[dict[str, Any], dict[Path, str], str]:
    headings = [heading for heading in find_headings(text) if heading.command == "section"]
    if not headings:
        raise TexSectionsError(f"no top-level \\section found in {source.relative_to(project_root)}")
    first = headings[0]
    body_end = find_body_end(text, first.start)
    prefix = text[: first.start]
    body = text[first.start:body_end]
    suffix = text[body_end:]
    body_sections = [heading for heading in find_headings(body) if heading.command == "section"]

    manifest = {
        "version": MANIFEST_VERSION,
        "entry": source.name,
        "root": root_name,
        "driver": source.name,
        "editable_preview": f"{PREVIEW_DIR}/{root_name}.editable.tex",
        "clean_preview": f"{PREVIEW_DIR}/{root_name}.clean.tex",
        "backbone": [
            {
                "id": "front-matter",
                "title": "Front matter",
                "file": "backbone/front-matter.tex",
                "position": "before-sections",
            },
            {
                "id": "back-matter",
                "title": "Back matter",
                "file": "backbone/back-matter.tex",
                "position": "after-sections",
            },
        ],
        "sections": [],
    }
    files: dict[Path, str] = {
        Path("backbone/front-matter.tex"): ensure_final_newline(prefix),
        Path("backbone/back-matter.tex"): ensure_final_newline(suffix),
    }
    used_section_ids: set[str] = set()

    for index, heading in enumerate(body_sections):
        next_start = body_sections[index + 1].start if index + 1 < len(body_sections) else len(body)
        section_block = body[heading.start:next_start]
        section_id = unique_id(slugify(heading.title, "section"), used_section_ids)
        used_section_ids.add(section_id)
        section_file = Path(section_id) / "main.tex"
        section = {
            "id": section_id,
            "title": heading.title,
            "file": section_file.as_posix(),
            "subsections": [],
        }

        subsection_headings = [
            sub for sub in find_headings(section_block) if sub.command == "subsection"
        ]
        if subsection_headings:
            section_intro = section_block[: subsection_headings[0].start]
            files[section_file] = ensure_final_newline(section_intro)
            used_sub_ids: set[str] = set()
            for sub_index, sub_heading in enumerate(subsection_headings):
                sub_next = (
                    subsection_headings[sub_index + 1].start
                    if sub_index + 1 < len(subsection_headings)
                    else len(section_block)
                )
                sub_id = unique_id(slugify(sub_heading.title, "subsection"), used_sub_ids)
                used_sub_ids.add(sub_id)
                sub_file = Path(section_id) / f"{sub_id}.tex"
                files[sub_file] = ensure_final_newline(section_block[sub_heading.start:sub_next])
                section["subsections"].append(
                    {"id": sub_id, "title": sub_heading.title, "file": sub_file.as_posix()}
                )
        else:
            files[section_file] = ensure_final_newline(section_block)
        manifest["sections"].append(section)

    driver_text = render_driver_block(manifest)
    return manifest, files, driver_text


def load_manifest(project_root: Path, root_name: str) -> dict[str, Any]:
    path = project_root / root_name / MANIFEST_NAME
    if not path.exists():
        raise TexSectionsError(f"manifest not found: {root_name}/{MANIFEST_NAME}")
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise TexSectionsError(f"invalid manifest JSON: {path}: {exc}") from exc
    validate_manifest_shape(manifest)
    return manifest


def save_manifest(project_root: Path, manifest: dict[str, Any]) -> None:
    path = project_root / manifest["root"] / MANIFEST_NAME
    atomic_write_text(path, json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")


def validate_manifest_shape(manifest: dict[str, Any]) -> None:
    required = {
        "version",
        "entry",
        "root",
        "driver",
        "editable_preview",
        "clean_preview",
        "sections",
    }
    missing = sorted(required - set(manifest))
    if missing:
        raise TexSectionsError(f"manifest missing fields: {', '.join(missing)}")
    if manifest["version"] != MANIFEST_VERSION:
        raise TexSectionsError(f"unsupported manifest version: {manifest['version']}")
    if not isinstance(manifest["sections"], list):
        raise TexSectionsError("manifest sections must be a list")
    if "backbone" in manifest and not isinstance(manifest["backbone"], list):
        raise TexSectionsError("manifest backbone must be a list")


def iter_document_items(
    manifest: dict[str, Any],
) -> Iterable[tuple[dict[str, Any], str, dict[str, Any] | None]]:
    for fragment in manifest.get("backbone", []):
        if fragment.get("position") == "before-sections":
            yield fragment, "fragment", None
    for section in manifest["sections"]:
        yield section, "section", None
        for subsection in section.get("subsections", []):
            yield subsection, "subsection", section
    for fragment in manifest.get("backbone", []):
        if fragment.get("position") == "after-sections":
            yield fragment, "fragment", None


def render_driver_block(manifest: dict[str, Any]) -> str:
    root = manifest["root"]
    lines = [
        f'% TEXSPLIT-CONTENT-BEGIN root="{root}"',
        "% Generated by scripts/tex_sections.py. Do not edit this block manually.",
    ]
    for item, _kind, _section in iter_document_items(manifest):
        lines.append(rf"\input{{{without_tex_suffix(root, item['file'])}}}")
    lines.append("% TEXSPLIT-CONTENT-END")
    return "\n".join(lines) + "\n"


def sync_driver(project_root: Path, manifest: dict[str, Any]) -> None:
    driver_path = project_root / manifest["driver"]
    driver_text = driver_path.read_text(encoding="utf-8") if driver_path.exists() else ""
    atomic_write_text(driver_path, replace_content_block(driver_text, manifest))


def replace_content_block(driver_text: str, manifest: dict[str, Any]) -> str:
    block = render_driver_block(manifest).rstrip()
    if manifest.get("backbone"):
        return block + "\n"
    region = find_content_region(driver_text)
    if region:
        start, end = region
        return f"{driver_text[:start].rstrip()}\n\n{block}\n\n{driver_text[end:].lstrip()}".rstrip() + "\n"
    insert_at = find_driver_insert_position(driver_text)
    if insert_at is None:
        return f"{driver_text.rstrip()}\n\n{block}\n"
    return (
        f"{driver_text[:insert_at].rstrip()}\n\n{block}\n\n{driver_text[insert_at:].lstrip()}".rstrip()
        + "\n"
    )


def find_content_region(text: str) -> tuple[int, int] | None:
    begin = None
    for match in CONTENT_BEGIN_RE.finditer(text):
        begin = match.start()
        break
    if begin is None:
        return None
    end_match = CONTENT_END_RE.search(text, begin)
    if not end_match:
        raise TexSectionsError("driver has TEXSPLIT-CONTENT-BEGIN without END")
    line_end = text.find("\n", end_match.end())
    end = len(text) if line_end == -1 else line_end + 1
    return begin, end


def find_driver_insert_position(text: str) -> int | None:
    positions = [text.find(marker) for marker in TAIL_MARKERS if text.find(marker) != -1]
    if positions:
        return min(positions)
    end_doc = text.find(r"\end{document}")
    return end_doc if end_doc != -1 else None


def render_document(project_root: Path, manifest: dict[str, Any], *, markers: bool) -> str:
    body = render_body(project_root, manifest, markers=markers).rstrip()
    begin_line = f'% TEXSPLIT-CONTENT-BEGIN root="{manifest["root"]}"'
    end_line = "% TEXSPLIT-CONTENT-END"
    replacement = f"{begin_line}\n{body}\n{end_line}\n"
    if manifest.get("backbone"):
        return replacement.rstrip() + "\n"

    driver_path = project_root / manifest["driver"]
    if not driver_path.exists():
        raise TexSectionsError(f"driver not found: {manifest['driver']}")
    driver = driver_path.read_text(encoding="utf-8")
    region = find_content_region(driver)
    if not region:
        driver = replace_content_block(driver, manifest)
        region = find_content_region(driver)
    assert region is not None
    start, end = region
    return f"{driver[:start].rstrip()}\n\n{replacement}\n{driver[end:].lstrip()}".rstrip() + "\n"


def make_preview_paths_portable(text: str) -> str:
    """Let generated preview files compile from project root or preview/ cwd."""
    text = re.sub(
        r"(?m)^(?P<indent>\s*)\\bibliography\s*\{(?P<files>[^}]*)\}(?P<tail>[^\n]*)$",
        replace_bibliography_for_preview,
        text,
    )
    text = re.sub(
        r"(?m)^(?P<indent>\s*)\\graphicspath\s*\{(?P<body>.*)\}(?P<tail>\s*(?:%.*)?)$",
        replace_graphicspath_for_preview,
        text,
    )
    return text


def replace_bibliography_for_preview(match: re.Match[str]) -> str:
    indent = match.group("indent")
    files = match.group("files")
    tail = match.group("tail")
    entries = [entry.strip() for entry in files.split(",") if entry.strip()]
    if not entries:
        return match.group(0)
    fallback = ",".join(preview_parent_path(entry) for entry in entries)
    probe = bibliography_probe(entries[0])
    return (
        f"{indent}\\IfFileExists{{{probe}}}"
        f"{{\\bibliography{{{files}}}}}"
        f"{{\\bibliography{{{fallback}}}}}{tail}"
    )


def replace_graphicspath_for_preview(match: re.Match[str]) -> str:
    body = match.group("body")
    paths = re.findall(r"\{([^{}]+)\}", body)
    if not paths:
        return match.group(0)
    merged: list[str] = []
    for path in paths:
        if path not in merged:
            merged.append(path)
        parent = preview_parent_path(path)
        if parent != path and parent not in merged:
            merged.append(parent)
    formatted = "".join(f"{{{path}}}" for path in merged)
    return f"{match.group('indent')}\\graphicspath{{{formatted}}}{match.group('tail')}"


def bibliography_probe(entry: str) -> str:
    stripped = entry.strip()
    return stripped if stripped.endswith(".bib") else f"{stripped}.bib"


def preview_parent_path(path: str) -> str:
    stripped = path.strip()
    if (
        not stripped
        or stripped.startswith("../")
        or stripped.startswith("/")
        or re.match(r"^[A-Za-z]+:", stripped)
    ):
        return stripped
    if stripped.startswith("./"):
        stripped = stripped[2:]
    return f"../{stripped}"


def render_body(project_root: Path, manifest: dict[str, Any], *, markers: bool) -> str:
    chunks: list[str] = []
    for item, kind, _section in iter_document_items(manifest):
        chunks.append(render_source_file(project_root, manifest, item, kind, markers))
    return "\n".join(chunk.rstrip() for chunk in chunks if chunk is not None) + "\n"


def render_source_file(
    project_root: Path,
    manifest: dict[str, Any],
    item: dict[str, Any],
    kind: str,
    markers: bool,
) -> str:
    rel = Path(manifest["root"]) / item["file"]
    path = project_root / rel
    content = ensure_final_newline(path.read_text(encoding="utf-8"))
    if not markers:
        return content
    marker_id = f"{kind}:{item['id']}"
    digest = sha256_text(content)
    return (
        f'% TEXSPLIT:BEGIN file="{rel.as_posix()}" id="{marker_id}" hash="{digest}"\n'
        f"{content}"
        f'% TEXSPLIT:END file="{rel.as_posix()}" id="{marker_id}"\n'
    )


def parse_preview(text: str) -> ParsedPreview:
    lines = text.splitlines(keepends=True)
    blocks: list[Block] = []
    new_blocks: list[NewBlock] = []
    prefix_parts: list[str] = []
    suffix_parts: list[str] = []
    current: dict[str, Any] | None = None
    content_parts: list[str] = []
    seen_first_block = False
    after_last_block_index = 0

    for line_number, line in enumerate(lines, start=1):
        stripped = line.rstrip("\n")
        begin = SOURCE_BEGIN_RE.match(stripped)
        new_section = NEW_SECTION_RE.match(stripped)
        new_subsection = NEW_SUBSECTION_RE.match(stripped)
        if current is None and (begin or new_section or new_subsection):
            current = {
                "kind": "source" if begin else ("new-section" if new_section else "new-subsection"),
                "start_line": line_number,
                "match": begin or new_section or new_subsection,
            }
            content_parts = []
            if not seen_first_block:
                seen_first_block = True
                prefix_parts = lines[: line_number - 1]
            continue
        if current is not None:
            if current["kind"] == "source":
                end = SOURCE_END_RE.match(stripped)
                if end:
                    match = current["match"]
                    if end.group("file") != match.group("file") or end.group("id") != match.group("id"):
                        raise TexSectionsError(f"mismatched marker end at line {line_number}")
                    blocks.append(
                        Block(
                            file=match.group("file"),
                            marker_id=match.group("id"),
                            marker_hash=match.group("hash"),
                            content="".join(content_parts),
                            start_line=current["start_line"],
                        )
                    )
                    current = None
                    after_last_block_index = line_number
                    continue
            else:
                if NEW_END_RE.match(stripped):
                    match = current["match"]
                    new_blocks.append(
                        NewBlock(
                            kind=current["kind"],
                            title=match.group("title"),
                            item_id=match.groupdict().get("id"),
                            section_id=match.groupdict().get("section"),
                            content="".join(content_parts),
                            start_line=current["start_line"],
                        )
                    )
                    current = None
                    after_last_block_index = line_number
                    continue
            content_parts.append(line)
            continue
    if current is not None:
        raise TexSectionsError(f"unclosed TEXSPLIT marker at line {current['start_line']}")
    if seen_first_block:
        suffix_parts = lines[after_last_block_index:]
    return ParsedPreview(blocks=blocks, new_blocks=new_blocks, prefix="".join(prefix_parts), suffix="".join(suffix_parts))


def apply_preview(
    project_root: Path,
    manifest: dict[str, Any],
    parsed: ParsedPreview,
    args: argparse.Namespace,
) -> None:
    duplicate = duplicate_values([block.file for block in parsed.blocks])
    if duplicate:
        raise TexSectionsError(f"duplicate source markers: {', '.join(duplicate)}")
    source_items = source_items_by_project_file(manifest)
    block_by_file = {block.file: block for block in parsed.blocks}
    missing = [path for path in source_items if path not in block_by_file]
    if missing and not args.delete_missing:
        for path in missing:
            print(f"WARNING: source block missing from editable preview: {path}", file=sys.stderr)

    changed_files: list[str] = []
    title_changes: list[str] = []
    manifest_changed = False
    target_only = args.only
    if target_only is not None:
        target_only = project_relative_source_path(manifest, target_only)

    backup_paths: list[Path] = [project_root / manifest["root"] / MANIFEST_NAME]
    pending_writes: list[tuple[Path, str]] = []
    pending_archives: list[Path] = []
    driver_text: str | None = None
    for block in parsed.blocks:
        if target_only and block.file != target_only:
            continue
        if block.file not in source_items:
            raise TexSectionsError(f"preview references unknown source file: {block.file}")
        path = project_root / block.file
        old = ensure_final_newline(path.read_text(encoding="utf-8"))
        old_hash = sha256_text(old)
        if old_hash != block.marker_hash and not args.force:
            raise TexSectionsError(
                "source file changed after editable preview was generated: "
                f"{block.file}. Use --force to overwrite."
            )
        content = ensure_final_newline(block.content)
        item, kind, section = source_items[block.file]
        if kind in {"section", "subsection"}:
            heading_title = first_heading_title(content, "section" if kind == "section" else "subsection")
            if heading_title and heading_title != item["title"]:
                title_changes.append(f"{block.marker_id}: {item['title']} -> {heading_title}")
                item["title"] = heading_title
                manifest_changed = True
        if content != old:
            backup_paths.append(path)
            changed_files.append(block.file)
            pending_writes.append((path, content))

    if args.adopt_new:
        new_backup_paths, new_writes = adopt_new_blocks(
            project_root, manifest, parsed.new_blocks, dry_run=args.dry_run
        )
        backup_paths.extend(new_backup_paths)
        pending_writes.extend(new_writes)
        manifest_changed = manifest_changed or bool(parsed.new_blocks)
    elif parsed.new_blocks:
        for block in parsed.new_blocks:
            print(
                f"WARNING: ignored new preview block at line {block.start_line}; "
                "rerun with --adopt-new",
                file=sys.stderr,
            )

    order_changed = update_order_from_blocks(manifest, parsed.blocks, adopt=args.adopt_order)
    if order_changed:
        manifest_changed = True
    if missing and args.delete_missing:
        missing_backup_paths, archives = delete_missing_blocks(
            project_root, manifest, missing, dry_run=args.dry_run
        )
        backup_paths.extend(missing_backup_paths)
        pending_archives.extend(archives)
        manifest_changed = True
    elif detect_order_change(manifest, parsed.blocks):
        print("WARNING: detected source block order changes; rerun with --adopt-order", file=sys.stderr)

    if args.allow_driver:
        backup_paths.append(project_root / manifest["driver"])
        driver_text = render_driver_from_preview(project_root, manifest, parsed)

    if manifest_changed or changed_files:
        backup_paths.append(project_root / manifest["driver"])

    if changed_files or manifest_changed or args.allow_driver:
        if not args.dry_run:
            backup_dir = create_backup(project_root, backup_paths)
            print(f"backup: {backup_dir.relative_to(project_root)}")
    if args.dry_run:
        print("dry-run: no files written")
    else:
        for path, content in pending_writes:
            atomic_write_text(path, content)
        for path in pending_archives:
            archive_path(project_root, manifest["root"], path)
        if manifest_changed:
            save_manifest(project_root, manifest)
        if args.allow_driver and driver_text is not None:
            atomic_write_text(project_root / manifest["driver"], driver_text)
        elif changed_files or manifest_changed:
            sync_driver(project_root, manifest)
    for path in changed_files:
        print(f"updated: {path}")
    for change in title_changes:
        print(f"updated title: {change}")
    if not changed_files and not title_changes and not manifest_changed:
        print("apply-preview: no changes")


def adopt_new_blocks(
    project_root: Path,
    manifest: dict[str, Any],
    new_blocks: list[NewBlock],
    *,
    dry_run: bool,
) -> tuple[list[Path], list[tuple[Path, str]]]:
    backups = [project_root / manifest["root"] / MANIFEST_NAME]
    writes: list[tuple[Path, str]] = []
    for block in new_blocks:
        if block.kind == "new-section":
            section_id = unique_id(
                block.item_id or slugify(block.title, "section"), existing_section_ids(manifest)
            )
            section = {
                "id": section_id,
                "title": block.title,
                "file": f"{section_id}/main.tex",
                "subsections": [],
            }
            manifest["sections"].append(section)
            path = project_root / manifest["root"] / section["file"]
        else:
            if not block.section_id:
                raise TexSectionsError(f"new subsection at line {block.start_line} has no section")
            section = find_section(manifest, block.section_id)
            subsection_id = unique_id(
                block.item_id or slugify(block.title, "subsection"),
                existing_subsection_ids(section),
            )
            subsection = {
                "id": subsection_id,
                "title": block.title,
                "file": f"{Path(section['file']).parent.as_posix()}/{subsection_id}.tex",
            }
            section["subsections"].append(subsection)
            path = project_root / manifest["root"] / subsection["file"]
        backups.append(path)
        writes.append((path, ensure_final_newline(block.content)))
        print(f"adopted new {block.kind.replace('new-', '')}: {path.relative_to(project_root)}")
    return backups, writes


def delete_missing_blocks(
    project_root: Path,
    manifest: dict[str, Any],
    missing: list[str],
    *,
    dry_run: bool,
) -> tuple[list[Path], list[Path]]:
    backups = [project_root / manifest["root"] / MANIFEST_NAME]
    archives: list[Path] = []
    for project_rel in missing:
        item = remove_manifest_item_by_project_file(manifest, project_rel)
        path = project_root / project_rel
        backups.append(path)
        if path.exists():
            archives.append(path)
        print(f"removed missing block from manifest: {item['id']}")
    return backups, archives


def render_driver_from_preview(
    project_root: Path, manifest: dict[str, Any], parsed: ParsedPreview
) -> str:
    if manifest.get("backbone"):
        return render_driver_block(manifest)
    driver_path = project_root / manifest["driver"]
    old_driver = driver_path.read_text(encoding="utf-8") if driver_path.exists() else ""
    new_driver = parsed.prefix.rstrip() + "\n\n" + render_driver_block(manifest).rstrip() + "\n\n" + parsed.suffix.lstrip()
    if r"\begin{document}" in old_driver and r"\begin{document}" not in new_driver:
        raise TexSectionsError("--allow-driver would remove \\begin{document}")
    if r"\end{document}" in old_driver and r"\end{document}" not in new_driver:
        raise TexSectionsError("--allow-driver would remove \\end{document}")
    return new_driver.rstrip() + "\n"


def update_driver_from_preview(
    project_root: Path, manifest: dict[str, Any], parsed: ParsedPreview
) -> None:
    atomic_write_text(
        project_root / manifest["driver"],
        render_driver_from_preview(project_root, manifest, parsed),
    )


def check_manifest_project(project_root: Path, manifest: dict[str, Any]) -> CheckResult:
    result = CheckResult(errors=[], warnings=[])
    try:
        validate_manifest_semantics(project_root, manifest, result)
        check_driver_sync(project_root, manifest, result)
        clean = render_document(project_root, manifest, markers=False)
        check_environment_counts(manifest["root"], clean, result)
        check_end_document_is_terminal(manifest["root"], clean, result)
        check_preview_stale(project_root, manifest, result)
    except TexSectionsError as exc:
        result.error(str(exc))
    return result


def validate_manifest_semantics(project_root: Path, manifest: dict[str, Any], result: CheckResult) -> None:
    section_ids: set[str] = set()
    backbone_ids: set[str] = set()
    referenced: set[Path] = set()
    root_dir = project_root / manifest["root"]
    for fragment in manifest.get("backbone", []):
        fragment_id = fragment.get("id")
        if not fragment_id:
            result.error("backbone fragment missing id")
        elif fragment_id in backbone_ids:
            result.error(f"duplicate backbone id: {fragment_id}")
        else:
            backbone_ids.add(fragment_id)
        if fragment.get("position") not in {"before-sections", "after-sections"}:
            result.error(f"invalid backbone position for {fragment_id}: {fragment.get('position')}")
        check_item_file(project_root, manifest, fragment, "fragment", result, referenced)
    for section in manifest["sections"]:
        if section["id"] in section_ids:
            result.error(f"duplicate section id: {section['id']}")
        section_ids.add(section["id"])
        check_item_file(project_root, manifest, section, "section", result, referenced)
        sub_ids: set[str] = set()
        for subsection in section.get("subsections", []):
            if subsection["id"] in sub_ids:
                result.error(f"duplicate subsection id in {section['id']}: {subsection['id']}")
            sub_ids.add(subsection["id"])
            check_item_file(project_root, manifest, subsection, "subsection", result, referenced)
    if root_dir.exists():
        for tex in root_dir.rglob("*.tex"):
            if tex.name == MANIFEST_NAME:
                continue
            try:
                rel = tex.relative_to(root_dir)
            except ValueError:
                continue
            if rel not in referenced:
                result.warn(f"orphan file found: {manifest['root']}/{rel.as_posix()}")


def check_item_file(
    project_root: Path,
    manifest: dict[str, Any],
    item: dict[str, Any],
    kind: str,
    result: CheckResult,
    referenced: set[Path],
) -> None:
    rel = Path(item["file"])
    referenced.add(rel)
    path = project_root / manifest["root"] / rel
    if not path.exists():
        result.error(f"missing file: {manifest['root']}/{rel.as_posix()}")
        return
    if kind == "fragment":
        return
    expected_command = "section" if kind == "section" else "subsection"
    title = first_heading_title(path.read_text(encoding="utf-8"), expected_command)
    if title is None:
        result.error(f"missing \\{expected_command} heading: {manifest['root']}/{rel.as_posix()}")
    elif title != item["title"]:
        result.error(
            "heading mismatch: "
            f"{manifest['root']}/{rel.as_posix()} manifest={item['title']!r} file={title!r}"
        )


def check_driver_sync(project_root: Path, manifest: dict[str, Any], result: CheckResult) -> None:
    driver_path = project_root / manifest["driver"]
    if not driver_path.exists():
        result.error(f"driver not found: {manifest['driver']}")
        return
    text = driver_path.read_text(encoding="utf-8")
    region = find_content_region(text)
    if not region:
        result.error("driver managed TEXSPLIT-CONTENT block missing")
        return
    start, end = region
    current = text[start:end].strip()
    expected = render_driver_block(manifest).strip()
    if current != expected:
        result.error("driver managed block is out of sync; run sync-driver")
    if manifest.get("backbone") and (text[:start] + text[end:]).strip():
        result.error("driver contains non-split content outside managed TEXSPLIT-CONTENT block")


def check_preview_stale(project_root: Path, manifest: dict[str, Any], result: CheckResult) -> None:
    preview = project_root / manifest["editable_preview"]
    if not preview.exists():
        result.warn(f"editable preview not found: {manifest['editable_preview']}")
        return
    try:
        parsed = parse_preview(preview.read_text(encoding="utf-8"))
    except TexSectionsError as exc:
        result.error(f"editable preview marker error: {exc}")
        return
    duplicate = duplicate_values([block.file for block in parsed.blocks])
    for value in duplicate:
        result.error(f"duplicate marker for source file: {value}")
    source_items = source_items_by_project_file(manifest)
    for block in parsed.blocks:
        path = project_root / block.file
        if block.file not in source_items:
            result.warn(f"preview references unknown source file: {block.file}")
            continue
        if not path.exists():
            result.error(f"preview references missing source file: {block.file}")
            continue
        if sha256_text(ensure_final_newline(path.read_text(encoding="utf-8"))) != block.marker_hash:
            result.warn(f"editable preview is stale for {block.file}")


def reconcile_report(project_root: Path, manifest: dict[str, Any]) -> list[str]:
    result = check_manifest_project(project_root, manifest)
    lines = [f"Document: {manifest['root']}"]
    lines.extend(f"WARNING: {warning}" for warning in result.warnings)
    lines.extend(f"ERROR: {error}" for error in result.errors)
    if len(lines) == 1:
        lines.append("OK: no drift detected")
    return lines


def update_manifest_titles_from_files(project_root: Path, manifest: dict[str, Any]) -> bool:
    changed = False
    for section in manifest["sections"]:
        title = first_heading_title((project_root / manifest["root"] / section["file"]).read_text(encoding="utf-8"), "section")
        if title and title != section["title"]:
            section["title"] = title
            changed = True
        for subsection in section.get("subsections", []):
            sub_title = first_heading_title((project_root / manifest["root"] / subsection["file"]).read_text(encoding="utf-8"), "subsection")
            if sub_title and sub_title != subsection["title"]:
                subsection["title"] = sub_title
                changed = True
    return changed


def update_file_titles_from_manifest(project_root: Path, manifest: dict[str, Any]) -> bool:
    changed = False
    for section in manifest["sections"]:
        path = project_root / manifest["root"] / section["file"]
        changed = replace_heading_in_file(path, "section", section["title"]) or changed
        for subsection in section.get("subsections", []):
            sub_path = project_root / manifest["root"] / subsection["file"]
            changed = replace_heading_in_file(sub_path, "subsection", subsection["title"]) or changed
    return changed


def adopt_orphans(project_root: Path, manifest: dict[str, Any]) -> bool:
    root_dir = project_root / manifest["root"]
    known = {Path(item[0]) for item in source_items_by_doc_file(manifest)}
    changed = False
    for tex in sorted(root_dir.rglob("*.tex")):
        rel = tex.relative_to(root_dir)
        if rel in known:
            continue
        if tex.name == "main.tex" and first_heading_title(tex.read_text(encoding="utf-8"), "section"):
            title = first_heading_title(tex.read_text(encoding="utf-8"), "section")
            assert title is not None
            section_id = unique_id(slugify(title, "section"), existing_section_ids(manifest))
            manifest["sections"].append(
                {"id": section_id, "title": title, "file": rel.as_posix(), "subsections": []}
            )
            changed = True
        elif first_heading_title(tex.read_text(encoding="utf-8"), "subsection"):
            parent = rel.parent.as_posix()
            section = next(
                (sec for sec in manifest["sections"] if Path(sec["file"]).parent.as_posix() == parent),
                None,
            )
            if section:
                title = first_heading_title(tex.read_text(encoding="utf-8"), "subsection")
                assert title is not None
                subsection_id = unique_id(slugify(title, "subsection"), existing_subsection_ids(section))
                section["subsections"].append(
                    {"id": subsection_id, "title": title, "file": rel.as_posix()}
                )
                changed = True
    return changed


def find_headings(text: str) -> list[Heading]:
    headings: list[Heading] = []
    for match in HEADING_START_RE.finditer(text):
        open_brace = match.end() - 1
        close_brace = find_matching_brace(text, open_brace)
        if close_brace is None:
            raise TexSectionsError(f"unclosed heading brace near byte offset {match.start()}")
        headings.append(
            Heading(
                command=match.group("command"),
                title=text[open_brace + 1 : close_brace].strip(),
                start=match.start(),
                end=close_brace + 1,
            )
        )
    return headings


def find_matching_brace(text: str, open_brace: int) -> int | None:
    depth = 0
    index = open_brace
    while index < len(text):
        char = text[index]
        if char == "\\":
            index += 2
            continue
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return index
        index += 1
    return None


def find_body_end(text: str, first_section_start: int) -> int:
    positions = []
    for marker in TAIL_MARKERS:
        position = text.find(marker, first_section_start)
        if position != -1:
            positions.append(position)
    return min(positions) if positions else len(text)


def first_heading_title(text: str, command: str) -> str | None:
    for heading in find_headings(text):
        if heading.command == command:
            return heading.title
    return None


def read_required_heading(path: Path, command: str) -> str:
    title = first_heading_title(path.read_text(encoding="utf-8"), command)
    if title is None:
        raise TexSectionsError(f"{path} does not contain \\{command} heading")
    return title


def replace_heading_in_file(path: Path, command: str, title: str) -> bool:
    text = path.read_text(encoding="utf-8")
    headings = [heading for heading in find_headings(text) if heading.command == command]
    if not headings:
        raise TexSectionsError(f"{path} does not contain \\{command} heading")
    heading = headings[0]
    replacement = f"\\{command}{{{title}}}"
    new_text = text[: heading.start] + replacement + text[heading.end :]
    if new_text == text:
        return False
    atomic_write_text(path, new_text)
    return True


def expand_inputs(project_root: Path, source: Path, stack: tuple[Path, ...] = ()) -> str:
    source = source.resolve()
    if source in stack:
        chain = " -> ".join(path.relative_to(project_root).as_posix() for path in (*stack, source))
        raise TexSectionsError(f"cyclic input detected: {chain}")
    text = source.read_text(encoding="utf-8")
    lines: list[str] = []
    for line in text.splitlines(keepends=True):
        active, comment = split_active_comment(line)
        stripped = active.strip()

        def replace(match: re.Match[str]) -> str:
            target = match.group("target").strip()
            resolved = resolve_input(project_root, source, target)
            if resolved is None:
                raise TexSectionsError(
                    f"missing input {target!r} referenced from {source.relative_to(project_root)}"
                )
            return expand_inputs(project_root, resolved, (*stack, source))

        if re.fullmatch(r"\\input\s*\{[^}]+\}", stripped):
            expanded = INPUT_RE.sub(replace, stripped)
            if line.endswith("\n") and not expanded.endswith("\n"):
                expanded += "\n"
            lines.append(expanded)
        else:
            lines.append(INPUT_RE.sub(replace, active) + comment)
    return "".join(lines)


def resolve_input(project_root: Path, source: Path, target: str) -> Path | None:
    raw = Path(target)
    bases = [raw] if raw.is_absolute() else [project_root / raw, source.parent / raw]
    candidates: list[Path] = []
    for base in bases:
        candidates.append(base)
        if base.suffix != ".tex":
            candidates.append(base.with_suffix(".tex"))
    seen: set[Path] = set()
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if resolved.exists() and resolved.is_file():
            return resolved
    return None


def split_active_comment(line: str) -> tuple[str, str]:
    escaped = False
    for index, char in enumerate(line):
        if char == "\\":
            escaped = not escaped
            continue
        if char == "%" and not escaped:
            return line[:index], line[index:]
        escaped = False
    return line, ""


def check_environment_counts(name: str, text: str, result: CheckResult) -> None:
    for env in COMMON_ENVIRONMENTS:
        escaped = re.escape(env)
        begin = len(re.findall(rf"\\begin\{{{escaped}\}}", text))
        end = len(re.findall(rf"\\end\{{{escaped}\}}", text))
        if begin != end:
            result.error(f"{name}: environment {env} has {begin} begin(s) and {end} end(s)")


def check_end_document_is_terminal(name: str, text: str, result: CheckResult) -> None:
    matches = list(re.finditer(r"\\end\{document\}", text))
    if not matches:
        result.error(f"{name}: missing \\end{{document}}")
        return
    if len(matches) > 1:
        result.error(f"{name}: multiple \\end{{document}} commands")
        return
    tail_lines = [
        line.strip()
        for line in text[matches[0].end() :].splitlines()
        if line.strip() and not line.strip().startswith("% TEXSPLIT")
    ]
    if tail_lines:
        result.error(f"{name}: non-whitespace content appears after \\end{{document}}")


def new_section(manifest: dict[str, Any], title: str) -> dict[str, Any]:
    section_id = unique_id(slugify(title, "section"), existing_section_ids(manifest))
    return {
        "id": section_id,
        "title": title,
        "file": f"{section_id}/main.tex",
        "subsections": [],
    }


def new_subsection(
    manifest: dict[str, Any], section: dict[str, Any], title: str
) -> dict[str, Any]:
    subsection_id = unique_id(slugify(title, "subsection"), existing_subsection_ids(section))
    section_dir = Path(section["file"]).parent
    return {
        "id": subsection_id,
        "title": title,
        "file": (section_dir / f"{subsection_id}.tex").as_posix(),
    }


def section_text(title: str) -> str:
    return f"\\section{{{title}}}\n\n% TODO\n"


def subsection_text(title: str) -> str:
    return f"\\subsection{{{title}}}\n\n% TODO\n"


def insert_item(
    items: list[dict[str, Any]], item: dict[str, Any], *, before: str | None, after: str | None
) -> None:
    if before:
        index = find_item_index(items, before)
        items.insert(index, item)
    elif after:
        index = find_item_index(items, after)
        items.insert(index + 1, item)
    else:
        items.append(item)


def move_item(items: list[dict[str, Any]], item_id: str, *, before: str | None, after: str | None) -> None:
    index = find_item_index(items, item_id)
    item = items.pop(index)
    insert_item(items, item, before=before, after=after)


def find_item_index(items: list[dict[str, Any]], item_id: str) -> int:
    for index, item in enumerate(items):
        if item["id"] == item_id:
            return index
    raise TexSectionsError(f"id not found: {item_id}")


def find_section(manifest: dict[str, Any], section_id: str) -> dict[str, Any]:
    return manifest["sections"][find_item_index(manifest["sections"], section_id)]


def find_section_with_index(
    manifest: dict[str, Any], section_id: str
) -> tuple[int, dict[str, Any]]:
    index = find_item_index(manifest["sections"], section_id)
    return index, manifest["sections"][index]


def find_subsection(section: dict[str, Any], subsection_id: str) -> dict[str, Any]:
    return section["subsections"][find_item_index(section["subsections"], subsection_id)]


def find_subsection_with_index(
    section: dict[str, Any], subsection_id: str
) -> tuple[int, dict[str, Any]]:
    index = find_item_index(section["subsections"], subsection_id)
    return index, section["subsections"][index]


def rename_section_path(project_root: Path, manifest: dict[str, Any], section: dict[str, Any]) -> None:
    root_dir = project_root / manifest["root"]
    old_dir = root_dir / Path(section["file"]).parent
    new_id = unique_id(slugify(section["title"], "section"), existing_section_ids(manifest) - {section["id"]})
    new_dir = root_dir / new_id
    if old_dir == new_dir:
        return
    move_path(old_dir, new_dir)
    section["file"] = f"{new_id}/main.tex"
    for subsection in section.get("subsections", []):
        subsection["file"] = f"{new_id}/{Path(subsection['file']).name}"


def rename_subsection_path(
    project_root: Path,
    manifest: dict[str, Any],
    section: dict[str, Any],
    subsection: dict[str, Any],
) -> None:
    old = project_root / manifest["root"] / subsection["file"]
    new_id = unique_id(
        slugify(subsection["title"], "subsection"),
        existing_subsection_ids(section) - {subsection["id"]},
    )
    new = old.parent / f"{new_id}.tex"
    if old == new:
        return
    move_path(old, new)
    subsection["file"] = new.relative_to(project_root / manifest["root"]).as_posix()


def update_order_from_blocks(
    manifest: dict[str, Any], blocks: list[Block], *, adopt: bool
) -> bool:
    if not adopt:
        return False
    existing_sections = {section["id"]: section for section in manifest["sections"]}
    section_order: list[str] = []
    subsection_order: dict[str, list[str]] = {}
    current_section: str | None = None
    subsection_to_section = {
        subsection["id"]: section["id"]
        for section in manifest["sections"]
        for subsection in section.get("subsections", [])
    }
    for block in blocks:
        kind, item_id = block.marker_id.split(":", 1)
        if kind == "section" and item_id in existing_sections:
            current_section = item_id
            if item_id not in section_order:
                section_order.append(item_id)
            subsection_order.setdefault(item_id, [])
        elif kind == "subsection":
            owner = current_section or subsection_to_section.get(item_id)
            if owner:
                subsection_order.setdefault(owner, [])
                if item_id not in subsection_order[owner]:
                    subsection_order[owner].append(item_id)
    changed = False
    if set(section_order) == set(existing_sections) and section_order != [s["id"] for s in manifest["sections"]]:
        manifest["sections"] = [existing_sections[item_id] for item_id in section_order]
        changed = True
    for section in manifest["sections"]:
        order = subsection_order.get(section["id"])
        if not order:
            continue
        existing = {sub["id"]: sub for sub in section.get("subsections", [])}
        if set(order) == set(existing) and order != [sub["id"] for sub in section["subsections"]]:
            section["subsections"] = [existing[item_id] for item_id in order]
            changed = True
    return changed


def detect_order_change(manifest: dict[str, Any], blocks: list[Block]) -> bool:
    clone = json.loads(json.dumps(manifest))
    return update_order_from_blocks(clone, blocks, adopt=True)


def source_items_by_project_file(
    manifest: dict[str, Any],
) -> dict[str, tuple[dict[str, Any], str, dict[str, Any] | None]]:
    root = manifest["root"]
    mapping: dict[str, tuple[dict[str, Any], str, dict[str, Any] | None]] = {}
    for item, kind, section in iter_document_items(manifest):
        mapping[f"{root}/{item['file']}"] = (item, kind, section)
    return mapping


def source_items_by_doc_file(manifest: dict[str, Any]) -> list[tuple[str, dict[str, Any], str]]:
    items: list[tuple[str, dict[str, Any], str]] = []
    for item, kind, _section in iter_document_items(manifest):
        items.append((item["file"], item, kind))
    return items


def remove_manifest_item_by_project_file(manifest: dict[str, Any], project_rel: str) -> dict[str, Any]:
    root = manifest["root"]
    for index, fragment in enumerate(list(manifest.get("backbone", []))):
        if f"{root}/{fragment['file']}" == project_rel:
            return manifest["backbone"].pop(index)
    for section_index, section in enumerate(list(manifest["sections"])):
        if f"{root}/{section['file']}" == project_rel:
            return manifest["sections"].pop(section_index)
        for sub_index, subsection in enumerate(list(section.get("subsections", []))):
            if f"{root}/{subsection['file']}" == project_rel:
                return section["subsections"].pop(sub_index)
    raise TexSectionsError(f"source file is not in manifest: {project_rel}")


def existing_section_ids(manifest: dict[str, Any]) -> set[str]:
    return {section["id"] for section in manifest["sections"]}


def existing_subsection_ids(section: dict[str, Any]) -> set[str]:
    return {subsection["id"] for subsection in section.get("subsections", [])}


def slugify(title: str, fallback: str) -> str:
    text = re.sub(r"\\[a-zA-Z]+(?:\[[^\]]*\])?(?:\{[^}]*\})?", " ", title)
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()
    text = text.lower().replace("'", "")
    text = re.sub(r"[^a-z0-9]+", "-", text)
    text = re.sub(r"-{2,}", "-", text).strip("-")
    return text or fallback


def unique_id(base: str, existing: set[str]) -> str:
    candidate = base
    suffix = 2
    while candidate in existing:
        candidate = f"{base}-{suffix}"
        suffix += 1
    return candidate


def without_tex_suffix(root: str, rel_file: str) -> str:
    path = Path(root) / rel_file
    if path.suffix == ".tex":
        path = path.with_suffix("")
    return path.as_posix()


def project_relative_source_path(manifest: dict[str, Any], path: str) -> str:
    if path.startswith(f"{manifest['root']}/"):
        return path
    return f"{manifest['root']}/{path}"


def normalize_doc_relative_path(manifest: dict[str, Any], path: str) -> Path:
    rel = Path(path)
    if rel.parts and rel.parts[0] == manifest["root"]:
        rel = Path(*rel.parts[1:])
    if rel.is_absolute() or ".." in rel.parts:
        raise TexSectionsError(f"path must stay inside {manifest['root']}/: {path}")
    return rel


def ensure_final_newline(text: str) -> str:
    return text if text.endswith("\n") else f"{text}\n"


def sha256_text(text: str) -> str:
    return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()


def duplicate_values(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    dup: set[str] = set()
    for value in values:
        if value in seen:
            dup.add(value)
        seen.add(value)
    return sorted(dup)


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
        temp = Path(handle.name)
    os.replace(temp, path)


def timestamp() -> str:
    return _datetime.datetime.now().strftime("%Y%m%d-%H%M%S")


def create_backup(project_root: Path, paths: Iterable[Path]) -> Path:
    backup = project_root / BACKUP_DIR / timestamp()
    suffix = 1
    while backup.exists():
        backup = project_root / BACKUP_DIR / f"{timestamp()}-{suffix}"
        suffix += 1
    backup.mkdir(parents=True, exist_ok=True)
    seen: set[Path] = set()
    for path in paths:
        if not path.exists():
            continue
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        rel = resolved.relative_to(project_root)
        target = backup / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if resolved.is_dir():
            shutil.copytree(resolved, target)
        else:
            shutil.copy2(resolved, target)
    return backup


def backup_for_manifest_write(project_root: Path, manifest: dict[str, Any]) -> Path:
    paths = [project_root / manifest["root"] / MANIFEST_NAME, project_root / manifest["driver"]]
    return create_backup(project_root, paths)


def archive_path(project_root: Path, root_name: str, path: Path) -> None:
    if not path.exists():
        return
    rel = path.relative_to(project_root)
    target = project_root / ARCHIVE_DIR / timestamp() / rel
    target.parent.mkdir(parents=True, exist_ok=True)
    move_path(path, target)


def move_path(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        raise TexSectionsError(f"target already exists: {target}")
    shutil.move(str(source), str(target))


def remove_path(path: Path) -> None:
    if not path.exists():
        return
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


if __name__ == "__main__":
    raise SystemExit(main())
