"""Discovery and loading for the optional compiled event backend."""

from __future__ import annotations

import csv
import hashlib
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from junction_filter import JunctionEvidence
from apadetect_native import executable_path


DEFAULT_NATIVE_EVENTS = executable_path()


def _native_executable_loads(path: Path) -> bool:
    """Return whether the OS can start the optional native executable."""
    try:
        completed = subprocess.run(
            [str(path)],
            text=True,
            capture_output=True,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return completed.returncode not in (126, 127) and completed.returncode >= 0


def resolve_event_backend(
    requested: str,
    native_path: Path | None = None,
) -> tuple[str, Path]:
    path = (native_path or DEFAULT_NATIVE_EVENTS).resolve()
    available = path.is_file() and os.access(path, os.X_OK)
    if requested == "python":
        return "python", path
    if requested == "native":
        if not available:
            raise FileNotFoundError(
                f"Native event backend is not executable: {path}"
            )
        return "native", path
    if requested != "auto":
        raise ValueError(f"Unsupported event backend: {requested}")
    return ("native" if available and _native_executable_loads(path) else "python"), path


def native_binary_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_meta(path: Path) -> dict[str, int]:
    with path.open() as handle:
        return {
            row[0]: int(row[1])
            for row in csv.reader(handle, delimiter="\t")
        }


def collect_coverage_events_native(
    executable: Path,
    bam_path: Path,
    chromosome: str,
    strand_mode: str,
):
    with tempfile.TemporaryDirectory(prefix="apadetect_events_") as temporary:
        root = Path(temporary)
        completed = subprocess.run(
            [
                str(executable),
                str(bam_path),
                chromosome,
                strand_mode,
                str(root),
            ],
            text=True,
            capture_output=True,
        )
        if completed.returncode:
            raise RuntimeError(
                f"Native event backend failed with exit {completed.returncode}: "
                f"{completed.stderr.strip()}"
            )
        meta = _read_meta(root / "meta.tsv")
        events = {
            strand: (
                np.fromfile(root / f"{name}.positions.bin", dtype=np.int64),
                np.fromfile(root / f"{name}.deltas.bin", dtype=np.int32),
            )
            for strand, name in (("+", "plus"), ("-", "minus"))
        }
        expected = {
            "+": meta["plus_events"],
            "-": meta["minus_events"],
        }
        for strand in ("+", "-"):
            if (
                len(events[strand][0]) != expected[strand]
                or len(events[strand][1]) != expected[strand]
            ):
                raise RuntimeError(
                    f"Native event output is truncated for strand {strand}"
                )
        evidence = JunctionEvidence()
        with (root / "junctions.tsv").open() as handle:
            for row in csv.reader(handle, delimiter="\t"):
                evidence.records[(int(row[0]), int(row[1]), row[2])] = [
                    int(row[3]),
                    int(row[4]),
                    int(row[5]),
                ]
        if len(evidence.records) != meta["junctions"]:
            raise RuntimeError("Native junction output is truncated")
        return meta["chromosome_length"], events, evidence
