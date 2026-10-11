#!/usr/bin/env python3
"""Resolve a run directory without allowing writes outside the workflow root."""

from __future__ import annotations

from pathlib import Path


def _relative(value: str, label: str) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"{label} must be a relative path without '..': {value}")
    return path


def _reject_symlinks(root: Path, target: Path) -> None:
    current = root
    for part in target.relative_to(root).parts:
        current = current / part
        if current.exists() and current.is_symlink():
            raise ValueError(f"Managed output path contains a symlink: {current}")


def managed_run_dir(root: Path, output_root: str, output_dir: str, mode: str) -> Path:
    root = root.resolve()
    managed = (root / _relative(output_root, "output_root")).resolve()
    run = (managed / _relative(output_dir, "output_dir")).resolve()
    if managed != root and root not in managed.parents:
        raise ValueError("output_root escapes the workflow root")
    if run != managed and managed not in run.parents:
        raise ValueError("output_dir escapes output_root")
    _reject_symlinks(root, managed)
    _reject_symlinks(root, run)
    if mode not in {"fresh", "resume"}:
        raise ValueError("output_mode must be 'fresh' or 'resume'")
    if mode == "fresh" and run.exists() and any(run.iterdir()):
        raise ValueError(f"Fresh run refuses non-empty output directory: {run}")
    return run
