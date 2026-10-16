#!/usr/bin/env python3
"""Stage an explicit, compact public release without publishing or copying caches.

Nothing is overwritten. Media, the recorded scene and supporting diagnostic
folders must be selected explicitly. The generated manifest contains hashes of
the distributed bytes, including the unmodified upstream solver source.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parent
UPSTREAM_COMMIT = "9986522d08b4d887e5ef9a53e12cecae4e376d22"
SKIP_PARTS = {".git", ".venv", "__pycache__", ".mpl-cache", ".pytest_cache", ".vscode", ".DS_Store"}
ROOT_SUFFIXES = {".py", ".sh", ".md", ".toml", ".txt", ".text"}
DIAGNOSTIC_NAMES = {"results.json", "protocol.json", "profiles.npz", "assessment.json", "summary.json", "grid-assessment.json", "cap512-diagnosis.json", "cap512-step-check.json", "RESULTS.md", "README.md", "run-result.json", "config.toml"}
DIAGNOSTIC_SUFFIXES = {".png", ".pdf"}
MAX_GIT_BYTES = 100_000_000


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        h = hashlib.file_digest(stream, "sha256")
    return h.hexdigest()


def walk_files(folder: Path):
    for path in sorted(folder.rglob("*")):
        relative = path.relative_to(folder)
        if any(part in SKIP_PARTS for part in relative.parts):
            continue
        if path.is_symlink():
            raise ValueError(f"Symlink must not enter the release: {path}")
        if path.is_file():
            yield path, relative


def make_selection(args):
    selected: dict[Path, Path] = {}

    def add(source: Path, target: Path):
        source = source.absolute()
        if source.is_symlink() or not source.is_file():
            raise ValueError(f"Expected regular input file: {source}")
        if target.is_absolute() or ".." in target.parts:
            raise ValueError(f"Invalid target path: {target}")
        if source.stat().st_size >= MAX_GIT_BYTES:
            raise ValueError(f"GitHub regular-file limit exceeded: {source}")
        if target in selected and selected[target] != source:
            raise ValueError(f"Conflicting release path: {target}")
        selected[target] = source

    for source in sorted(ROOT.iterdir()):
        if source.is_file() and (source.suffix in ROOT_SUFFIXES or source.name in {".gitignore", "LICENSE", "NOTICE"}):
            add(source, Path(source.name))
    for name in ("src", "tests", "research", "fonts", "data", "release"):
        if not (ROOT / name).exists():
            continue
        for source, relative in walk_files(ROOT / name):
            add(source, Path(name) / relative)
    upstream = ROOT / "third_party/coffee-src"
    for source, relative in walk_files(upstream):
        add(source, Path("third_party/coffee-src") / relative)
    for source in sorted((ROOT / "third_party").iterdir()):
        if source.is_file() and (source.suffix in {".patch", ".json"}):
            add(source, Path("third_party") / source.name)
    for source in args.media:
        add(source, Path("media") / source.name)
    if args.scene:
        for source, relative in walk_files(args.scene):
            if source.name in {"scene.npz", "scene.json"}:
                add(source, Path("scene") / relative)
    for folder in args.diagnostic:
        for source, relative in walk_files(folder):
            if source.name in DIAGNOSTIC_NAMES or source.suffix in DIAGNOSTIC_SUFFIXES:
                add(source, Path("diagnostics") / folder.name / relative)
    if args.snapshot:
        if args.snapshot.stat().st_size > 25_000_000:
            raise ValueError("Optional field snapshot exceeds the compact 25 MB allowance")
        add(args.snapshot, Path("recorded-inputs") / args.snapshot.name)
    for source in args.record:
        add(source, Path("recorded-inputs") / source.name)
    return selected


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=ROOT.parent / "sideratlas-pulsar-v1")
    ap.add_argument("--media", type=Path, action="append", default=[], help="Final video or preview file; repeatable")
    ap.add_argument("--scene", type=Path, help="Directory containing recorded scene.npz and scene.json")
    ap.add_argument("--diagnostic", type=Path, action="append", default=[], help="Explicit diagnostic folder; excludes raw trajectories")
    ap.add_argument("--snapshot", type=Path, help="Optional compact canonical field input, at most 25 MB")
    ap.add_argument("--record", type=Path, action="append", default=[], help="Other explicitly selected provenance input")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    selection = make_selection(args)
    total = sum(p.stat().st_size for p in selection.values())
    largest = sorted(((p.stat().st_size, str(target)) for target, p in selection.items()), reverse=True)[:10]
    summary = dict(file_count=len(selection), total_bytes=total, largest_files=largest,
                   output=str(args.output.absolute()), dry_run=args.dry_run)
    if args.dry_run:
        print(json.dumps(summary, indent=2))
        return
    if not args.media or not args.scene:
        ap.error("A release requires explicitly selected --media and --scene; use --dry-run for source-only inspection")
    output = args.output.absolute()
    if output.exists():
        raise FileExistsError(f"Release directory already exists: {output}")
    output.mkdir(parents=True)
    records = []
    for target, source in sorted(selection.items()):
        destination = output / target
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        records.append(dict(path=str(target), size_bytes=destination.stat().st_size, sha256=digest(destination)))
    manifest = dict(version="1.0.0", repository="https://github.com/Loganpis/sideratlas-pulsar",
        scientific_scope="Vela-informed visualization with disclosed approximation and observational disagreement",
        quantitative_vela_reconstruction=False, upstream=dict(url="https://github.com/fizban007/CoffeeGPU", commit=UPSTREAM_COMMIT),
        excluded="Environments, compiled binaries, build trees, raw run histories and historical output collections",
        files=records, file_count=len(records), total_bytes=sum(r["size_bytes"] for r in records))
    (output / "release-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary | {"manifest": str(output / "release-manifest.json")}, indent=2))


if __name__ == "__main__":
    main()
