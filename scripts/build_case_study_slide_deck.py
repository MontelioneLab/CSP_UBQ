#!/usr/bin/env python3
"""Build a PowerPoint of case-study PNGs for CSV rows (optionally filtered by holo_pdb)."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

from PIL import Image
from pptx import Presentation
from pptx.util import Inches, Pt


_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from scripts.target_resolution import canonical_output_dir_name  # noqa: E402

DEFAULT_SEARCH_ROOTS = (
    "outputs",
)


def _repo_root() -> Path:
    return _REPO


def _canonical_dir(row: Dict[str, str]) -> str | None:
    holo = (row.get("holo_pdb") or "").strip()
    apo = (row.get("apo_bmrb") or "").strip()
    if not holo or not apo:
        return None
    return canonical_output_dir_name(holo, apo)


def _parse_holo_pdbs(raw: str | None) -> set[str] | None:
    if not raw:
        return None
    items = {s.strip().upper() for s in raw.split(",") if s.strip()}
    return items or None


def _find_case_insensitive(directory: Path, filename: str) -> Path | None:
    """Return directory/filename if present, else a case-insensitive match."""
    exact = directory / filename
    if exact.is_file():
        return exact
    target = filename.lower()
    try:
        for p in directory.iterdir():
            if p.is_file() and p.name.lower() == target:
                return p
    except FileNotFoundError:
        return None
    return None


def resolve_case_study_png(
    repo: Path,
    row: Dict[str, str],
    *,
    version: int = 1,
    search_roots: Sequence[str] = DEFAULT_SEARCH_ROOTS,
) -> Path | None:
    holo_pdb = (row.get("holo_pdb") or "").strip()
    logical_dir = _canonical_dir(row)
    if not holo_pdb or logical_dir is None:
        return None
    if version == 1:
        name = f"{holo_pdb}_case_study.png"
    elif version == 2:
        name = f"{holo_pdb}_case_study_2.png"
    elif version == 3:
        name = f"{holo_pdb.upper()}_case_study_3.png"
    else:
        raise ValueError(f"Unsupported case-study version: {version}")

    for sub in search_roots:
        p = _find_case_insensitive(repo / sub / logical_dir, name)
        if p is not None:
            return p
    return None


def add_image_slide(prs: Presentation, png_path: Path) -> None:
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    im = Image.open(png_path).convert("RGB")
    iw, ih = im.size
    sw = float(prs.slide_width)
    sh = float(prs.slide_height)
    aspect_img = iw / ih
    aspect_slide = sw / sh
    if aspect_img > aspect_slide:
        width = sw
        height = sw / aspect_img
    else:
        height = sh
        width = sh * aspect_img
    left = int((sw - width) / 2)
    top = int((sh - height) / 2)
    slide.shapes.add_picture(str(png_path), left, top, width=int(width), height=int(height))


def add_missing_slide(
    prs: Presentation,
    pdb: str,
    logical_dir: str,
    *,
    version: int,
    search_roots: Iterable[str],
) -> None:
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    box = slide.shapes.add_textbox(Inches(0.5), Inches(3), Inches(12), Inches(2))
    tf = box.text_frame
    if version == 1:
        suffix = "_case_study.png"
    elif version == 2:
        suffix = "_case_study_2.png"
    else:
        suffix = "_case_study_3.png"
    tf.text = f"holo_pdb: {pdb} (case_study v{version})"
    p = tf.paragraphs[0]
    p.font.size = Pt(28)
    p.font.bold = True
    roots = ", ".join(f"{r}/{logical_dir}/" for r in search_roots)
    p2 = tf.add_paragraph()
    p2.text = f"Figure not found (expected {pdb}{suffix} under {roots})."
    p2.font.size = Pt(18)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--csv",
        type=Path,
        default=_repo_root() / "data" / "CSP_UBQ_ph0.5_temp5C.csv",
        help="CSV with holo_pdb / apo_bmrb columns",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=_repo_root() / "tmp" / "slides" / "ph05_temp5C_case_study" / "output.pptx",
        help="Output .pptx path",
    )
    ap.add_argument(
        "--holo-pdbs",
        type=str,
        default=None,
        help="Comma-separated holo_pdb IDs to include (case-insensitive). "
        "When set, every matching apo pair from the CSV is included.",
    )
    ap.add_argument(
        "--include-v2",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include case_study_2.png slides after each v1 slide (default: true). Ignored if --version is set.",
    )
    ap.add_argument(
        "--version",
        type=int,
        default=None,
        choices=[1, 2, 3],
        help="Only include this case-study version (1, 2, or 3). Overrides --include-v2.",
    )
    ap.add_argument(
        "--search-roots",
        type=str,
        default=",".join(DEFAULT_SEARCH_ROOTS),
        help="Comma-separated output root dirs relative to repo (search order).",
    )
    ap.add_argument(
        "--skip-missing",
        action="store_true",
        help="Omit rows whose case-study PNG is missing (default: insert a placeholder slide).",
    )
    args = ap.parse_args()
    repo = _repo_root()
    search_roots = [s.strip() for s in args.search_roots.split(",") if s.strip()]
    holo_filter = _parse_holo_pdbs(args.holo_pdbs)
    if args.version is not None:
        versions = [args.version]
    else:
        versions = [1, 2] if args.include_v2 else [1]

    args.out.parent.mkdir(parents=True, exist_ok=True)

    with args.csv.open(newline="", encoding="utf-8") as f:
        rows: List[Dict[str, str]] = list(csv.DictReader(f))

    prs = Presentation()
    prs.slide_width = Inches(13.333333)
    prs.slide_height = Inches(7.5)

    n_rows = 0
    for row in rows:
        pdb = (row.get("holo_pdb") or "").strip()
        if not pdb:
            continue
        if holo_filter is not None and pdb.upper() not in holo_filter:
            continue
        n_rows += 1
        logical = _canonical_dir(row) or pdb
        for version in versions:
            png = resolve_case_study_png(
                repo, row, version=version, search_roots=search_roots
            )
            if png is not None:
                add_image_slide(prs, png)
            elif not args.skip_missing:
                add_missing_slide(
                    prs, pdb, logical, version=version, search_roots=search_roots
                )

    prs.save(args.out)
    print(
        f"Wrote {len(prs.slides)} slides from {n_rows} CSV rows "
        f"(versions={versions}) -> {args.out.resolve()}"
    )


if __name__ == "__main__":
    main()
