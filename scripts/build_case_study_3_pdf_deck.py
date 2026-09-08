#!/usr/bin/env python3
"""
Build multi-page PDFs of case_study_3 PNGs for one or more CSP CSV pair lists.

Each page is one target's ``{HOLO}_case_study_3.png``. Output basenames include the
CSV stem so decks are distinguishable, e.g.:

  outputs/case_study_decks/case_study_3_CSP_UBQ_ph0.5_temp5C.pdf
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from PIL import Image

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from scripts.target_resolution import canonical_output_dir_name  # noqa: E402

DEFAULT_SEARCH_ROOTS = ("outputs",)
DEFAULT_CSVS = (
    "data/CSP_UBQ.csv",
    "data/CSP_UBQ_ph0.5_temp5C.csv",
    "data/CSP_UBQ_ph0.5_temp5C_ideal_sequence_match.csv",
    "data/CSP_UBQ_ph0.5_temp5C_non_ideal_sequence_match.csv",
)


def _find_case_insensitive(directory: Path, filename: str) -> Optional[Path]:
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


def resolve_case_study_3_png(
    repo: Path,
    row: Dict[str, str],
    *,
    search_roots: Sequence[str] = DEFAULT_SEARCH_ROOTS,
) -> Optional[Path]:
    holo = (row.get("holo_pdb") or "").strip()
    apo = (row.get("apo_bmrb") or "").strip()
    if not holo or not apo:
        return None
    logical = canonical_output_dir_name(holo, apo)
    name = f"{holo.upper()}_case_study_3.png"
    for sub in search_roots:
        p = _find_case_insensitive(repo / sub / logical, name)
        if p is not None:
            return p
        # Also try original holo casing in filename.
        p = _find_case_insensitive(repo / sub / logical, f"{holo}_case_study_3.png")
        if p is not None:
            return p
    return None


def _prepare_page(png_path: Path, *, max_width: int) -> Image.Image:
    im = Image.open(png_path).convert("RGB")
    if im.width > max_width:
        new_h = max(1, int(round(im.height * (max_width / im.width))))
        im = im.resize((max_width, new_h), Image.Resampling.LANCZOS)
    return im


def build_pdf_for_csv(
    csv_path: Path,
    out_pdf: Path,
    *,
    repo: Path,
    search_roots: Sequence[str],
    max_width: int,
) -> Tuple[int, int, List[str]]:
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    pages: List[Image.Image] = []
    missing: List[str] = []
    for row in rows:
        holo = (row.get("holo_pdb") or "").strip()
        apo = (row.get("apo_bmrb") or "").strip()
        label = canonical_output_dir_name(holo, apo) if holo and apo else (holo or "?")
        png = resolve_case_study_3_png(repo, row, search_roots=search_roots)
        if png is None:
            missing.append(label)
            continue
        pages.append(_prepare_page(png, max_width=max_width))

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    if not pages:
        raise RuntimeError(f"No case_study_3 PNGs found for {csv_path}")

    first, rest = pages[0], pages[1:]
    first.save(out_pdf, "PDF", save_all=True, append_images=rest, resolution=150.0)
    for im in pages:
        im.close()
    return len(rows), len(pages), missing


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--csv",
        action="append",
        default=None,
        help="CSP CSV path (repeatable). Default: the four standard sets.",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=_REPO / "outputs" / "case_study_decks",
        help="Directory for output PDFs",
    )
    ap.add_argument(
        "--search-roots",
        default=",".join(DEFAULT_SEARCH_ROOTS),
        help="Comma-separated output roots relative to repo (search order)",
    )
    ap.add_argument(
        "--max-width",
        type=int,
        default=2400,
        help="Downscale pages wider than this (pixels) before PDF write",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    csv_paths = [Path(p) for p in (args.csv or DEFAULT_CSVS)]
    search_roots = [s.strip() for s in args.search_roots.split(",") if s.strip()]

    n_fail = 0
    for csv_path in csv_paths:
        if not csv_path.is_file():
            print(f"ERROR: missing CSV {csv_path}", file=sys.stderr)
            n_fail += 1
            continue
        stem = csv_path.stem
        out_pdf = args.out_dir / f"case_study_3_{stem}.pdf"
        try:
            n_rows, n_pages, missing = build_pdf_for_csv(
                csv_path,
                out_pdf,
                repo=_REPO,
                search_roots=search_roots,
                max_width=args.max_width,
            )
        except Exception as exc:
            print(f"ERROR: {csv_path}: {exc}", file=sys.stderr)
            n_fail += 1
            continue
        print(
            f"Wrote {out_pdf}  pages={n_pages}/{n_rows} rows  "
            f"missing={len(missing)}"
        )
        if missing:
            preview = ", ".join(missing[:8])
            more = "" if len(missing) <= 8 else f" (+{len(missing) - 8} more)"
            print(f"  missing PNGs: {preview}{more}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
