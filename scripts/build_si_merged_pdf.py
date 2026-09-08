#!/usr/bin/env python3
"""
Compile SI table/figure TeX masters and merge into SI_documents/SI merged.pdf.

Requires latexmk (or pdflatex) and pypdf.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent

TABLES_MAIN = "main"
TABLES_OUT = "CSP_UBQ_SUPPL_TABLES.pdf"
TOC_MASTER = "si_toc"
TOC_OUT = "SI_TOC.pdf"

TABLE_MASTERS: list[tuple[str, str]] = [
    ("si_table_st1", "CSP_UBQ_TABLE.pdf"),
    ("si_tables_ec", "CSP_UBQ_SUPPL_TABLES_EC_CLASSES.pdf"),
    ("si_tables_scope", "CSP_UBQ_SUPPL_TABLES_SCOPE.pdf"),
    ("si_tables_domain", "CSP_UBQ_SUPPL_TABLES_DOMAIN.pdf"),
    ("si_tables_dissimilar", "CSP_dissimilar_conditions.pdf"),
]

C_TERM_MASTER = "si_c_term_figures"
C_TERM_OUT = "SI C term.pdf"
N_TERM = "SI N term.pdf"
SUPPLEMENTARY_TEXT = "Supplementary Text.pdf"
REFERENCES = "SI_C_term_references.pdf"
MERGED = "SI merged.pdf"

# Heuristic SF16 page index (0-based) in a prior C-term PDF when PNG is missing.
DEFAULT_SF16_PAGE_INDEX = 6
# N-term page 0 = title; page 1 = stale TOC (dropped); pages 2+ = narrative.
N_TERM_TITLE_PAGES = 1
N_TERM_DROP_AFTER_TITLE = 1


def _which_latex() -> list[str]:
    latexmk = shutil.which("latexmk")
    if latexmk:
        return [latexmk, "-pdf", "-interaction=nonstopmode", "-halt-on-error"]
    pdflatex = shutil.which("pdflatex")
    if pdflatex:
        return [pdflatex, "-interaction=nonstopmode", "-halt-on-error"]
    texbin = Path("/Library/TeX/texbin")
    if (texbin / "latexmk").is_file():
        return [str(texbin / "latexmk"), "-pdf", "-interaction=nonstopmode", "-halt-on-error"]
    if (texbin / "pdflatex").is_file():
        return [str(texbin / "pdflatex"), "-interaction=nonstopmode", "-halt-on-error"]
    raise FileNotFoundError(
        "Neither latexmk nor pdflatex found. Install MacTeX/TeX Live before building SI PDFs."
    )


def _compile_tex(master_stem: str, tex_dir: Path) -> Path:
    cmd_base = _which_latex()
    tex_path = tex_dir / f"{master_stem}.tex"
    if not tex_path.is_file():
        raise FileNotFoundError(f"Missing TeX master: {tex_path}")

    if Path(cmd_base[0]).name == "latexmk":
        subprocess.run(cmd_base + [tex_path.name], cwd=str(tex_dir), check=True)
    else:
        cmd = cmd_base + [tex_path.name]
        subprocess.run(cmd, cwd=str(tex_dir), check=True)
        subprocess.run(cmd, cwd=str(tex_dir), check=True)

    pdf_path = tex_dir / f"{master_stem}.pdf"
    if not pdf_path.is_file():
        raise RuntimeError(f"LaTeX did not produce {pdf_path}")
    return pdf_path


def _ensure_references_pdf(si_dir: Path, c_term_pdf: Path) -> Path:
    """Ensure SI_C_term_references.pdf exists (last page of a prior C-term)."""
    from pypdf import PdfReader, PdfWriter

    refs = si_dir / REFERENCES
    if refs.is_file():
        return refs
    if not c_term_pdf.is_file():
        raise FileNotFoundError(
            f"Missing references PDF ({refs}) and no C-term to extract from ({c_term_pdf})."
        )
    reader = PdfReader(str(c_term_pdf))
    if not reader.pages:
        raise RuntimeError(f"Empty C-term PDF: {c_term_pdf}")
    writer = PdfWriter()
    writer.add_page(reader.pages[-1])
    with refs.open("wb") as f:
        writer.write(f)
    print(f"[SI PDF] Extracted references page -> {refs}")
    return refs


def _extract_sf16_png_from_c_term(
    c_term_pdf: Path,
    out_png: Path,
    page_index: int,
) -> bool:
    if out_png.is_file():
        return True
    if not c_term_pdf.is_file():
        return False

    pdftoppm = shutil.which("pdftoppm")
    if pdftoppm:
        stem = out_png.with_suffix("")
        subprocess.run(
            [
                pdftoppm,
                "-png",
                "-f",
                str(page_index + 1),
                "-l",
                str(page_index + 1),
                "-singlefile",
                str(c_term_pdf),
                str(stem),
            ],
            check=True,
        )
        return out_png.is_file()

    try:
        from pypdf import PdfReader, PdfWriter
    except ImportError:
        return False

    reader = PdfReader(str(c_term_pdf))
    if page_index < 0 or page_index >= len(reader.pages):
        return False
    writer = PdfWriter()
    writer.add_page(reader.pages[page_index])
    tmp_pdf = out_png.with_suffix(".pdf")
    with tmp_pdf.open("wb") as f:
        writer.write(f)
    print(f"[SI PDF] SF16 PNG unavailable; saved page PDF to {tmp_pdf}")
    return False


def _build_c_term_pdf(
    tex_dir: Path,
    si_dir: Path,
    figures_dir: Path,
    *,
    sf16_page_index: int,
) -> Path:
    from pypdf import PdfReader, PdfWriter

    existing_c = si_dir / C_TERM_OUT
    refs_pdf = _ensure_references_pdf(si_dir, existing_c)

    sf16_png = figures_dir / "SF16_pdb_search.png"
    if not sf16_png.is_file() and existing_c.is_file():
        _extract_sf16_png_from_c_term(existing_c, sf16_png, sf16_page_index)

    fig_pdf = _compile_tex(C_TERM_MASTER, tex_dir)
    writer = PdfWriter()
    for page in PdfReader(str(fig_pdf)).pages:
        writer.add_page(page)

    sf16_pdf = figures_dir / "SF16_pdb_search.pdf"
    if not sf16_png.is_file() and sf16_pdf.is_file():
        insert_at = min(6, len(writer.pages))
        sf16_pages = PdfReader(str(sf16_pdf)).pages
        pages = list(writer.pages)
        new_writer = PdfWriter()
        for i, page in enumerate(pages):
            if i == insert_at:
                for sp in sf16_pages:
                    new_writer.add_page(sp)
            new_writer.add_page(page)
        if insert_at >= len(pages):
            for sp in sf16_pages:
                new_writer.add_page(sp)
        writer = new_writer

    # Append dedicated references only (never re-append old SF19).
    for page in PdfReader(str(refs_pdf)).pages:
        writer.add_page(page)
    print(f"[SI PDF] Appended references from {refs_pdf.name}")

    out = si_dir / C_TERM_OUT
    with out.open("wb") as f:
        writer.write(f)
    return out


def _merge_n_term_with_toc(si_dir: Path, tex_dir: Path) -> list:
    """Return pages: N-term title + TOC + Supplementary Text + leftover N-term narrative."""
    from pypdf import PdfReader

    n_term = si_dir / N_TERM
    if not n_term.is_file():
        raise FileNotFoundError(f"Missing N-term PDF: {n_term}")
    supp_text = si_dir / SUPPLEMENTARY_TEXT
    if not supp_text.is_file():
        raise FileNotFoundError(f"Missing Supplementary Text PDF: {supp_text}")

    print(f"[SI PDF] Compiling {TOC_MASTER}.tex -> {TOC_OUT}")
    toc_built = _compile_tex(TOC_MASTER, tex_dir)
    toc_dest = si_dir / TOC_OUT
    shutil.copy2(toc_built, toc_dest)

    n_reader = PdfReader(str(n_term))
    toc_reader = PdfReader(str(toc_dest))
    text_reader = PdfReader(str(supp_text))
    pages = []
    # Title page(s)
    for i in range(min(N_TERM_TITLE_PAGES, len(n_reader.pages))):
        pages.append(n_reader.pages[i])
    # Fresh TOC, then Supplementary Text
    pages.extend(toc_reader.pages)
    pages.extend(text_reader.pages)
    # Any leftover N-term narrative after dropped stale TOC page(s)
    start = N_TERM_TITLE_PAGES + N_TERM_DROP_AFTER_TITLE
    for i in range(start, len(n_reader.pages)):
        pages.append(n_reader.pages[i])
    print(
        f"[SI PDF] N-term splice: keep pages 1-{N_TERM_TITLE_PAGES}, "
        f"insert TOC ({len(toc_reader.pages)} p), "
        f"insert {SUPPLEMENTARY_TEXT} ({len(text_reader.pages)} p), "
        f"keep N-term pages {start+1}-{len(n_reader.pages) if len(n_reader.pages) > start else 'none'}"
    )
    return pages


def build_si_merged(
    *,
    repo_root: Path,
    figures_dir: Path,
    si_dir: Path,
    tex_dir: Path,
    sf16_page_index: int = DEFAULT_SF16_PAGE_INDEX,
    skip_c_term: bool = False,
    chunk_tables: bool = False,
) -> Path:
    try:
        from pypdf import PdfReader, PdfWriter
    except ImportError as exc:
        raise SystemExit(
            "pypdf is required for SI PDF merge. Install with: pip install pypdf"
        ) from exc

    _which_latex()

    if chunk_tables:
        for stem, dest_name in TABLE_MASTERS:
            print(f"[SI PDF] Compiling {stem}.tex -> {dest_name}")
            pdf_path = _compile_tex(stem, tex_dir)
            shutil.copy2(pdf_path, si_dir / dest_name)
        table_pdfs = [si_dir / dest for _, dest in TABLE_MASTERS]
    else:
        print(f"[SI PDF] Compiling {TABLES_MAIN}.tex -> {TABLES_OUT}")
        pdf_path = _compile_tex(TABLES_MAIN, tex_dir)
        dest = si_dir / TABLES_OUT
        shutil.copy2(pdf_path, dest)
        table_pdfs = [dest]

    if not skip_c_term:
        print("[SI PDF] Building C-term figures PDF")
        _build_c_term_pdf(
            tex_dir,
            si_dir,
            figures_dir,
            sf16_page_index=sf16_page_index,
        )

    writer = PdfWriter()
    for page in _merge_n_term_with_toc(si_dir, tex_dir):
        writer.add_page(page)

    for path in table_pdfs:
        reader = PdfReader(str(path))
        print(f"[SI PDF] Merging {path.name} ({len(reader.pages)} pages)")
        for page in reader.pages:
            writer.add_page(page)

    c_term = si_dir / C_TERM_OUT
    if not c_term.is_file():
        raise FileNotFoundError(f"Missing C-term PDF: {c_term}")
    c_reader = PdfReader(str(c_term))
    print(f"[SI PDF] Merging {c_term.name} ({len(c_reader.pages)} pages)")
    for page in c_reader.pages:
        writer.add_page(page)

    merged = si_dir / MERGED
    with merged.open("wb") as f:
        writer.write(f)

    root_copy = repo_root / "SI_merged.pdf"
    shutil.copy2(merged, root_copy)
    print(f"[SI PDF] Wrote {merged} and {root_copy} ({len(writer.pages)} pages)")
    return merged


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", type=Path, default=_REPO)
    ap.add_argument("--figures-dir", type=Path, default=None)
    ap.add_argument("--si-dir", type=Path, default=None)
    ap.add_argument("--tex-dir", type=Path, default=None)
    ap.add_argument(
        "--sf16-page-index",
        type=int,
        default=DEFAULT_SF16_PAGE_INDEX,
        help="0-based page index of SF16 in a previous C-term PDF.",
    )
    ap.add_argument(
        "--skip-c-term",
        action="store_true",
        help="Reuse existing SI C term.pdf instead of rebuilding from figures.",
    )
    ap.add_argument(
        "--chunk-tables",
        action="store_true",
        help="Compile legacy per-section table masters instead of main.tex.",
    )
    # Kept for CLI compatibility; ignored (references PDF is used instead).
    ap.add_argument("--ref-pages-from-end", type=int, default=0, help=argparse.SUPPRESS)
    args = ap.parse_args()

    repo = args.repo_root.resolve()
    figures = (args.figures_dir or repo / "figures").resolve()
    si_dir = (args.si_dir or repo / "SI_documents").resolve()
    tex_dir = (args.tex_dir or si_dir / "tex").resolve()

    try:
        build_si_merged(
            repo_root=repo,
            figures_dir=figures,
            si_dir=si_dir,
            tex_dir=tex_dir,
            sf16_page_index=args.sf16_page_index,
            skip_c_term=args.skip_c_term,
            chunk_tables=args.chunk_tables,
        )
    except subprocess.CalledProcessError as exc:
        print(f"[SI PDF] LaTeX failed with exit {exc.returncode}", file=sys.stderr)
        return exc.returncode or 1
    except Exception as exc:  # noqa: BLE001
        print(f"[SI PDF] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
