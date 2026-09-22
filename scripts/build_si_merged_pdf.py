#!/usr/bin/env python3
"""
Compile SI table/figure TeX masters and merge into SI_documents/SI merged.pdf.

With --combined, write Supplemental_Information.pdf instead of SI_merged.pdf:
the same SI tables/figures/references with a combined TOC and SI Text, plus
All_Case_Studies.pdf after the References.

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
TOC_COMBINED_MASTER = "si_toc_combined"
TOC_COMBINED_OUT = "SI_TOC_combined.pdf"
TITLE_MASTER = "si_title"
TITLE_OUT = "SI_title.pdf"
SUPP_TEXT_MASTER = "si_supplementary_text"
SUPP_TEXT_COMBINED_MASTER = "si_supplementary_text_combined"
SUPP_TEXT_COMBINED_OUT = "Supplementary Text combined.pdf"
CASE_STUDIES_PDF = "All_Case_Studies.pdf"
COMBINED_OUT = "Supplemental_Information.pdf"

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
REFERENCES_MASTER = "si_references"
MERGED = "SI merged.pdf"

# Heuristic SF19 (PDB Advanced Search) page index (0-based) in a prior C-term PDF when PNG is missing.
# C-term pages are S9–S18 then PDB as S19 (insert index 10).
DEFAULT_SF19_PAGE_INDEX = 10
DEFAULT_SF19_INSERT_INDEX = 10
# Compiled title page replaces N-term page 0. N-term page 1 = stale TOC (dropped);
# pages 2+ = leftover narrative if present.
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


def _ensure_references_pdf(si_dir: Path, c_term_pdf: Path, tex_dir: Path | None = None) -> Path:
    """Compile si_references.tex when present; else reuse or extract the static PDF."""
    from pypdf import PdfReader, PdfWriter

    refs = si_dir / REFERENCES
    tex_root = tex_dir if tex_dir is not None else si_dir / "tex"
    refs_master = tex_root / f"{REFERENCES_MASTER}.tex"
    if refs_master.is_file():
        print(f"[SI PDF] Compiling {REFERENCES_MASTER}.tex -> {REFERENCES}")
        built = _compile_tex(REFERENCES_MASTER, tex_root)
        shutil.copy2(built, refs)
        return refs
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


def _extract_sf19_png_from_c_term(
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
    page_text = (reader.pages[page_index].extract_text() or "").lower()
    if "advanced search" not in page_text and "pdb advanced" not in page_text:
        print(
            "[SI PDF] SF19 PNG unavailable; prior C-term page "
            f"{page_index} is not a PDB Advanced Search page — not extracting."
        )
        return False
    writer = PdfWriter()
    writer.add_page(reader.pages[page_index])
    tmp_pdf = out_png.with_suffix(".pdf")
    with tmp_pdf.open("wb") as f:
        writer.write(f)
    print(f"[SI PDF] SF19 PNG unavailable; saved page PDF to {tmp_pdf}")
    return False


def _pdf_looks_like_pdb_search(path: Path) -> bool:
    """True if a fallback SF19 PDF is actually the PDB Advanced Search page."""
    try:
        from pypdf import PdfReader
    except ImportError:
        return False
    try:
        text = (PdfReader(str(path)).pages[0].extract_text() or "").lower()
    except Exception:
        return False
    if "advanced search" in text or "pdb advanced" in text:
        return True
    if "buffer" in text and ("threshold" in text or "sweep" in text):
        return False
    if "offset grid" in text:
        return False
    return False


def _footer_overlay_page(width: float, height: float, number: int):
    """Blank overlay with a white footer patch and a centered page number."""
    from pypdf import PageObject
    from pypdf.generic import DecodedStreamObject, DictionaryObject, NameObject

    overlay = PageObject.create_blank_page(width=width, height=height)
    font = DictionaryObject()
    font[NameObject("/Type")] = NameObject("/Font")
    font[NameObject("/Subtype")] = NameObject("/Type1")
    font[NameObject("/BaseFont")] = NameObject("/Helvetica")
    font[NameObject("/Encoding")] = NameObject("/WinAnsiEncoding")
    fonts = DictionaryObject()
    fonts[NameObject("/F1")] = font
    resources = DictionaryObject()
    resources[NameObject("/Font")] = fonts
    overlay[NameObject("/Resources")] = resources

    label = str(number)
    font_size = 10.0
    text_w = 0.556 * font_size * len(label)
    y = 22.0
    x_text = width / 2.0 - text_w / 2.0
    rect_w = max(52.0, text_w + 18.0)
    rect_h = 16.0
    x_rect = width / 2.0 - rect_w / 2.0
    y_rect = y - 4.0
    content = (
        f"q 1 1 1 rg {x_rect:.2f} {y_rect:.2f} {rect_w:.2f} {rect_h:.2f} re f Q\n"
        f"BT /F1 {font_size:.1f} Tf {x_text:.2f} {y:.2f} Td ({label}) Tj ET\n"
    )
    stream = DecodedStreamObject()
    stream.set_data(content.encode("latin-1"))
    overlay.replace_contents(stream)
    return overlay


def _stamp_continuous_page_numbers(writer) -> None:
    """Number pages 1..N in the footer and set matching PDF page labels."""
    from pypdf.constants import PageLabelStyle

    n_pages = len(writer.pages)
    for i, page in enumerate(writer.pages, start=1):
        box = page.mediabox
        overlay = _footer_overlay_page(float(box.width), float(box.height), i)
        page.merge_page(overlay)
    if n_pages:
        writer.set_page_label(0, n_pages - 1, style=PageLabelStyle.DECIMAL, start=1)
    print(f"[SI PDF] Stamped continuous page numbers 1-{n_pages}")


def _build_c_term_pdf(
    tex_dir: Path,
    si_dir: Path,
    figures_dir: Path,
    *,
    sf19_page_index: int,
) -> Path:
    from pypdf import PdfReader, PdfWriter

    existing_c = si_dir / C_TERM_OUT
    refs_pdf = _ensure_references_pdf(si_dir, existing_c, tex_dir)

    sf19_png = figures_dir / "SF19_pdb_search.png"
    sf19_pdf = figures_dir / "SF19_pdb_search.pdf"
    if not sf19_png.is_file() and not sf19_pdf.is_file() and existing_c.is_file():
        _extract_sf19_png_from_c_term(existing_c, sf19_png, sf19_page_index)

    fig_pdf = _compile_tex(C_TERM_MASTER, tex_dir)
    writer = PdfWriter()
    for page in PdfReader(str(fig_pdf)).pages:
        writer.add_page(page)

    if not sf19_png.is_file() and sf19_pdf.is_file() and _pdf_looks_like_pdb_search(sf19_pdf):
        insert_at = min(DEFAULT_SF19_INSERT_INDEX, len(writer.pages))
        sf19_pages = PdfReader(str(sf19_pdf)).pages
        pages = list(writer.pages)
        new_writer = PdfWriter()
        for i, page in enumerate(pages):
            if i == insert_at:
                for sp in sf19_pages:
                    new_writer.add_page(sp)
            new_writer.add_page(page)
        if insert_at >= len(pages):
            for sp in sf19_pages:
                new_writer.add_page(sp)
        writer = new_writer
    elif not sf19_png.is_file():
        print(
            "[SI PDF] SF19 PDB-search PNG/PDF missing or not a PDB Advanced Search page; "
            "skipping S19 splice (TOC still lists it)."
        )

    # Append dedicated references only (never re-append old SF19).
    for page in PdfReader(str(refs_pdf)).pages:
        writer.add_page(page)
    print(f"[SI PDF] Appended references from {refs_pdf.name}")

    out = si_dir / C_TERM_OUT
    with out.open("wb") as f:
        writer.write(f)
    return out


def _merge_n_term_with_toc(
    si_dir: Path,
    tex_dir: Path,
    *,
    toc_master: str = TOC_MASTER,
    toc_out_name: str = TOC_OUT,
    supp_text: Path | None = None,
) -> list:
    """Return pages: compiled title + TOC + Supplementary Text + leftover N-term narrative."""
    from pypdf import PdfReader

    print(f"[SI PDF] Compiling {TITLE_MASTER}.tex -> {TITLE_OUT}")
    title_built = _compile_tex(TITLE_MASTER, tex_dir)
    title_dest = si_dir / TITLE_OUT
    shutil.copy2(title_built, title_dest)

    if supp_text is None:
        print(f"[SI PDF] Compiling {SUPP_TEXT_MASTER}.tex -> {SUPPLEMENTARY_TEXT}")
        text_built = _compile_tex(SUPP_TEXT_MASTER, tex_dir)
        supp_text = si_dir / SUPPLEMENTARY_TEXT
        shutil.copy2(text_built, supp_text)
    if not supp_text.is_file():
        raise FileNotFoundError(f"Missing Supplementary Text PDF: {supp_text}")

    print(f"[SI PDF] Compiling {toc_master}.tex -> {toc_out_name}")
    toc_built = _compile_tex(toc_master, tex_dir)
    toc_dest = si_dir / toc_out_name
    shutil.copy2(toc_built, toc_dest)

    n_term = si_dir / N_TERM
    toc_reader = PdfReader(str(toc_dest))
    text_reader = PdfReader(str(supp_text))
    title_reader = PdfReader(str(title_dest))
    pages = list(title_reader.pages)
    pages.extend(toc_reader.pages)
    pages.extend(text_reader.pages)
    leftover_end = "none"
    if n_term.is_file():
        n_reader = PdfReader(str(n_term))
        start = N_TERM_TITLE_PAGES + N_TERM_DROP_AFTER_TITLE
        for i in range(start, len(n_reader.pages)):
            pages.append(n_reader.pages[i])
        leftover_end = (
            f"{start + 1}-{len(n_reader.pages)}" if len(n_reader.pages) > start else "none"
        )
    print(
        f"[SI PDF] Title splice: compiled {TITLE_OUT} ({len(title_reader.pages)} p), "
        f"insert TOC ({len(toc_reader.pages)} p), "
        f"insert {supp_text.name} ({len(text_reader.pages)} p), "
        f"keep N-term leftover pages {leftover_end}"
    )
    return pages


def build_si_merged(
    *,
    repo_root: Path,
    figures_dir: Path,
    si_dir: Path,
    tex_dir: Path,
    sf19_page_index: int = DEFAULT_SF19_PAGE_INDEX,
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
            sf19_page_index=sf19_page_index,
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

    _stamp_continuous_page_numbers(writer)

    merged = si_dir / MERGED
    with merged.open("wb") as f:
        writer.write(f)

    root_copy = repo_root / "SI_merged.pdf"
    shutil.copy2(merged, root_copy)
    print(f"[SI PDF] Wrote {merged} and {root_copy} ({len(writer.pages)} pages)")
    return merged


def build_supplemental_information(
    *,
    repo_root: Path,
    si_dir: Path,
    tex_dir: Path,
    case_studies_pdf: Path | None = None,
) -> Path:
    """Merge variant TOC/SI text + existing SI parts + All_Case_Studies.pdf.

    Does not rewrite SI_merged.pdf. Reuses existing tables and C-term PDFs.
    """
    try:
        from pypdf import PdfReader, PdfWriter
    except ImportError as exc:
        raise SystemExit(
            "pypdf is required for SI PDF merge. Install with: pip install pypdf"
        ) from exc

    _which_latex()

    tables = si_dir / TABLES_OUT
    c_term = si_dir / C_TERM_OUT
    if not tables.is_file():
        raise FileNotFoundError(
            f"Missing tables PDF: {tables}. Run without --combined first."
        )
    if not c_term.is_file():
        raise FileNotFoundError(
            f"Missing C-term PDF: {c_term}. Run without --combined first."
        )

    case_pdf = (case_studies_pdf or repo_root / CASE_STUDIES_PDF).resolve()
    if not case_pdf.is_file():
        raise FileNotFoundError(f"Missing case-studies PDF: {case_pdf}")

    print(
        f"[SI PDF] Compiling {SUPP_TEXT_COMBINED_MASTER}.tex -> {SUPP_TEXT_COMBINED_OUT}"
    )
    text_built = _compile_tex(SUPP_TEXT_COMBINED_MASTER, tex_dir)
    text_dest = si_dir / SUPP_TEXT_COMBINED_OUT
    shutil.copy2(text_built, text_dest)

    writer = PdfWriter()
    for page in _merge_n_term_with_toc(
        si_dir,
        tex_dir,
        toc_master=TOC_COMBINED_MASTER,
        toc_out_name=TOC_COMBINED_OUT,
        supp_text=text_dest,
    ):
        writer.add_page(page)

    tables_reader = PdfReader(str(tables))
    print(f"[SI PDF] Merging {tables.name} ({len(tables_reader.pages)} pages)")
    for page in tables_reader.pages:
        writer.add_page(page)

    c_reader = PdfReader(str(c_term))
    print(f"[SI PDF] Merging {c_term.name} ({len(c_reader.pages)} pages)")
    for page in c_reader.pages:
        writer.add_page(page)

    case_reader = PdfReader(str(case_pdf))
    print(f"[SI PDF] Appending {case_pdf.name} ({len(case_reader.pages)} pages)")
    for page in case_reader.pages:
        writer.add_page(page)

    _stamp_continuous_page_numbers(writer)

    merged = si_dir / COMBINED_OUT
    with merged.open("wb") as f:
        writer.write(f)

    root_copy = repo_root / COMBINED_OUT
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
        "--sf19-page-index",
        type=int,
        default=DEFAULT_SF19_PAGE_INDEX,
        help="0-based page index of SF19 (PDB Advanced Search) in a previous C-term PDF.",
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
    ap.add_argument(
        "--combined",
        action="store_true",
        help=(
            "Build Supplemental_Information.pdf from existing SI parts plus "
            "All_Case_Studies.pdf. Does not rewrite SI_merged.pdf."
        ),
    )
    ap.add_argument(
        "--case-studies-pdf",
        type=Path,
        default=None,
        help="Path to All_Case_Studies.pdf (default: <repo-root>/All_Case_Studies.pdf).",
    )
    # Kept for CLI compatibility; ignored (references PDF is used instead).
    ap.add_argument("--ref-pages-from-end", type=int, default=0, help=argparse.SUPPRESS)
    args = ap.parse_args()

    repo = args.repo_root.resolve()
    figures = (args.figures_dir or repo / "figures").resolve()
    si_dir = (args.si_dir or repo / "SI_documents").resolve()
    tex_dir = (args.tex_dir or si_dir / "tex").resolve()

    try:
        if args.combined:
            build_supplemental_information(
                repo_root=repo,
                si_dir=si_dir,
                tex_dir=tex_dir,
                case_studies_pdf=args.case_studies_pdf,
            )
        else:
            build_si_merged(
                repo_root=repo,
                figures_dir=figures,
                si_dir=si_dir,
                tex_dir=tex_dir,
                sf19_page_index=args.sf19_page_index,
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
