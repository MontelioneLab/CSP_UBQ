#!/usr/bin/env python3
"""
Create a LaTeX table from CSP_UBQ.csv and confusion_matrix_per_system.csv.

Columns: apo_pdb, holo_pdb, apo_bmrb, holo_bmrb, f1, mcc.
Data is joined on holo_pdb = system_id. By default only rows with F1/MCC
in confusion_matrix_per_system.csv are included (inner join). Use --all to
include every CSP_UBQ row (left join), with missing values left blank.

To populate F1/MCC for more targets: ensure each target has
outputs/<holo_pdb>/master_alignment.csv (from the pipeline and/or
scripts/merge_csv.py), then run:
  python scripts/confusion_matrix_analysis.py --outputs outputs/
This overwrites outputs/confusion_matrix_per_system.csv. Then regenerate
this LaTeX table.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# SI Tables S2–S13 / S14 and SI Figs S1–S12 / S13 (see SELECTION_NUMBERING)
SELECTION_NUMBERING: dict[str, dict[str, str]] = {
    "hydrolases": {"st_id": "ST2", "sf_id": "SF1", "title": "Hydrolase Receptors", "file_token": "hydrolases"},
    "isomerases": {"st_id": "ST3", "sf_id": "SF2", "title": "Isomerase Receptors", "file_token": "isomerases"},
    "oxidoreductases": {"st_id": "ST4", "sf_id": "SF3", "title": "Oxidoreductase Receptors", "file_token": "oxidoreductases"},
    "transferases": {"st_id": "ST5", "sf_id": "SF4", "title": "Transferase Receptors", "file_token": "transferases"},
    "translocases": {"st_id": "ST6", "sf_id": "SF5", "title": "Translocase Receptors", "file_token": "translocases"},
    "all_alpha_proteins": {"st_id": "ST7", "sf_id": "SF6", "title": "All Alpha Receptors", "file_token": "all_alpha_proteins"},
    "all_beta_proteins": {"st_id": "ST8", "sf_id": "SF7", "title": "All Beta Receptors", "file_token": "all_beta_proteins"},
    "alpha_and_beta_proteins_(a+b)": {"st_id": "ST9", "sf_id": "SF8", "title": "Alpha and Beta (a+b) Receptors", "file_token": "alpha_and_beta_proteins_a_plus_b"},
    "alpha_and_beta_proteins_a_plus_b": {"st_id": "ST9", "sf_id": "SF8", "title": "Alpha and Beta (a+b) Receptors", "file_token": "alpha_and_beta_proteins_a_plus_b"},
    "cbp": {"st_id": "ST10", "sf_id": "SF9", "title": "CBP Domain Receptors", "file_token": "CBP"},
    "bet_et": {"st_id": "ST11", "sf_id": "SF10", "title": "BET-ET Domain Receptors", "file_token": "BET_ET"},
    "tfiih": {"st_id": "ST12", "sf_id": "SF11", "title": "TFIIH Domain Receptors", "file_token": "TFIIH"},
    "ubiquitin": {"st_id": "ST13", "sf_id": "SF12", "title": "Ubiquitin Domain Receptors", "file_token": "ubiquitin"},
    "dissimilar_apo_holo_conditions": {
        "st_id": "ST14",
        "sf_id": "SF13",
        "title": "targets with dissimilar apo/holo experimental conditions",
        "file_token": "dissimilar_apo_holo_conditions",
    },
}


def escape_underscores(value):
    """Escape underscores for LaTeX (replace _ with \\_)."""
    if isinstance(value, str):
        return value.replace("_", r"\_")
    return value


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Create LaTeX table from CSP and confusion matrix CSVs."
    )
    parser.add_argument(
        "--csp-csv",
        type=Path,
        default=Path("data/CSP_UBQ.csv"),
        help="Path to CSP_UBQ.csv (default: data/CSP_UBQ.csv).",
    )
    parser.add_argument(
        "--confusion-csv",
        type=Path,
        default=Path("outputs") / "confusion_matrix_per_system.csv",
        help="Path to confusion_matrix_per_system.csv (default: outputs/confusion_matrix_per_system.csv).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs") / "csp_confusion_latex_table.tex",
        help="Path to output .tex file (default: outputs/csp_confusion_latex_table.tex).",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        dest="include_all",
        help="Include all CSP_UBQ rows (left join); leave missing F1/MCC blank.",
    )
    parser.add_argument(
        "--targets-csv",
        type=Path,
        default=None,
        help="Optional CSV with a 'holo_pdb' column to filter rows.",
    )
    parser.add_argument(
        "--targets",
        type=str,
        default=None,
        help="Optional comma-separated holo_pdb IDs to filter rows.",
    )
    parser.add_argument(
        "--selection-name",
        type=str,
        default=None,
        help="Optional custom selection name to include in the LaTeX caption.",
    )
    parser.add_argument(
        "--include-figure",
        action="store_true",
        default=True,
        dest="include_figure",
        help="Include the corresponding SI figure block in the output (default: True).",
    )
    parser.add_argument(
        "--no-include-figure",
        action="store_false",
        dest="include_figure",
        help="Omit the figure block from the output.",
    )
    return parser.parse_args(argv)


def _bmrb_to_str(val) -> str:
    """Convert BMRB ID to string without decimal (e.g. 4700.0 -> '4700')."""
    if pd.isna(val):
        return ""
    s = str(val).strip()
    if s in ("", "nan", "NaN"):
        return ""
    try:
        f = float(s)
        if f == int(f):
            return str(int(f))
    except (ValueError, TypeError):
        pass
    return s


def _normalize_pdb_id(val) -> str:
    """Normalize PDB IDs and repair known spreadsheet-scientific-notation artifacts."""
    if pd.isna(val):
        return ""
    s = str(val).strip().upper()
    if s in ("", "NAN"):
        return ""
    # Known CSV artifact: 1OO9 was auto-converted to scientific notation in spreadsheet software.
    aliases = {
        "1.00E+91": "1OO9",
    }
    return aliases.get(s, s)


def _build_system_mapping(confusion_df: pd.DataFrame, outputs_root: Path) -> pd.DataFrame:
    """
    Build a mapping from system_id -> (holo_pdb, apo_bmrb, holo_bmrb, f1_score, mcc).

    For systems with suffixed IDs (e.g. 2M0G_1, 2M0G_2), this inspects
    outputs/<system_id>/master_alignment.csv to recover the apo/holo BMRB IDs
    and the underlying holo_pdb. This lets us distinguish multiple apo BMRB
    partners for the same holo PDB when building the LaTeX table.
    """
    records: list[dict] = []
    for _, row in confusion_df.iterrows():
        system_id = _normalize_pdb_id(row["system_id"])
        system_dir = outputs_root / system_id
        ma_path = system_dir / "master_alignment.csv"

        holo_pdb = ""
        apo_bmrb = ""
        holo_bmrb = ""
        if ma_path.exists():
            try:
                ma_df = pd.read_csv(ma_path, dtype=str)
                if not ma_df.empty:
                    ma_row = ma_df.iloc[0]
                    holo_pdb = _normalize_pdb_id(ma_row.get("holo_pdb", ""))
                    apo_bmrb = _bmrb_to_str(ma_row.get("apo_bmrb", ""))
                    holo_bmrb = _bmrb_to_str(ma_row.get("holo_bmrb", ""))
            except Exception:
                # Fall back to base system_id if anything goes wrong reading alignment.
                pass

        if not holo_pdb:
            # Fallback: assume base PDB code is the part before any underscore.
            base = system_id.split("_", 1)[0]
            holo_pdb = _normalize_pdb_id(base)

        records.append(
            {
                "system_id": system_id,
                "holo_pdb": holo_pdb,
                "apo_bmrb": apo_bmrb,
                "holo_bmrb": holo_bmrb,
                "f1_score": row.get("f1_score"),
                "mcc": row.get("mcc"),
            }
        )

    return pd.DataFrame.from_records(records)


def _load_allowed_targets(targets_csv: Path | None, targets_str: str | None) -> set[str] | None:
    """Load optional holo_pdb filters from CSV and/or comma-separated string."""
    allowed: set[str] = set()
    if targets_csv is not None:
        if not targets_csv.exists():
            raise FileNotFoundError(f"Targets CSV not found: {targets_csv}")
        tdf = pd.read_csv(targets_csv, dtype=str)
        if "holo_pdb" not in tdf.columns:
            raise ValueError(f"Targets CSV must include 'holo_pdb': {targets_csv}")
        allowed.update(
            _normalize_pdb_id(v)
            for v in tdf["holo_pdb"].astype(str).tolist()
            if str(v).strip()
        )
    if targets_str:
        allowed.update(
            _normalize_pdb_id(v)
            for v in targets_str.split(",")
            if str(v).strip()
        )
    return allowed or None


def load_and_merge(
    csp_path: Path,
    confusion_path: Path,
    include_all: bool = False,
    allowed_targets: set[str] | None = None,
) -> pd.DataFrame:
    """
    Load both CSVs and join CSP_UBQ rows to confusion-matrix rows.

    For holo_pdb values that have multiple apo partners (e.g. 2M0G_1, 2M0G_2),
    we disambiguate using apo_bmrb by reading outputs/<system_id>/master_alignment.csv
    and mapping on (holo_pdb, apo_bmrb).
    """
    # Keep IDs as text so PDB codes like 1E91 are not parsed as scientific notation.
    csp_df = pd.read_csv(csp_path, dtype=str)
    confusion_df = pd.read_csv(confusion_path)

    for col in ("apo_pdb", "holo_pdb", "apo_bmrb", "holo_bmrb"):
        if col in csp_df.columns:
            csp_df[col] = csp_df[col].astype(str).str.strip()
    for col in ("apo_pdb", "holo_pdb"):
        if col in csp_df.columns:
            csp_df[col] = csp_df[col].apply(_normalize_pdb_id)
    # BMRB IDs as integer-style strings (no '.0')
    for col in ("apo_bmrb", "holo_bmrb"):
        if col in csp_df.columns:
            csp_df[col] = csp_df[col].apply(_bmrb_to_str)
    if allowed_targets is not None:
        csp_df = csp_df[csp_df["holo_pdb"].isin(allowed_targets)].copy()

    confusion_df["system_id"] = confusion_df["system_id"].apply(_normalize_pdb_id)
    outputs_root = confusion_path.parent
    system_map = _build_system_mapping(confusion_df, outputs_root=outputs_root)

    # First, precise join on (holo_pdb, apo_bmrb) so multi-apo systems like 2M0G
    # map to the correct suffixed system_id row.
    merged = pd.merge(
        csp_df,
        system_map[["holo_pdb", "apo_bmrb", "f1_score", "mcc"]],
        on=["holo_pdb", "apo_bmrb"],
        how="left" if include_all else "left",
    )

    # Optional fallback: for rows still missing F1/MCC, try a simpler match on holo_pdb
    # using systems that do not carry an apo_bmrb (e.g. aggregate-only systems).
    mask_missing = merged["f1_score"].isna() & merged["mcc"].isna()
    if mask_missing.any():
        fallback_map = system_map[system_map["apo_bmrb"] == ""].drop_duplicates(
            subset=["holo_pdb"]
        )
        if not fallback_map.empty:
            fb = pd.merge(
                merged.loc[mask_missing],
                fallback_map[["holo_pdb", "f1_score", "mcc"]],
                on="holo_pdb",
                how="left",
                suffixes=("", "_fb"),
            )
            for col in ("f1_score", "mcc"):
                fb[col] = fb[col].where(~fb[f"{col}_fb"].notna(), fb[f"{col}_fb"])
                fb = fb.drop(columns=[f"{col}_fb"])
            merged.update(fb)

    if not include_all:
        merged = merged[~(merged["f1_score"].isna() & merged["mcc"].isna())]

    out = merged[["apo_pdb", "holo_pdb", "apo_bmrb", "holo_bmrb", "f1_score", "mcc"]].copy()
    out = out.rename(columns={"f1_score": "f1"})
    out["_mcc_sort"] = pd.to_numeric(out["mcc"], errors="coerce")
    out = out.sort_values(
        by=["_mcc_sort", "holo_pdb", "apo_pdb"],
        ascending=[False, True, True],
        na_position="last",
    ).drop(columns=["_mcc_sort"]).reset_index(drop=True)
    return out


def _is_blank(val) -> bool:
    """True if value is missing, empty, or the string nan/NAN."""
    if pd.isna(val):
        return True
    s = str(val).strip().upper()
    return s == "" or s == "NAN"


def format_cell(val, escape: bool, numeric: bool):
    """Format a generic cell value for LaTeX; blank if missing."""
    if _is_blank(val):
        return ""
    if numeric:
        try:
            return f"{float(val):.2f}"
        except (ValueError, TypeError):
            return str(val)
    s = str(val).strip()
    if escape:
        s = escape_underscores(s)
    return s


def format_pdb_cell(val, blank_if_empty: bool = False) -> str:
    """Format a PDB ID: uppercase, blank if empty/NAN, else \\href to RCSB."""
    if _is_blank(val):
        return ""
    pdb_id = str(val).strip().upper()
    url = f"https://www.rcsb.org/structure/{pdb_id}"
    return f"\\href{{{url}}}{{{pdb_id}}}"


def format_bmrb_cell(val, blank_if_empty: bool = False) -> str:
    """Format a BMRB ID as string (no decimal); blank if empty/NAN, else \\href to BMRB."""
    if _is_blank(val):
        return ""
    # BMRB IDs are numeric strings: avoid 4700.0 from float, use integer string
    try:
        f = float(val)
        if f == int(f):
            bmrb_id = str(int(f))
        else:
            bmrb_id = str(val).strip()
    except (ValueError, TypeError):
        bmrb_id = str(val).strip()
    url = f"https://bmrb.io/data_library/summary/index.php?bmrbId={bmrb_id}"
    return f"\\href{{{url}}}{{{bmrb_id}}}"


def _format_data_row(row: pd.Series | None) -> list[str]:
    """Return 6 formatted cells for one side of the two-column table."""
    if row is None:
        return ["", "", "", "", "", ""]
    return [
        format_pdb_cell(row["apo_pdb"], blank_if_empty=True),
        format_pdb_cell(row["holo_pdb"]),
        format_bmrb_cell(row["apo_bmrb"], blank_if_empty=True),
        format_bmrb_cell(row["holo_bmrb"]),
        format_cell(row["f1"], escape=False, numeric=True),
        format_cell(row["mcc"], escape=False, numeric=True),
    ]


def build_latex(df: pd.DataFrame, caption: str, label: str) -> str:
    """Build one longtable with paired left/right row blocks per line."""
    split_idx = (len(df) + 1) // 2
    left_df = df.iloc[:split_idx].reset_index(drop=True)
    right_df = df.iloc[split_idx:].reset_index(drop=True)
    total_lines = max(len(left_df), len(right_df))

    header = (
        "% In your document preamble: \\usepackage{booktabs}, \\usepackage{hyperref}, \\usepackage{longtable}, \\usepackage{caption}\n"
        "\\scriptsize\n"
        "\\setlength{\\tabcolsep}{3pt}\n"
        "\\begin{longtable}{llllrr@{\\hspace{0.8cm}}llllrr}\n"
        "\\captionsetup{labelformat=empty}\n"
        f"\\caption{{{caption}}}\n"
        f"\\label{{{label}}} \\\\\n"
        "\\toprule\n"
        "apo\\_pdb & holo\\_pdb & apo\\_bmrb & holo\\_bmrb & F1 & MCC & "
        "apo\\_pdb & holo\\_pdb & apo\\_bmrb & holo\\_bmrb & F1 & MCC \\\\\n"
        "\\midrule\n"
        "\\endfirsthead\n"
        "\\caption[]{(continued)} \\\\\n"
        "\\toprule\n"
        "apo\\_pdb & holo\\_pdb & apo\\_bmrb & holo\\_bmrb & F1 & MCC & "
        "apo\\_pdb & holo\\_pdb & apo\\_bmrb & holo\\_bmrb & F1 & MCC \\\\\n"
        "\\midrule\n"
        "\\endhead\n"
        "\\bottomrule\n"
        "\\endfoot\n"
        "\\bottomrule\n"
        "\\endlastfoot\n"
    )

    body_lines: list[str] = []
    for i in range(total_lines):
        left_row = left_df.iloc[i] if i < len(left_df) else None
        right_row = right_df.iloc[i] if i < len(right_df) else None
        left_cells = _format_data_row(left_row)
        right_cells = _format_data_row(right_row)
        body_lines.append(" & ".join(left_cells + right_cells) + " \\\\\n")

    return header + "".join(body_lines) + "\\end{longtable}\n"


def _sf_id_to_si_label(sf_id: str) -> str:
    """Convert SF1 -> S1, SF2 -> S2, etc."""
    if sf_id and sf_id.upper().startswith("SF"):
        return "S" + sf_id[2:]
    return sf_id or ""


def build_plot_block(
    selection_slug: str,
    selection_label: str,
    selection_file_token: str,
    sf_id: str | None = None,
    sf_title: str | None = None,
) -> str:
    """
    Build a LaTeX figure block with the single supplementary plot.

    Expects plot file named: <sf_id>_<selection_file_token>.png
    e.g. SF1_hydrolases.png
    """
    if sf_id:
        png_file = f"{sf_id}_{selection_file_token}.png"
    else:
        png_file = f"SF1_{selection_file_token}.png"
    si_label = _sf_id_to_si_label(sf_id or "SF1")
    title_part = sf_title if (sf_id and sf_title) else escape_underscores(selection_label)
    caption_text = f"\\textbf{{SI Fig. {si_label}. Confusion Matrix Histograms for {title_part}}}"
    fig_label = f"fig:{sf_id.lower()}_{selection_slug}" if sf_id else f"fig:sf1_{selection_slug}"
    return (
        "\n"
        "\\begin{figure}[htbp]\n"
        "\\centering\n"
        f"\\includegraphics[width=\\linewidth]{{{png_file}}}\n"
        "\\captionsetup{labelformat=empty}\n"
        f"\\caption{{{caption_text}}}\n"
        f"\\label{{{fig_label}}}\n"
        "\\end{figure}\n"
    )


def _normalize_selection_key(selection_name: str) -> str:
    return selection_name.strip().lower().replace("-", "_").replace(" ", "_")


def main(argv=None) -> int:
    args = parse_args(argv)
    csp_path = args.csp_csv.resolve()
    confusion_path = args.confusion_csv.resolve()
    out_path = args.output.resolve()

    if not csp_path.exists():
        print(f"Error: CSP CSV not found: {csp_path}", file=sys.stderr)
        return 1
    if not confusion_path.exists():
        print(f"Error: Confusion CSV not found: {confusion_path}", file=sys.stderr)
        return 1

    try:
        allowed_targets = _load_allowed_targets(args.targets_csv, args.targets)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    df = load_and_merge(
        csp_path,
        confusion_path,
        include_all=args.include_all,
        allowed_targets=allowed_targets,
    )
    if df.empty:
        print("No rows after merge; nothing to write.", file=sys.stderr)
        return 1

    if args.selection_name:
        selection_raw = str(args.selection_name).strip()
        selection_slug = _normalize_selection_key(selection_raw)
        numbering = SELECTION_NUMBERING.get(selection_slug)
        if numbering:
            si_num = numbering["st_id"][2:]  # ST2 -> S2, ST10 -> S10
            caption = f"\\textbf{{SI Table S{si_num}. Data for {numbering['title']}}}"
            label = f"tab:{numbering['st_id'].lower()}"
        else:
            selection_label = escape_underscores(selection_raw)
            caption = (
                f"CSP-UBQ systems ({selection_label}): "
                "PDB/BMRB identifiers and per-system F1 and MCC."
            )
            label = f"tab:csp_confusion_{selection_slug}"
    else:
        caption = "Data for all receptors in this study"
        label = "tab:st13"
    latex = build_latex(df, caption=caption, label=label)

    # Append the corresponding SI figure block when we have a selection with numbering
    if args.include_figure and args.selection_name:
        selection_slug = _normalize_selection_key(str(args.selection_name).strip())
        numbering = SELECTION_NUMBERING.get(selection_slug)
        if numbering:
            plot_block = build_plot_block(
                selection_slug=selection_slug,
                selection_label=numbering["title"],
                selection_file_token=numbering["file_token"],
                sf_id=numbering["sf_id"],
                sf_title=numbering["title"],
            )
            latex = latex + plot_block

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(latex, encoding="utf-8")
    print(f"Wrote {len(df)} rows to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
