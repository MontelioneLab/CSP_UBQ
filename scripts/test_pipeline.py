"""
Tests for every step of the CSP pipeline.

Requires: pytest and all pipeline dependencies (numpy, mdtraj, matplotlib, requests, etc.)

Run: pytest scripts/test_pipeline.py
Run integration tests: pytest scripts/test_pipeline.py -m integration
"""

from __future__ import annotations

import csv
import os
import sys
import tempfile
from pathlib import Path
from typing import Optional
from unittest.mock import MagicMock, patch

# Ensure project root is on path for both pytest and direct execution
_project_root = Path(__file__).resolve().parent.parent
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

import pytest

# Support running as script or module
try:
    from .config import paths, ensure_directories
    from .pipeline import (
        lookup_all_rows_from_holo_pdb,
        get_duplicate_index_for_row,
        process_row,
    )
    from .align import align_global
    from .bmrb_io import fetch_bmrb, parse_sequence_and_shifts_from_saveframes
    from .rcsb_io import fetch_pdb, parse_pdb_sequences
    from .csp import (
        CSPResult,
        compute_csp_A,
        compute_csp_multiple_saveframes,
        compute_csp_multiple_saveframes_ca,
        compute_threshold_with_outlier_removal,
        percentile_linear,
        _apply_hn_significance_masks,
        _compute_basic_stats,
        _floor_primary_hn_cutoff,
    )
    from .sasa_analysis import compute_sasa_occlusion, write_occlusion_analysis_csv
    from .interaction_analysis import compute_interaction_filter, write_interaction_analysis_csv
    from .interaction_analysis import compute_ca_distance_filter, write_ca_distance_csv
    from .interaction_analysis import (
        compute_second_shell_filter,
        write_second_shell_csv,
        second_shell_column_name,
    )
    from .merge_csv import (
        merge_all_csv_files,
        compute_classification,
        parse_optional_bool,
        has_recorded_csp,
        filter_recorded_csp_dataframe,
        is_binding_site_row,
    )
    from .analyze_targets_single_atom_shifts import compute_1d_metrics_for_target
    from .analyze_targets import load_alignment, compute_f1_score
    from .create_fig_3 import collect_distance_categories
    from .create_si_fig_s13 import collect_distance_categories as collect_s13_distance_categories
    from .create_si_fig_s14 import collect_distance_categories as collect_s14_distance_categories
    from .HSQC_visualize import (
        plot_hsqc_variants,
        resolve_hsqc_offsets,
        format_offset_annotation,
        _classify_residue,
        _classification_color_map,
        _EXCLUDED_LINE_COLOR,
    )
    from .visualize import (
        plot_csp_classification_bars,
        write_pymol_color_csp_mask_script,
        write_pymol_occlusion_script,
        write_pymol_delta_sasa_script,
    )
except Exception:
    from scripts.config import paths, ensure_directories
    from scripts.pipeline import (
        lookup_all_rows_from_holo_pdb,
        get_duplicate_index_for_row,
        process_row,
    )
    from scripts.align import align_global
    from scripts.bmrb_io import fetch_bmrb, parse_sequence_and_shifts_from_saveframes
    from scripts.rcsb_io import fetch_pdb, parse_pdb_sequences
    from scripts.csp import (
        CSPResult,
        compute_csp_A,
        compute_csp_multiple_saveframes,
        compute_csp_multiple_saveframes_ca,
        compute_threshold_with_outlier_removal,
        percentile_linear,
        _apply_hn_significance_masks,
        _compute_basic_stats,
        _floor_primary_hn_cutoff,
    )
    from scripts.sasa_analysis import compute_sasa_occlusion, write_occlusion_analysis_csv
    from scripts.interaction_analysis import compute_interaction_filter, write_interaction_analysis_csv
    from scripts.interaction_analysis import compute_ca_distance_filter, write_ca_distance_csv
    from scripts.interaction_analysis import (
        compute_second_shell_filter,
        write_second_shell_csv,
        second_shell_column_name,
    )
    from scripts.merge_csv import (
        merge_all_csv_files,
        compute_classification,
        parse_optional_bool,
        has_recorded_csp,
        filter_recorded_csp_dataframe,
        is_binding_site_row,
    )
    from scripts.analyze_targets_single_atom_shifts import compute_1d_metrics_for_target
    from scripts.analyze_targets import load_alignment, compute_f1_score
    from scripts.create_fig_3 import collect_distance_categories
    from scripts.create_si_fig_s13 import collect_distance_categories as collect_s13_distance_categories
    from scripts.create_si_fig_s14 import collect_distance_categories as collect_s14_distance_categories
    from scripts.HSQC_visualize import (
        plot_hsqc_variants,
        resolve_hsqc_offsets,
        format_offset_annotation,
        _classify_residue,
        _classification_color_map,
        _EXCLUDED_LINE_COLOR,
    )
    from scripts.visualize import (
        plot_csp_classification_bars,
        write_pymol_color_csp_mask_script,
        write_pymol_occlusion_script,
        write_pymol_delta_sasa_script,
    )


# --- Fixtures ---

@pytest.fixture
def temp_dir():
    """Temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as d:
        yield d


@pytest.fixture
def sample_csv_path(temp_dir):
    """Create a sample CSP_UBQ.csv for lookup tests."""
    csv_path = os.path.join(temp_dir, "data/CSP_UBQ.csv")
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["apo_bmrb", "holo_bmrb", "apo_pdb", "holo_pdb"])
        w.writerow(["18251", "4700", "", "1cf4"])
        w.writerow(["34688", "4516", "7QCX", "1d5g"])
        w.writerow(["6268", "5327", "1U2N", "1l8c"])
        w.writerow(["6268", "5327", "1U2N", "1r8u"])  # first of duplicate holo_pdb 1r8u
        w.writerow(["6268", "5987", "1U2N", "1r8u"])  # second of duplicate holo_pdb 1r8u
        w.writerow(["51289", "5480", "", "1sy9"])
    return csv_path


@pytest.fixture
def minimal_pdb_content():
    """Minimal PDB content with ATOM records for a few residues."""
    return """HEADER    TEST
ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  MET A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   MET A   1       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   MET A   1       3.000   0.000   0.000  1.00  0.00           O
ATOM      5  N   GLY A   2       4.000   0.000   0.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       5.000   0.000   0.000  1.00  0.00           C
ATOM      7  C   GLY A   2       6.000   0.000   0.000  1.00  0.00           C
ATOM      8  O   GLY A   2       7.000   0.000   0.000  1.00  0.00           O
ATOM      9  N   ALA B   1       8.000   0.000   0.000  1.00  0.00           N
ATOM     10  CA  ALA B   1       9.000   0.000   0.000  1.00  0.00           C
END
"""


@pytest.fixture
def minimal_nmr_star_v3():
    """Minimal NMR-STAR v3 format for BMRB parse test (includes α-H)."""
    return """save_assigned_chemical_shifts_1
_Assigned_chem_shift_list.Sf_category    assigned_chemical_shifts
loop_
_Atom_chem_shift.Seq_ID
_Atom_chem_shift.Comp_ID
_Atom_chem_shift.Atom_ID
_Atom_chem_shift.Val
1 MET N 120.5
1 MET H 8.41
1 MET CA 55.0
1 MET HA 4.52
2 GLY N 108.0
2 GLY H 8.2
2 GLY CA 44.0
2 GLY HA2 3.90
2 GLY HA3 3.95
stop_
save_
"""


@pytest.fixture
def minimal_nmr_star_v21():
    """Minimal NMR-STAR v2.1-style saveframe with α-H."""
    return """save_assigned_chem_shift_list_1
_Saveframe_category   assigned_chemical_shifts
loop_
_Residue_seq_code
_Residue_label
_Atom_name
_Chem_shift_value
1 MET N 120.5
1 MET H 8.41
1 MET CA 55.0
1 MET HA 4.52
stop_
save_
"""


# --- 1. Pure Unit Tests ---

class TestLookupAllRowsFromHoloPdb:
    """Tests for lookup_all_rows_from_holo_pdb."""

    def test_lookup_single_match(self, sample_csv_path):
        rows = lookup_all_rows_from_holo_pdb(sample_csv_path, "1cf4")
        assert len(rows) == 1
        assert rows[0]["apo_bmrb"] == "18251"
        assert rows[0]["holo_bmrb"] == "4700"
        assert rows[0]["holo_pdb"] == "1cf4"

    def test_lookup_duplicate_holo_pdb(self, sample_csv_path):
        rows = lookup_all_rows_from_holo_pdb(sample_csv_path, "1r8u")
        assert len(rows) == 2
        assert rows[0]["apo_bmrb"] == "6268" and rows[0]["holo_bmrb"] == "5327"
        assert rows[1]["apo_bmrb"] == "6268" and rows[1]["holo_bmrb"] == "5987"

    def test_lookup_not_found_raises(self, sample_csv_path):
        with pytest.raises(ValueError, match="No entry found for holo_pdb ID: X999"):
            lookup_all_rows_from_holo_pdb(sample_csv_path, "X999")


class TestGetDuplicateIndexForRow:
    """Tests for get_duplicate_index_for_row."""

    def test_no_duplicates_returns_none(self, sample_csv_path):
        row = {"apo_bmrb": "18251", "holo_bmrb": "4700", "holo_pdb": "1cf4"}
        assert get_duplicate_index_for_row(sample_csv_path, row) is None

    def test_duplicates_returns_index(self, sample_csv_path):
        row1 = {"apo_bmrb": "6268", "holo_bmrb": "5327", "holo_pdb": "1r8u"}
        row2 = {"apo_bmrb": "6268", "holo_bmrb": "5987", "holo_pdb": "1r8u"}
        assert get_duplicate_index_for_row(sample_csv_path, row1) == 1
        assert get_duplicate_index_for_row(sample_csv_path, row2) == 2


class TestAlignGlobal:
    """Tests for align_global."""

    def test_identical_sequences(self):
        a, b, mapping, score = align_global("ACGT", "ACGT")
        assert a == "ACGT"
        assert b == "ACGT"
        assert len(mapping) == 4
        assert score > 0

    def test_sequences_with_gaps(self):
        a, b, mapping, score = align_global("ACG", "ACGT")
        assert len(mapping) >= 3
        assert score > 0

    def test_empty_sequence(self):
        a, b, mapping, score = align_global("", "ACGT")
        assert a == "----" or a == ""
        assert len(mapping) == 0


class TestEnsureDirectories:
    """Tests for ensure_directories."""

    def test_creates_directories(self, temp_dir):
        out_dir = os.path.join(temp_dir, "test_outputs")
        ensure_directories(out_dir)
        assert os.path.isdir(out_dir)
        assert os.path.isdir(paths.cs_cache_dir)
        assert os.path.isdir(paths.pdb_cache_dir)


# --- 2. File-Based Tests ---

class TestParsePdbSequences:
    """Tests for parse_pdb_sequences."""

    def test_parse_minimal_pdb(self, temp_dir, minimal_pdb_content):
        pdb_path = os.path.join(temp_dir, "test.pdb")
        with open(pdb_path, "w") as f:
            f.write(minimal_pdb_content)
        chains = parse_pdb_sequences(pdb_path)
        assert "A" in chains
        assert "B" in chains
        assert "M" in chains["A"]  # MET
        assert "G" in chains["A"]  # GLY
        assert "A" in chains["B"]  # ALA


class TestParseSequenceAndShiftsFromSaveframes:
    """Tests for parse_sequence_and_shifts_from_saveframes."""

    def test_parse_minimal_nmr_star(self, temp_dir, minimal_nmr_star_v3):
        star_path = os.path.join(temp_dir, "test_3.str")
        with open(star_path, "w") as f:
            f.write(minimal_nmr_star_v3)
        result = parse_sequence_and_shifts_from_saveframes(star_path)
        assert isinstance(result, list)
        assert len(result) == 1
        seq, H, N, CA, HA, _sid_min, _name = result[0]
        assert seq == "MG"
        assert HA[1] == pytest.approx(4.52)
        assert HA[2] == pytest.approx(3.925)
        csv_path = os.path.join(temp_dir, "parsed", "test_3_residue_shifts.csv")
        assert os.path.isfile(csv_path)
        with open(csv_path, "r") as f:
            lines = f.read().strip().splitlines()
        assert lines[0] == "seq_id,aa,H,N,CA,HA"
        assert "4.520000" in lines[1]
        assert "3.925000" in lines[2]

    def test_parse_minimal_nmr_star_v21(self, temp_dir, minimal_nmr_star_v21):
        star_path = os.path.join(temp_dir, "case_21.str")
        with open(star_path, "w") as f:
            f.write(minimal_nmr_star_v21)
        result = parse_sequence_and_shifts_from_saveframes(star_path)
        assert len(result) == 1
        seq, _H, _N, _CA, HA, _sid_min, _name = result[0]
        assert seq == "M"
        assert HA[1] == pytest.approx(4.52)
        per_sf_csv = os.path.join(
            temp_dir,
            "parsed",
            "case_21_assigned-chem-shift-list-1_residue_shifts.csv",
        )
        assert os.path.isfile(per_sf_csv)

    def test_parse_multi_entity_v3_splits_chains(self, temp_dir):
        """Multi-entity v3 lists must not interleave overlapping Seq_IDs into a chimera."""
        fixture = """save_assigned_chemical_shifts_1
_Assigned_chem_shift_list.Sf_category    assigned_chemical_shifts
loop_
_Atom_chem_shift.Entity_ID
_Atom_chem_shift.Seq_ID
_Atom_chem_shift.Comp_ID
_Atom_chem_shift.Atom_ID
_Atom_chem_shift.Auth_asym_ID
_Atom_chem_shift.Val
1 1 MET N A 120.5
1 1 MET H A 8.41
1 1 MET CA A 55.0
1 1 MET HA A 4.52
1 2 ALA N A 122.0
1 2 ALA H A 8.20
1 2 ALA CA A 52.0
1 2 ALA HA A 4.30
2 1 GLY N B 110.0
2 1 GLY H B 8.50
2 1 GLY CA B 45.5
2 2 SER N B 116.0
2 2 SER H B 8.40
2 2 SER CA B 58.0
stop_
save_
"""
        star_path = os.path.join(temp_dir, "multi_entity_3.str")
        with open(star_path, "w") as f:
            f.write(fixture)
        result = parse_sequence_and_shifts_from_saveframes(star_path)
        assert len(result) == 2
        by_name = {name: seq for seq, _H, _N, _CA, _HA, _sid, name in result}
        assert any("chain_A" in n for n in by_name)
        assert any("chain_B" in n for n in by_name)
        seq_a = next(seq for seq, *_rest, name in result if "chain_A" in name)
        seq_b = next(seq for seq, *_rest, name in result if "chain_B" in name)
        assert seq_a == "MA"
        assert seq_b == "GS"
        # Chimera from interleaved Seq_IDs would be MGAS / MGSA-like, not MA + GS
        assert not any(s.startswith("MG") for s in by_name.values())

    def test_seqid_spaced_sequence_preserves_missing_residues(self):
        """Missing Seq_IDs must become X placeholders, not collapse the register."""
        from scripts.bmrb_io import _sequence_and_shifts_from_residue_map

        residue_to_atoms = {
            (1, "MET"): {"H": [8.1], "N": [120.0]},
            (2, "ALA"): {"H": [8.2], "N": [121.0]},
            (5, "SER"): {"H": [8.5], "N": [115.0]},
        }
        seq, H, N, _CA, _HA, seq_id_min = _sequence_and_shifts_from_residue_map(
            residue_to_atoms
        )
        assert seq == "MAXXS"
        assert seq_id_min == 1
        assert len(seq) == 5
        assert H[1] == pytest.approx(8.1)
        assert H[2] == pytest.approx(8.2)
        assert 3 not in H and 4 not in H
        assert H[5] == pytest.approx(8.5)
        assert N[5] == pytest.approx(115.0)

    def test_seqid_spaced_polymer_fill_does_not_desync_shifts(self):
        """Polymer fills gaps but shifts stay keyed by Seq_ID, not densified CS order."""
        from scripts.bmrb_io import _sequence_and_shifts_from_residue_map

        residue_to_atoms = {
            (1, "MET"): {"H": [8.1], "N": [120.0]},
            (2, "ALA"): {"H": [8.2], "N": [121.0]},
            (5, "SER"): {"H": [8.5], "N": [115.0]},
        }
        polymer = "MAGLS"  # positions 3–4 = G,L fill missing Seq_IDs
        seq, H, N, _CA, _HA, seq_id_min = _sequence_and_shifts_from_residue_map(
            residue_to_atoms, polymer_sequence=polymer
        )
        assert seq == "MAGLS"
        assert seq_id_min == 1
        # SER shift must sit at Seq_ID 5, not densified pos 3
        assert H[5] == pytest.approx(8.5)
        assert 3 not in H
        assert N[5] == pytest.approx(115.0)

    def test_align_by_seqid_offset_exact_matches_only(self):
        """Seq_ID offset alignment pairs only identical residues (7SFT-like N-term tag)."""
        from scripts.align import align_by_seqid_offset

        # Apo starts at Seq_ID 2 with a tag; shared domain begins at apo 23 / holo 1
        apo = "LXXSATTXXMXXASGSGASGSGGAHKVRAGG"
        holo = "GGAHKVRAGGP"
        aligned_a, aligned_h, mapping, score, offset = align_by_seqid_offset(
            apo, holo, apo_seq_id_min=2, holo_seq_id_min=1
        )
        assert offset == 22
        # mapping stores (apo_sid, holo_sid)
        assert mapping[0] == (23, 1)
        assert all(
            apo[asid - 2] == holo[hsid - 1] and apo[asid - 2] not in {"X", "-"}
            for asid, hsid in mapping
        )
        # Exact matched columns only in mapping — no X↔real pairs
        assert (33, 11) not in mapping  # apo X vs holo P
        assert score == 2 * len(mapping)

    def test_align_by_seqid_offset_co_columns_xx(self):
        """Consecutive shared X fillers stay co-column (no X/- zig-zag)."""
        from scripts.align import align_by_seqid_offset

        # Shared real AAs frame a run of X/X fillers (2BN5-like).
        apo = "AKQ" + ("X" * 8) + "GAD"
        holo = "AKQ" + ("X" * 8) + "GAD"
        aligned_a, aligned_h, mapping, score, offset = align_by_seqid_offset(
            apo, holo, apo_seq_id_min=1, holo_seq_id_min=1
        )
        assert offset == 0
        assert "X-" not in aligned_a and "-X" not in aligned_a
        assert "X-" not in aligned_h and "-X" not in aligned_h
        assert "XXXXXXXX" in aligned_a
        assert "XXXXXXXX" in aligned_h
        xx_cols = sum(1 for a, h in zip(aligned_a, aligned_h) if a == "X" and h == "X")
        assert xx_cols == 8
        assert all(apo[asid - 1] != "X" and holo[hsid - 1] != "X" for asid, hsid in mapping)
        assert score == 2 * len(mapping)

    @pytest.mark.integration
    def test_parse_bmrb_34638_splits_receptor_and_peptide(self):
        """BMRB 34638 (7ovc holo): receptor chain A must not merge with peptide chain B."""
        star_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "CS_Lists",
            "34638_3.str",
        )
        if not os.path.isfile(star_path):
            pytest.skip("CS_Lists/34638_3.str not present")
        result = parse_sequence_and_shifts_from_saveframes(star_path)
        assert len(result) >= 2
        seq_a = next(
            (seq for seq, _H, _N, _CA, _HA, name in result if "chain_A" in name),
            None,
        )
        assert seq_a is not None, f"expected chain_A in {[n for *_, n in result]}"
        assert seq_a.startswith("MADEATRR")
        assert not seq_a.startswith("MGADSEV")
        assert 160 <= len(seq_a) <= 170
        assert len(seq_a) != 190


class TestMergeAllCsvFiles:
    """Tests for merge_all_csv_files."""

    def test_merge_creates_master_csv(self, temp_dir):
        # Create minimal csp_table.csv
        csp_path = os.path.join(temp_dir, "csp_table.csv")
        with open(csp_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["apo_bmrb", "holo_bmrb", "holo_pdb", "chain", "apo_resi", "apo_aa", "holo_resi", "holo_aa",
                        "H_apo", "N_apo", "H_holo", "N_holo", "dH", "dN", "csp_A", "significant"])
            w.writerow(["18251", "4700", "1cf4", "A", "1", "M", "1", "M", "8.41", "122.0", "8.26", "115.7",
                        "0.04", "-1.5", "0.31", "1"])

        # Create occlusion_analysis.csv with residue_number, residue_name
        occ_path = os.path.join(temp_dir, "occlusion_analysis.csv")
        with open(occ_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["residue_number", "residue_name", "chain_id", "delta_sasa", "is_occluded"])
            w.writerow(["1", "MET", "A", "0.02", "True"])

        # Create interaction_filter.csv
        int_path = os.path.join(temp_dir, "interaction_filter.csv")
        with open(int_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["residue_number", "residue_name", "chain_id", "has_hbond"])
            w.writerow(["1", "MET", "A", "False"])

        # Create ca_distance_filter.csv
        ca_path = os.path.join(temp_dir, "ca_distance_filter.csv")
        with open(ca_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["residue_number", "residue_name", "chain_id", "passes_filter"])
            w.writerow(["1", "MET", "A", "False"])

        master_path = os.path.join(temp_dir, "master_alignment.csv")
        merge_all_csv_files(temp_dir, master_path)

        assert os.path.exists(master_path)
        with open(master_path, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) >= 1
        assert "holo_resi" in rows[0] or "sequential_position" in str(rows[0])

    def test_merge_includes_is_second_shell_6A(self, temp_dir):
        csp_path = os.path.join(temp_dir, "csp_table.csv")
        with open(csp_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["apo_bmrb", "holo_bmrb", "holo_pdb", "chain", "apo_resi", "apo_aa", "holo_resi", "holo_aa",
                        "H_apo", "N_apo", "H_holo", "N_holo", "dH", "dN", "csp_A", "significant"])
            w.writerow(["18251", "4700", "1cf4", "A", "1", "M", "1", "M", "8.41", "122.0", "8.26", "115.7",
                        "0.04", "-1.5", "0.31", "1"])
            w.writerow(["18251", "4700", "1cf4", "A", "2", "G", "2", "G", "8.20", "110.0", "8.10", "109.0",
                        "0.10", "1.0", "0.20", "0"])

        for filename, headers, rows in [
            (
                "occlusion_analysis.csv",
                ["residue_number", "residue_name", "chain_id", "delta_sasa", "is_occluded"],
                [["1", "MET", "A", "0.02", "True"], ["2", "GLY", "A", "0.0", "False"]],
            ),
            (
                "interaction_filter.csv",
                ["residue_number", "residue_name", "chain_id", "has_hbond"],
                [["1", "MET", "A", "False"], ["2", "GLY", "A", "False"]],
            ),
            (
                "ca_distance_filter.csv",
                ["residue_number", "residue_name", "chain_id", "passes_filter"],
                [["1", "MET", "A", "False"], ["2", "GLY", "A", "False"]],
            ),
            (
                "second_shell_filter.csv",
                [
                    "residue_number", "residue_name", "chain_id", "is_binding_site",
                    "min_distance_to_binding_site", "is_second_shell_6A", "distance_threshold",
                ],
                [
                    ["1", "MET", "A", "True", "0.0000", "False", "6.0000"],
                    ["2", "GLY", "A", "False", "4.0000", "True", "6.0000"],
                ],
            ),
        ]:
            path = os.path.join(temp_dir, filename)
            with open(path, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(headers)
                w.writerows(rows)

        master_path = os.path.join(temp_dir, "master_alignment.csv")
        merge_all_csv_files(temp_dir, master_path)

        with open(master_path, "r") as f:
            reader = csv.DictReader(f)
            master_rows = list(reader)
        assert "is_second_shell_6A" in reader.fieldnames
        by_resi = {r["holo_resi"]: r for r in master_rows}
        assert by_resi["1"]["is_second_shell_6A"] in ("False", "0", "false")
        assert by_resi["2"]["is_second_shell_6A"] in ("True", "1", "true")


class TestCompute1dMetricsForTarget:
    """Tests for compute_1d_metrics_for_target."""

    def test_compute_1d_metrics_creates_csv(self, temp_dir):
        csp_path = os.path.join(temp_dir, "csp_table.csv")
        with open(csp_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["apo_bmrb", "holo_bmrb", "holo_pdb", "chain", "apo_resi", "apo_aa", "holo_resi", "holo_aa",
                        "H_apo", "N_apo", "CA_apo", "H_holo", "N_holo", "CA_holo",
                        "H_offset", "N_offset", "CA_offset", "dH", "dN", "csp_A", "significant"])
            w.writerow(["18251", "4700", "1cf4", "A", "1", "M", "1", "M",
                        "8.41", "122.0", "55.6", "8.26", "115.7", "55.8",
                        "-0.05", "1.1", "0", "0.04", "-1.5", "0.31", "1"])

        compute_1d_metrics_for_target(Path(temp_dir))

        out_path = Path(temp_dir) / "1d_analysis.csv"
        assert out_path.exists()
        with open(out_path, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) >= 1


# --- 3. Mocked Network Tests ---

class TestFetchBmrbMocked:
    """Tests for fetch_bmrb with mocked network."""

    def test_fetch_bmrb_writes_cache(self, temp_dir):
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.content = b"save_test\nloop_\n_Atom_chem_shift.Val\n1.0\nsave_\n"

        with patch("scripts.bmrb_io.requests.get", return_value=mock_resp):
            result = fetch_bmrb("12345", cache_dir=temp_dir, force=True)
            assert os.path.exists(result)
            assert "12345" in result


class TestFetchPdbMocked:
    """Tests for fetch_pdb with mocked network."""

    @patch("scripts.rcsb_io.requests.get")
    def test_fetch_pdb_writes_cache(self, mock_get, temp_dir, minimal_pdb_content):
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.content = minimal_pdb_content.encode()
        mock_get.return_value = mock_resp

        result = fetch_pdb("1abc", cache_dir=temp_dir, force=True)
        assert os.path.exists(result)
        assert "1abc" in result


# --- 4. CSP Computation Tests ---

class TestComputeCspMultipleSaveframes:
    """Tests for compute_csp_multiple_saveframes."""

    def test_compute_csp_returns_results(self):
        # (sequence, H_shifts, N_shifts, CA_shifts, HA_shifts, saveframe_name)
        apo_seq = "MQT"
        holo_seq = "MQT"
        H_apo = {1: 8.41, 2: 8.50, 3: 8.22}
        N_apo = {1: 122.0, 2: 121.7, 3: 117.2}
        CA_apo = {1: 55.0, 2: 56.0, 3: 61.5}
        HA_apo: dict = {}
        H_holo = {1: 8.26, 2: 8.50, 3: 8.26}
        N_holo = {1: 115.7, 2: 121.7, 3: 115.7}
        CA_holo = {1: 55.8, 2: 56.0, 3: 60.6}
        HA_holo: dict = {}

        apo_sequences = [(apo_seq, H_apo, N_apo, CA_apo, HA_apo, "apo1")]
        holo_sequences = [(holo_seq, H_holo, N_holo, CA_holo, HA_holo, "holo1")]

        results = compute_csp_multiple_saveframes(
            apo_sequences, holo_sequences,
            "18251", "4700", "1cf4",
            referencing_method="mean",
        )
        assert isinstance(results, list)
        assert len(results) >= 1
        r = results[0]
        assert hasattr(r, "csp_A")
        assert hasattr(r, "holo_index")
        assert hasattr(r, "significant")

    def test_grid_referencing_paths_use_canonical_target_id(self, tmp_path):
        """Grid offset caches live under output_root/target_id/ (matches pipeline target_label)."""
        output_root = str(tmp_path)
        target_id = "1CF4_18251"
        slug = "h7_h10_n6_n12_c1"
        output_base_dir = output_root
        out_dir = os.path.join(output_base_dir, target_id) if target_id else None
        csv_path = os.path.join(out_dir, f"offset_grid_{slug}.csv") if out_dir else None
        assert csv_path == os.path.join(output_root, "1CF4_18251", f"offset_grid_{slug}.csv")


class TestCaseStudyViewPaths:
    def test_case_study_view_candidates_canonical_then_legacy_pdb(self, tmp_path):
        from scripts.case_study import _case_study_view_load_candidates

        d = str(tmp_path)
        cands = _case_study_view_load_candidates(d, "1cf4", "1CF4_18251")
        assert cands[0] == os.path.join(d, "1CF4_18251_case_study_view.json")
        assert os.path.join(d, "1cf4_case_study_view.json") in cands

    def test_hsqc_offset_panel_title_is_plain(self):
        from scripts.case_study import format_hsqc_offset_panel_title

        # Offsets belong in the HSQC legend only, not the case-study panel title.
        assert format_hsqc_offset_panel_title() == "Apo/Holo HSQC Offset"
        assert format_hsqc_offset_panel_title(-0.04, 1.2) == "Apo/Holo HSQC Offset"

    def test_read_applied_hn_offsets_from_csp_table(self, tmp_path):
        from scripts.case_study import read_applied_hn_offsets_from_csp_table

        csp = tmp_path / "csp_table.csv"
        csp.write_text(
            "apo_resi,H_offset,N_offset\n"
            "1,,\n"
            "2,-0.0400,1.2000\n",
            encoding="utf-8",
        )
        assert read_applied_hn_offsets_from_csp_table(str(tmp_path)) == (-0.04, 1.2)


# --- 5. SASA / Interaction / CA Distance Tests ---

HAS_MDTRAJ = False
try:
    import mdtraj
    HAS_MDTRAJ = True
except ImportError:
    pass

PDB_1CF4_PATH = None
_project_root = Path(__file__).resolve().parent.parent
_candidate = _project_root / "PDB_FILES" / "1cf4.pdb"
if _candidate.exists():
    PDB_1CF4_PATH = str(_candidate)


@pytest.mark.skipif(not HAS_MDTRAJ or not PDB_1CF4_PATH, reason="mdtraj or PDB file not available")
class TestComputeSasaOcclusion:
    """Tests for compute_sasa_occlusion."""

    def test_compute_sasa_returns_dict(self):
        result = compute_sasa_occlusion(
            PDB_1CF4_PATH,
            sasa_threshold=0.0,
            receptor_chain_id="A",
            ligand_chain_id="B",
        )
        assert "residue_info" in result
        assert isinstance(result["residue_info"], list)
        if result["residue_info"]:
            assert "residue_number" in result["residue_info"][0]
            assert "delta_sasa" in result["residue_info"][0]


@pytest.mark.skipif(not HAS_MDTRAJ or not PDB_1CF4_PATH, reason="mdtraj or PDB file not available")
class TestComputeInteractionFilter:
    """Tests for compute_interaction_filter."""

    def test_compute_interaction_returns_dict(self):
        result = compute_interaction_filter(
            PDB_1CF4_PATH,
            distance_threshold=4.5,
            receptor_chain_id="A",
            ligand_chain_id="B",
        )
        assert "residue_info" in result
        assert isinstance(result["residue_info"], list)


@pytest.mark.skipif(not HAS_MDTRAJ or not PDB_1CF4_PATH, reason="mdtraj or PDB file not available")
class TestComputeCaDistanceFilter:
    """Tests for compute_ca_distance_filter."""

    def test_compute_ca_distance_returns_dict(self):
        result = compute_ca_distance_filter(
            PDB_1CF4_PATH,
            distance_threshold=6.0,
            receptor_chain_id="A",
            ligand_chain_id="B",
        )
        assert "residue_info" in result
        assert isinstance(result["residue_info"], list)


class TestSecondShellHelpers:
    """Tests for second-shell naming and binding-site helpers (no structure I/O)."""

    def test_second_shell_column_name(self):
        assert second_shell_column_name(6.0) == "is_second_shell_6A"
        assert second_shell_column_name(6) == "is_second_shell_6A"
        assert second_shell_column_name(6.5) == "is_second_shell_6p5A"

    def test_is_binding_site_row_union(self):
        assert is_binding_site_row({"is_occluded_occlusion": "True"})
        assert is_binding_site_row({"passes_filter_distance": "1"})
        assert is_binding_site_row({"has_hbond_interaction": "True"})
        assert is_binding_site_row({"passes_sub_2A_filter_any_atom": "True"})
        assert not is_binding_site_row({})


@pytest.mark.skipif(not HAS_MDTRAJ, reason="mdtraj not available")
class TestSecondShellFilter:
    """Tests for exclusive second-shell residue mask."""

    def _write_shell_pdb(self, temp_dir: str) -> str:
        # Residue 1 at origin (binding site), residue 2 ~4 Å away (second shell),
        # residue 3 ~20 Å away (neither).
        content = """HEADER    SECOND SHELL TEST
ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  MET A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   MET A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   MET A   1       1.251   2.390   0.000  1.00  0.00           O
ATOM      5  N   GLY A   2       4.000   0.000   0.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       5.000   0.000   0.000  1.00  0.00           C
ATOM      7  C   GLY A   2       6.000   0.000   0.000  1.00  0.00           C
ATOM      8  O   GLY A   2       7.000   0.000   0.000  1.00  0.00           O
ATOM      9  N   ALA A   3      20.000   0.000   0.000  1.00  0.00           N
ATOM     10  CA  ALA A   3      21.000   0.000   0.000  1.00  0.00           C
ATOM     11  C   ALA A   3      22.000   0.000   0.000  1.00  0.00           C
ATOM     12  O   ALA A   3      23.000   0.000   0.000  1.00  0.00           O
END
"""
        path = os.path.join(temp_dir, "shell_test.pdb")
        with open(path, "w") as f:
            f.write(content)
        return path

    def test_exclusive_second_shell_logic(self, temp_dir):
        pdb_path = self._write_shell_pdb(temp_dir)
        result = compute_second_shell_filter(
            pdb_path,
            binding_site_residue_numbers=[1],
            distance_threshold=6.0,
            receptor_chain_id="A",
        )
        assert "error" not in result
        by_res = {info["residue_number"]: info for info in result["residue_info"]}
        assert by_res[1]["is_binding_site"] is True
        assert by_res[1]["is_second_shell_6A"] is False
        assert by_res[2]["is_binding_site"] is False
        assert by_res[2]["is_second_shell_6A"] is True
        assert by_res[3]["is_binding_site"] is False
        assert by_res[3]["is_second_shell_6A"] is False

    def test_empty_binding_site_all_false(self, temp_dir):
        pdb_path = self._write_shell_pdb(temp_dir)
        result = compute_second_shell_filter(
            pdb_path,
            binding_site_residue_numbers=[],
            distance_threshold=6.0,
            receptor_chain_id="A",
        )
        assert "error" not in result
        assert all(not info["is_second_shell_6A"] for info in result["residue_info"])
        assert all(not info["is_binding_site"] for info in result["residue_info"])

    def test_write_second_shell_csv(self, temp_dir):
        pdb_path = self._write_shell_pdb(temp_dir)
        result = compute_second_shell_filter(
            pdb_path,
            binding_site_residue_numbers=[1],
            distance_threshold=6.0,
            receptor_chain_id="A",
        )
        out_path = os.path.join(temp_dir, "second_shell_filter.csv")
        write_second_shell_csv(result["residue_info"], out_path, distance_threshold=6.0)
        with open(out_path, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert "is_second_shell_6A" in reader.fieldnames
        assert len(rows) == 3


# --- 5b. HN significance mask helpers ---

class TestHnSignificanceMasks:
    """Tests for alternate HN CSP significance boolean masks."""

    def _make_result(self, idx: int, csp_A: Optional[float]) -> CSPResult:
        return CSPResult(
            apo_index=idx,
            holo_index=idx,
            apo_aa="A",
            holo_aa="A",
            H_apo=8.0,
            N_apo=120.0,
            H_holo=8.0,
            N_holo=120.0,
            dH=0.0,
            dN=0.0,
            csp_A=csp_A,
            significant=None,
        )

    def test_percentile_linear_basic(self):
        assert percentile_linear([1.0], 90.0) == 1.0
        assert percentile_linear([0.0, 10.0], 50.0) == 5.0
        assert percentile_linear([0.0, 10.0], 90.0) == 9.0

    def test_sigma_masks_match_cleaned_mean_and_sd(self):
        # Values chosen so outlier removal with z=3 leaves a stable cleaned set.
        values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.50]
        stats = compute_threshold_with_outlier_removal(
            values, outlier_z=3.0, significance_z=0.0, max_iterations=10, max_outlier_fraction=0.2
        )
        results = [self._make_result(i + 1, v) for i, v in enumerate(values)]
        _apply_hn_significance_masks(
            results,
            cutoff=stats.threshold,
            final_stats=stats,
            csp_values=values,
        )

        mean0 = stats.mean
        thr1 = stats.mean + stats.sd
        thr2 = stats.mean + 2.0 * stats.sd

        thr_primary = max(0.05, mean0)
        for r, v in zip(results, values):
            assert r.significant_sigma_0 == (v >= mean0)
            assert r.significant_sigma_1 == (v >= thr1)
            assert r.significant_sigma_2 == (v >= thr2)
            assert r.significant_1sd == r.significant_sigma_1
            assert r.significant_2sd == r.significant_sigma_2
            # Primary uses floored cutoff; this helper receives cutoff as-is.
            assert r.significant == (v >= stats.threshold)
            assert r.significant_max_05_cleaned_mean == (v >= thr_primary)

    def test_percentile_and_ppm_masks(self):
        values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.20]
        stats = compute_threshold_with_outlier_removal(
            values, outlier_z=3.0, significance_z=0.0, max_iterations=10, max_outlier_fraction=0.2
        )
        results = [self._make_result(i + 1, v) for i, v in enumerate(values)]
        _apply_hn_significance_masks(
            results,
            cutoff=stats.threshold,
            final_stats=stats,
            csp_values=values,
        )

        top10 = percentile_linear(values, 90.0)
        top5 = percentile_linear(values, 95.0)
        assert results[-1].significant_top_10_percentile is True
        assert results[-1].significant_top_5_percentile is True
        assert results[0].significant_top_10_percentile is False
        assert all(r.significant_top_10_percentile == (r.csp_A >= top10) for r in results)
        assert all(r.significant_top_5_percentile == (r.csp_A >= top5) for r in results)

        assert results[0].significant_03_ppm is False  # 0.01
        assert results[2].significant_03_ppm is True   # 0.03
        assert results[3].significant_05_ppm is False  # 0.04
        assert results[4].significant_05_ppm is True   # 0.05
        assert results[8].significant_10_ppm is False  # 0.09
        assert results[9].significant_10_ppm is True   # 0.20

    def test_none_csp_leaves_masks_none(self):
        stats = compute_threshold_with_outlier_removal(
            [0.01, 0.02, 0.03],
            outlier_z=3.0,
            significance_z=0.0,
            max_iterations=10,
            max_outlier_fraction=0.2,
        )
        results = [self._make_result(1, None), self._make_result(2, 0.05)]
        _apply_hn_significance_masks(
            results,
            cutoff=stats.threshold,
            final_stats=stats,
            csp_values=[0.01, 0.02, 0.03],
        )
        assert results[0].significant is None
        assert results[0].significant_sigma_0 is None
        assert results[0].significant_03_ppm is None
        assert results[0].significant_raw_mean is None
        assert results[0].significant_max_05_cleaned_mean is None
        assert results[1].significant_03_ppm is True

    def test_raw_mean_and_floor_masks(self):
        # Tight cluster + extreme outlier so cleaned mean drops below 0.05 ppm.
        values = [0.01, 0.012, 0.014, 0.016, 0.018, 0.020, 0.022, 0.024, 0.026, 2.0]
        stats = compute_threshold_with_outlier_removal(
            values, outlier_z=2.0, significance_z=0.0, max_iterations=10, max_outlier_fraction=0.2
        )
        results = [self._make_result(i + 1, v) for i, v in enumerate(values)]
        _apply_hn_significance_masks(
            results,
            cutoff=stats.threshold,
            final_stats=stats,
            csp_values=values,
        )

        raw_mean, raw_sd = _compute_basic_stats(values)
        thr_raw_1 = raw_mean + raw_sd
        thr_raw_2 = raw_mean + 2.0 * raw_sd
        thr_floor = max(0.05, stats.mean)

        assert stats.outliers_removed >= 1
        assert stats.mean < raw_mean
        assert stats.mean < 0.05
        assert thr_floor == 0.05

        for r, v in zip(results, values):
            assert r.significant_raw_mean == (v >= raw_mean)
            assert r.significant_raw_1sd == (v >= thr_raw_1)
            assert r.significant_raw_2sd == (v >= thr_raw_2)
            assert r.significant_max_05_cleaned_mean == (v >= thr_floor)
            # Unfloored cutoff still passed here; primary floor is applied upstream.
            assert r.significant == (v >= stats.threshold)
            assert r.significant_sigma_0 == (v >= stats.mean)

    def test_floor_primary_hn_cutoff_helper(self):
        assert _floor_primary_hn_cutoff(0.02) == 0.05
        assert _floor_primary_hn_cutoff(0.05) == 0.05
        assert _floor_primary_hn_cutoff(0.12) == 0.12


class TestAaMismatchCspGate:
    """CSPs require identical apo/holo amino acids."""

    def test_aa_mismatch_excludes_csp(self):
        results = compute_csp_A(
            mapping=[(1, 1), (2, 2)],
            apo_seq="GL",
            holo_seq="EL",
            H_apo={1: 8.24, 2: 7.67},
            N_apo={1: 108.3, 2: 120.4},
            H_holo={1: 8.42, 2: 8.09},
            N_holo={1: 119.2, 2: 121.7},
            threshold_config={"absolute_cutoff": 0.05},
            enable_referencing=False,
        )
        assert results[0].excluded_aa_mismatch is True
        assert results[0].csp_A is None
        assert results[0].significant is None
        assert results[1].excluded_aa_mismatch is False
        assert results[1].csp_A is not None
        assert results[1].significant is True
        assert sum(1 for r in results if r.csp_A is not None) == 1


# --- 6. Visualization Tests ---

class TestPlotHsqcVariants:
    """Tests for plot_hsqc_variants."""

    def test_plot_hsqc_creates_file(self, temp_dir):
        results = [
            CSPResult(1, 1, "M", "M", 8.41, 122.0, 8.26, 115.7, 0.04, -1.5, 0.31, True,
                      H_holo_original=8.26, N_holo_original=115.7, H_offset=-0.05, N_offset=1.1),
        ]
        out_path = os.path.join(temp_dir, "hsqc.png")
        returned = plot_hsqc_variants(results, out_path, title="Test")
        assert os.path.exists(out_path)
        assert os.path.getsize(out_path) > 0
        assert returned == (-0.05, 1.1)

    def test_plot_hsqc_prefers_applied_offsets_over_grid(self, temp_dir):
        """Stored H_offset/N_offset must be used even if a grid search would differ."""
        results = [
            CSPResult(
                1,
                1,
                "M",
                "M",
                8.41,
                122.0,
                8.36,
                123.1,
                0.0,
                0.0,
                0.05,
                False,
                H_holo_original=8.41,
                N_holo_original=122.0,
                H_offset=-0.123,
                N_offset=0.456,
            ),
            CSPResult(
                2,
                2,
                "Q",
                "Q",
                8.10,
                120.0,
                8.05,
                119.5,
                0.0,
                0.0,
                0.05,
                False,
                H_holo_original=8.10,
                N_holo_original=120.0,
                H_offset=-0.123,
                N_offset=0.456,
            ),
        ]
        assert resolve_hsqc_offsets(results, "H", "N") == (-0.123, 0.456)
        out_path = os.path.join(temp_dir, "hsqc_applied.png")
        returned = plot_hsqc_variants(results, out_path, title="Applied offsets")
        assert returned == (-0.123, 0.456)
        assert format_offset_annotation("H", -0.123, "N", 0.456) == (
            "ΔH=-0.123 ppm, ΔN=+0.456 ppm"
        )

    def test_plot_hsqc_with_binding_results(self, temp_dir):
        results = [
            CSPResult(1, 1, "M", "M", 8.41, 122.0, 8.26, 115.7, 0.04, -1.5, 0.31, True,
                      H_holo_original=8.26, N_holo_original=115.7, H_offset=-0.05, N_offset=1.1),
            CSPResult(2, 2, "Q", "Q", 8.10, 120.0, 8.05, 119.5, 0.02, -0.2, 0.05, False,
                      H_holo_original=8.05, N_holo_original=119.5, H_offset=0.0, N_offset=0.0),
        ]
        binding_results = {
            "residue_info": [
                {"residue_number": 1, "residue_name": "MET", "has_hbond": True,
                 "has_charge_complement": False, "has_pi_contact": False,
                 "has_sasa_occlusion": False, "has_ca_distance": False,
                 "has_any_atom_sub_2A": False},
                {"residue_number": 2, "residue_name": "GLN", "has_hbond": False,
                 "has_charge_complement": False, "has_pi_contact": False,
                 "has_sasa_occlusion": False, "has_ca_distance": False,
                 "has_any_atom_sub_2A": False},
            ],
            "dataset_type": "union",
            "n_union_residues": 1,
        }
        out_path = os.path.join(temp_dir, "hsqc_classified.png")
        returned = plot_hsqc_variants(
            results,
            out_path,
            title="Test classified",
            binding_results=binding_results,
        )
        assert os.path.exists(out_path)
        assert os.path.getsize(out_path) > 0
        assert returned == (-0.05, 1.1)

    def test_missing_csp_binding_residue_is_excluded_gray_not_fn(self):
        """Residues with significant=None must not be labeled FN on HSQC overlays."""
        res = CSPResult(
            63,
            66,
            "Q",
            "Q",
            9.37,
            125.23,
            9.789,
            106.397,
            0.419,
            -18.833,
            None,
            None,
            H_holo_original=9.789,
            N_holo_original=106.397,
            excluded_large_dn=True,
        )
        cls = _classify_residue(
            res,
            position_map={66: 66},
            binding_lookup={66: True},
            significance_field="significant",
        )
        assert cls == "EXCLUDED"
        assert cls != "FN"
        assert _classification_color_map()[cls] == _EXCLUDED_LINE_COLOR

    def test_extreme_csp_z_is_excluded_gray(self):
        """Residues with csp_z > max_classification_csp_z use gray HSQC connectors."""
        res = CSPResult(
            55,
            55,
            "K",
            "K",
            8.0,
            120.0,
            8.1,
            126.3,
            0.1,
            6.3,
            0.63,
            True,
            H_holo_original=8.1,
            N_holo_original=126.3,
            z_score=264.0,
        )
        cls = _classify_residue(
            res,
            position_map={55: 55},
            binding_lookup={55: False},
            significance_field="significant",
        )
        assert cls == "EXCLUDED"
        assert _classification_color_map()[cls] == _EXCLUDED_LINE_COLOR


class TestPlotCspClassificationBars:
    """Tests for plot_csp_classification_bars."""

    def test_plot_classification_creates_file(self, temp_dir):
        results = [
            CSPResult(1, 1, "M", "M", 8.41, 122.0, 8.26, 115.7, 0.04, -1.5, 0.31, True,
                      significant_1sd=True, significant_2sd=False),
        ]
        interaction_results = {
            "residue_info": [
                {"residue_number": 1, "residue_name": "MET", "has_hbond": False,
                 "has_charge_complement": False, "has_pi_contact": False,
                 "has_sasa_occlusion": True, "has_ca_distance": False},
            ],
            "n_interacting_residues": 0,
            "n_union_residues": 1,
        }
        out_path = os.path.join(temp_dir, "classification.png")
        plot_csp_classification_bars(
            results, interaction_results, out_path,
            title="Test", significance_field="significant",
        )
        assert os.path.exists(out_path)
        assert os.path.getsize(out_path) > 0


class TestWritePymolScripts:
    """Tests for PyMOL script writers."""

    def test_write_color_csp_mask_script(self, temp_dir):
        results = [
            CSPResult(1, 1, "M", "M", 8.41, 122.0, 8.26, 115.7, 0.04, -1.5, 0.31, True),
        ]
        sasa_results = {
            "residue_info": [{"residue_number": 1, "residue_name": "MET", "delta_sasa": 0.02, "is_occluded": True}],
            "receptor_chain": "A",
            "ligand_chain": "B",
        }
        out_path = os.path.join(temp_dir, "color_csp_mask.pml")
        write_pymol_color_csp_mask_script(
            results, "1cf4", out_path,
            sasa_results=sasa_results,
            receptor_chain="A",
            ligand_chain="B",
            output_dir="./",
        )
        assert os.path.exists(out_path)
        with open(out_path, "r") as f:
            content = f.read()
        assert "cmd.load" in content or "load" in content

    def test_write_occlusion_script(self, temp_dir):
        sasa_results = {
            "residue_info": [{"residue_number": 1, "residue_name": "MET", "delta_sasa": 0.02, "is_occluded": True}],
            "receptor_chain": "A",
            "ligand_chain": "B",
        }
        out_path = os.path.join(temp_dir, "color_occlusion.pml")
        write_pymol_occlusion_script(
            sasa_results, "1cf4", out_path,
            receptor_chain="A",
            ligand_chain="B",
            output_dir="./",
        )
        assert os.path.exists(out_path)

    def test_write_delta_sasa_script(self, temp_dir):
        sasa_results = {
            "residue_info": [{"residue_number": 1, "residue_name": "MET", "delta_sasa": 0.02, "is_occluded": True}],
            "receptor_chain": "A",
            "ligand_chain": "B",
        }
        out_path = os.path.join(temp_dir, "delta_sasa.pml")
        write_pymol_delta_sasa_script(sasa_results, "1cf4", out_path, output_dir="./")
        assert os.path.exists(out_path)


# --- Missing CSP exclusion ---

class TestMissingCspClassification:
    """Residues without a recorded CSP must not receive TP/FP/TN/FN labels."""

    def _base_row(self, **overrides):
        row = {
            "significant": "0",
            "csp_A": "0.05",
            "is_occluded_occlusion": "False",
            "passes_filter_distance": "False",
            "has_hbond_interaction": "False",
            "has_charge_complement_interaction": "False",
            "has_pi_contact_interaction": "False",
            "passes_sub_2A_filter_any_atom": "False",
        }
        row.update(overrides)
        return row

    def test_parse_optional_bool_tri_state(self):
        assert parse_optional_bool("") is None
        assert parse_optional_bool(None) is None
        assert parse_optional_bool("nan") is None
        assert parse_optional_bool(0) is False
        assert parse_optional_bool("0") is False
        assert parse_optional_bool(1) is True
        assert parse_optional_bool("1") is True

    def test_has_recorded_csp_requires_numeric_csp_a(self):
        assert has_recorded_csp({"csp_A": "0.12", "significant": "0"})
        assert not has_recorded_csp({"csp_A": "", "significant": "0"})
        assert not has_recorded_csp({"csp_A": "", "significant": ""})
        # Fallback when csp_A column absent
        assert has_recorded_csp({"significant": "1"}, csp_column="missing_col")
        assert not has_recorded_csp({"significant": ""}, csp_column="missing_col")

    def test_compute_classification_empty_when_significant_missing(self):
        assert compute_classification(self._base_row(significant="", csp_A="")) == ""
        assert compute_classification(self._base_row(significant="")) == ""

    def test_compute_classification_empty_when_csp_z_above_max(self):
        assert (
            compute_classification(
                self._base_row(significant="1", csp_A="0.63", csp_z="264.0")
            )
            == ""
        )
        assert (
            compute_classification(
                self._base_row(significant="1", csp_A="0.63", csp_z="100.0")
            )
            == "FP"
        )

    def test_compute_classification_tn_when_recorded_nonsignificant(self):
        assert compute_classification(self._base_row(significant="0", csp_A="0.01")) == "TN"

    def test_compute_classification_fn_when_binding_and_nonsignificant(self):
        assert (
            compute_classification(
                self._base_row(significant="0", csp_A="0.01", is_occluded_occlusion="True")
            )
            == "FN"
        )

    def test_filter_recorded_csp_dataframe_drops_blank_csp(self):
        import pandas as pd

        df = pd.DataFrame(
            {
                "csp_A": ["0.1", "", "0.2"],
                "significant": ["1", "", "0"],
            }
        )
        filtered = filter_recorded_csp_dataframe(df)
        assert len(filtered) == 2
        assert list(filtered["csp_A"]) == ["0.1", "0.2"]


class TestMissingCspExcludedFromFiguresAndF1:
    """Figure 3 / F1 loaders must not count blank-CSP rows as TN/FN."""

    def test_collect_distance_categories_skips_missing_csp(self, temp_dir):
        import pandas as pd

        target_dir = Path(temp_dir) / "TEST_1"
        target_dir.mkdir()
        df = pd.DataFrame(
            {
                "significant": ["1", "0", ""],
                "csp_A": ["0.20", "0.01", ""],
                "min_ca_distance_distance": [3.0, 8.0, 4.0],
                "passes_filter_distance": [True, False, True],
                "has_charge_complement_interaction": [False, False, False],
                "has_pi_contact_interaction": [False, False, False],
                "has_hbond_interaction": [False, False, False],
                "is_occluded_occlusion": [False, False, False],
            }
        )
        df.to_csv(target_dir / "master_alignment.csv", index=False)

        tp, fp, fn, tn = collect_distance_categories(Path(temp_dir))
        # Row0: sig+binding → TP; Row1: nonsig+nonbinding → TN; Row2 missing CSP excluded
        assert len(tp) == 1 and tp[0] == pytest.approx(3.0)
        assert len(tn) == 1 and tn[0] == pytest.approx(8.0)
        assert fp == []
        assert fn == []

    def test_si_fig_s13_s14_match_figure_3_classification(self, temp_dir):
        """S13/S14 must use parse_optional_bool + recorded-CSP filter like Figure 3.

        Pandas stores ``significant`` as float64 (0.0/1.0/NaN) whenever any row is
        blank; the old ``str(value) in {"1", "true", ...}`` parser treated 1.0 as
        False and dumped those residues into FN.
        """
        import pandas as pd

        target_dir = Path(temp_dir) / "TEST_1"
        target_dir.mkdir()
        df = pd.DataFrame(
            {
                "holo_pdb": ["1abc", "1abc", "1abc"],
                "pdb_residue_number": [10, 11, 12],
                "significant": [1.0, 0.0, float("nan")],
                "csp_A": [0.20, 0.01, float("nan")],
                "min_ca_distance_distance": [3.0, 8.0, 4.0],
                "min_nn_distance_nn_distance": [3.5, 8.5, 4.5],
                "min_any_atom_distance_any_atom": [2.5, 7.5, 3.5],
                "passes_filter_distance": [True, False, True],
                "has_charge_complement_interaction": [False, False, False],
                "has_pi_contact_interaction": [False, False, False],
                "has_hbond_interaction": [False, False, False],
                "is_occluded_occlusion": [False, False, False],
            }
        )
        df.to_csv(target_dir / "master_alignment.csv", index=False)

        fig3 = collect_distance_categories(Path(temp_dir))
        s13 = collect_s13_distance_categories(Path(temp_dir))
        s14 = collect_s14_distance_categories(Path(temp_dir))
        for label, result in (("fig3", fig3), ("s13", s13), ("s14", s14)):
            tp, fp, fn, tn = result
            assert len(tp) == 1, label
            assert fp == []
            assert fn == []
            assert len(tn) == 1, label
        assert fig3[0][0] == pytest.approx(3.0)
        assert s13[0][0] == pytest.approx(3.5)
        assert s14[0][0] == pytest.approx(2.5)
        assert fig3[3][0] == pytest.approx(8.0)
        assert s13[3][0] == pytest.approx(8.5)
        assert s14[3][0] == pytest.approx(7.5)

    def test_load_alignment_and_f1_exclude_missing_csp(self, temp_dir):
        import pandas as pd

        alignment_path = Path(temp_dir) / "master_alignment.csv"
        df = pd.DataFrame(
            {
                "significant": ["1", "0", ""],
                "csp_A": ["0.20", "0.01", ""],
                "min_ca_distance_distance": [3.0, 8.0, 4.0],
                "passes_filter_distance": [True, False, True],
                "has_charge_complement_interaction": [False, False, False],
                "has_pi_contact_interaction": [False, False, False],
                "has_hbond_interaction": [False, False, False],
                "is_occluded_occlusion": [False, False, False],
            }
        )
        df.to_csv(alignment_path, index=False)

        loaded = load_alignment(alignment_path)
        assert len(loaded) == 2

        metrics = compute_f1_score(loaded)
        # TP=1 (sig+binding), TN=1 (nonsig+nonbinding); missing row must not inflate TN/FN
        assert metrics.true_positives == 1
        assert metrics.false_positives == 0
        assert metrics.false_negatives == 0
        assert metrics.total_rows == 2


class TestApoApoBufferAndEntities:
    """NMR-STAR v2.1 sample loops use _Mol_label / _Concentration_value."""

    def test_v21_buffer_and_entities(self, temp_dir):
        from scripts.find_apo_apo_matches import (
            extract_buffer_text,
            extract_entity_names,
            enrich_pairs_csv,
        )

        star = Path(temp_dir) / "1_21.str"
        star.write_text(
            """
data_1
save_assembly
   _Saveframe_category         molecular_system
   loop_
      _Mol_system_component_name
      _Mol_label
      'PROTEIN (UBIQUITIN)' $Ubiquitin
   stop_
save_

save_Ubiquitin
   _Saveframe_category                          monomeric_polymer
   _Mol_type                                    polymer
   _Mol_polymer_class                           protein
   _Name_common                                 ubiquitin
   _Residue_count                               76
save_

save_sample
   _Saveframe_category   sample
   loop_
      _Mol_label
      _Concentration_value
      _Concentration_value_units
      $Ubiquitin           . mM '[U-99% 13C; U-99% 15N]'
      'Phosphate Buffer' 50 mM 'natural abundance'
       D2O                7 %  'natural abundance'
   stop_
save_

save_sample_conditions_1
   _Saveframe_category   sample_conditions
   loop_
      _Variable_type
      _Variable_value
      _Variable_value_units
      pH 6.5 .
      temperature 298 K
      ionic_strength 50 mM
   stop_
save_
""",
            encoding="utf-8",
        )
        buf = extract_buffer_text(star)
        assert "Phosphate Buffer=50" in buf and "mM" in buf
        assert "D2O=7" in buf
        assert "PROTEIN (UBIQUITIN)" in buf or "Ubiquitin" in buf
        names = extract_entity_names(star)
        assert "ubiquitin" in names.lower()

        pairs = Path(temp_dir) / "apo_apo_pairs.csv"
        with pairs.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(
                f,
                fieldnames=["query_apo_bmrb", "match_apo_bmrb", "query_star_path", "match_star_path"],
            )
            w.writeheader()
            w.writerow(
                {
                    "query_apo_bmrb": "1",
                    "match_apo_bmrb": "1",
                    "query_star_path": str(star),
                    "match_star_path": str(star),
                }
            )
        dest = enrich_pairs_csv(pairs, cs_dir=Path(temp_dir), out_csv=Path(temp_dir) / "out.csv")
        with dest.open(newline="", encoding="utf-8") as f:
            row = next(csv.DictReader(f))
        assert row["query_pH"] == "6.5"
        assert "Phosphate Buffer=50" in row["query_buffer"]
        assert row["query_entities"]
        assert row["match_entities"] == row["query_entities"]


# --- 7. Integration Test ---

@pytest.mark.integration
@pytest.mark.skipif(
    not os.path.exists(os.path.join(Path(__file__).resolve().parent.parent, "data/CSP_UBQ.csv")),
    reason="CSP_UBQ.csv not found",
)
class TestProcessRowIntegration:
    """Integration test for process_row (requires network, slow)."""

    def test_process_row_1cf4(self, temp_dir):
        row = {"apo_bmrb": "18251", "holo_bmrb": "4700", "holo_pdb": "1cf4"}
        process_row(
            row,
            temp_dir,
            generate_case_study=False,
            receptor_msa_png=False,
        )
        tgt_dir = os.path.join(temp_dir, "1CF4_18251")
        assert os.path.exists(tgt_dir)
        assert os.path.exists(os.path.join(tgt_dir, "csp_table.csv"))
        assert os.path.exists(os.path.join(tgt_dir, "occlusion_analysis.csv"))
        assert os.path.exists(os.path.join(tgt_dir, "master_alignment.csv"))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
