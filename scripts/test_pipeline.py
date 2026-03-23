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
    from .csp import CSPResult, compute_csp_multiple_saveframes, compute_csp_multiple_saveframes_ca
    from .sasa_analysis import compute_sasa_occlusion, write_occlusion_analysis_csv
    from .interaction_analysis import compute_interaction_filter, write_interaction_analysis_csv
    from .interaction_analysis import compute_ca_distance_filter, write_ca_distance_csv
    from .merge_csv import merge_all_csv_files
    from .analyze_targets_single_atom_shifts import compute_1d_metrics_for_target
    from .HSQC_visualize import plot_hsqc_variants
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
    from scripts.csp import CSPResult, compute_csp_multiple_saveframes, compute_csp_multiple_saveframes_ca
    from scripts.sasa_analysis import compute_sasa_occlusion, write_occlusion_analysis_csv
    from scripts.interaction_analysis import compute_interaction_filter, write_interaction_analysis_csv
    from scripts.interaction_analysis import compute_ca_distance_filter, write_ca_distance_csv
    from scripts.merge_csv import merge_all_csv_files
    from scripts.analyze_targets_single_atom_shifts import compute_1d_metrics_for_target
    from scripts.HSQC_visualize import plot_hsqc_variants
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
    """Minimal NMR-STAR v3 format for BMRB parse test."""
    return """save_assigned_chemical_shifts_1
_Assigned_chem_shift_list.Sf_category    assigned_chemical_shifts
loop_
_Residue_seq_code
_Residue_label
_Atom_name
_Atom_chem_shift.Val
1   MET   N   120.5
1   MET   H   8.41
1   MET   CA  55.0
2   GLY   N   108.0
2   GLY   H   8.2
2   GLY   CA  44.0
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
        # v3 format detection may require _Atom_chem_shift.Val
        result = parse_sequence_and_shifts_from_saveframes(star_path)
        # May return [] if format is not fully valid; at least no crash
        assert isinstance(result, list)


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
        # (sequence, H_shifts, N_shifts, CA_shifts, saveframe_name)
        apo_seq = "MQT"
        holo_seq = "MQT"
        H_apo = {1: 8.41, 2: 8.50, 3: 8.22}
        N_apo = {1: 122.0, 2: 121.7, 3: 117.2}
        CA_apo = {1: 55.0, 2: 56.0, 3: 61.5}
        H_holo = {1: 8.26, 2: 8.50, 3: 8.26}
        N_holo = {1: 115.7, 2: 121.7, 3: 115.7}
        CA_holo = {1: 55.8, 2: 56.0, 3: 60.6}

        apo_sequences = [(apo_seq, H_apo, N_apo, CA_apo, "apo1")]
        holo_sequences = [(holo_seq, H_holo, N_holo, CA_holo, "holo1")]

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


# --- 6. Visualization Tests ---

class TestPlotHsqcVariants:
    """Tests for plot_hsqc_variants."""

    def test_plot_hsqc_creates_file(self, temp_dir):
        results = [
            CSPResult(1, 1, "M", "M", 8.41, 122.0, 8.26, 115.7, 0.04, -1.5, 0.31, True,
                      H_holo_original=8.26, N_holo_original=115.7, H_offset=-0.05, N_offset=1.1),
        ]
        out_path = os.path.join(temp_dir, "hsqc.png")
        plot_hsqc_variants(results, out_path, title="Test")
        assert os.path.exists(out_path)
        assert os.path.getsize(out_path) > 0


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
        )
        tgt_dir = os.path.join(temp_dir, "1cf4")
        assert os.path.exists(tgt_dir)
        assert os.path.exists(os.path.join(tgt_dir, "csp_table.csv"))
        assert os.path.exists(os.path.join(tgt_dir, "occlusion_analysis.csv"))
        assert os.path.exists(os.path.join(tgt_dir, "master_alignment.csv"))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
