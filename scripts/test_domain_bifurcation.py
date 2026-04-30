"""Unit tests for domain vs full-length bifurcation (offline, no blastp)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from scripts.domain_full_length_bifurcation import (
    _classify_row,
    _pick_uniprot_accession,
    process_dataset,
)
from scripts.uniprot_blast import BlastHit
from scripts.uniprot_io import normalize_uniprot_accession


@pytest.mark.parametrize(
    "holo_len,uni,ratio,expected",
    [(100, 200, 1.5, "domain"), (100, 100, 1.5, "full_length"), (100, 150, 1.5, "full_length")],
)
def test_classify_row(holo_len, uni, ratio, expected):
    assert _classify_row(holo_len, uni, ratio) == expected


def test_normalize_uniprot_accession():
    assert normalize_uniprot_accession("sp|P12345|NAME") == "P12345"
    assert normalize_uniprot_accession("P38398-2") == "P38398-2"


def test_pick_prefers_deposited(monkeypatch):
    hits = [
        BlastHit(
            subject_accession="P99999",
            pident=95.0,
            length=90,
            bitscore=180.0,
            qstart=1,
            qend=90,
            sstart=10,
            send=99,
            subject_title="noise",
        ),
        BlastHit(
            subject_accession="P12345",
            pident=92.0,
            length=85,
            bitscore=170.0,
            qstart=1,
            qend=85,
            sstart=5,
            send=89,
            subject_title="prefer",
        ),
    ]

    def _no_tax(_acc):
        return None

    monkeypatch.setattr(
        "scripts.domain_full_length_bifurcation.fetch_taxonomy_id",
        _no_tax,
    )

    bh, why = _pick_uniprot_accession(
        hits,
        holo_tax_id=9606,
        deposited=["P12345"],
        query_len=90,
        tax_cache={},
        min_pident=25.0,
        min_query_cov=0.4,
    )
    assert bh.subject_accession == "P12345"
    assert why == "matched_rcsb_reference"


@pytest.fixture
def fake_csv(tmp_path):
    inp = tmp_path / "in.csv"
    inp.write_text(
        "apo_bmrb,holo_bmrb,apo_pdb,holo_pdb,ec_classes,scope_fold_type\n"
        "99999,88888,,1d5g,hydrolases,all beta\n",
        encoding="utf-8",
    )
    rc = tmp_path / "rec.csv"
    rc.write_text(
        "apo_pdb,holo_pdb,apo_bmrb,holo_bmrb,receptor_chain_id\n"
        ",1d5g,99999,88888,A\n",
        encoding="utf-8",
    )
    return inp, rc


@patch("scripts.domain_full_length_bifurcation.run_blastp")
@patch("scripts.domain_full_length_bifurcation.fetch_holo_polymer_context")
@patch("scripts.domain_full_length_bifurcation._extract_sequence_bmrb_then_pdb")
def test_process_dataset_stubbed(
    mock_extract,
    mock_rcsb,
    mock_blast,
    fake_csv,
    tmp_path,
):
    inp, rc = fake_csv
    mock_extract.side_effect = [
        ("M" * 120, "bmrb_aggregate", False),
        ("M" * 118, "bmrb_aggregate", False),
    ]
    ctx = MagicMock()
    ctx.ncbi_taxonomy_id = None
    ctx.uniprot_accessions = []
    ctx.preferred_uniprot_accession = None
    ctx.entity_id = None
    mock_rcsb.return_value = ctx
    mock_blast.return_value = [
        BlastHit(
            subject_accession="Q12923",
            pident=97.0,
            length=118,
            bitscore=230.0,
            qstart=1,
            qend=118,
            sstart=400,
            send=517,
            subject_title="x",
        ),
    ]

    def _fake_fasta(acc):
        assert acc == "Q12923"
        return "M" * 500, 500

    with patch(
        "scripts.domain_full_length_bifurcation.fetch_fasta_sequence",
        _fake_fasta,
    ):
        counts = process_dataset(
            input_csv=inp,
            receptor_csv=rc,
            out_dir=tmp_path,
            basename="tst",
            cs_dir=tmp_path,
            pdb_cache_dir=str(tmp_path),
            fetch_bmrb=False,
            blast_db="stub",
            blast_remote=False,
            min_pident=20.0,
            min_query_cov=0.4,
            ratio=1.5,
            sleep_between_uniprot=0.0,
            dry_run_rest=False,
        )

    assert mock_blast.call_count == 1
    assert counts["domains"] == 1
    dom = (tmp_path / "tst_domains.csv").read_text(encoding="utf-8")
    assert "domain" in dom.lower() or "Q12923" in dom


@patch("scripts.domain_full_length_bifurcation.run_blastp")
@patch("scripts.domain_full_length_bifurcation.fetch_holo_polymer_context")
@patch("scripts.domain_full_length_bifurcation._extract_sequence_bmrb_then_pdb")
def test_process_dataset_rcsb_uniprot_skips_blast(
    mock_extract,
    mock_rcsb,
    mock_blast,
    fake_csv,
    tmp_path,
):
    inp, rc = fake_csv
    mock_extract.side_effect = [
        ("M" * 120, "bmrb_aggregate", False),
        ("M" * 118, "bmrb_aggregate", False),
    ]
    ctx = MagicMock()
    ctx.ncbi_taxonomy_id = 9606
    ctx.entity_id = "1"
    ctx.uniprot_accessions = ["Q12923"]
    ctx.preferred_uniprot_accession = "Q12923"
    mock_rcsb.return_value = ctx

    def _fake_fasta(acc):
        assert acc == "Q12923"
        return "M" * 500, 500

    with patch(
        "scripts.domain_full_length_bifurcation.fetch_fasta_sequence",
        _fake_fasta,
    ):
        counts = process_dataset(
            input_csv=inp,
            receptor_csv=rc,
            out_dir=tmp_path,
            basename="tst_rcsb",
            cs_dir=tmp_path,
            pdb_cache_dir=str(tmp_path),
            fetch_bmrb=False,
            blast_db="stub",
            blast_remote=False,
            min_pident=20.0,
            min_query_cov=0.4,
            ratio=1.5,
            sleep_between_uniprot=0.0,
            dry_run_rest=False,
        )

    mock_blast.assert_not_called()
    assert counts["domains"] == 1
    all_txt = (tmp_path / "tst_rcsb_bifurcation_all_rows.csv").read_text(encoding="utf-8")
    assert "rcsb" in all_txt.lower()
    assert "n/a" in all_txt.lower()


@patch("scripts.domain_full_length_bifurcation.run_blastp")
@patch("scripts.domain_full_length_bifurcation.fetch_holo_polymer_context")
@patch("scripts.domain_full_length_bifurcation._extract_sequence_bmrb_then_pdb")
def test_process_dataset_resume_skips_hot_path(
    mock_extract,
    mock_rcsb,
    mock_blast,
    fake_csv,
    tmp_path,
):
    inp, rc = fake_csv

    ctx = MagicMock()
    ctx.ncbi_taxonomy_id = None
    ctx.uniprot_accessions = []
    ctx.preferred_uniprot_accession = None
    ctx.entity_id = None
    mock_rcsb.return_value = ctx
    mock_blast.return_value = [
        BlastHit(
            subject_accession="Q12923",
            pident=97.0,
            length=118,
            bitscore=230.0,
            qstart=1,
            qend=118,
            sstart=400,
            send=517,
            subject_title="x",
        ),
    ]

    seq_calls = iter(
        [
            ("M" * 120, "bmrb_aggregate", False),
            ("M" * 118, "bmrb_aggregate", False),
        ]
    )

    def _extract_once(*args, **kwargs):
        return next(seq_calls)

    mock_extract.side_effect = _extract_once

    def _fake_fasta(acc):
        assert acc == "Q12923"
        return "M" * 500, 500

    with patch(
        "scripts.domain_full_length_bifurcation.fetch_fasta_sequence",
        _fake_fasta,
    ):
        process_dataset(
            input_csv=inp,
            receptor_csv=rc,
            out_dir=tmp_path,
            basename="resume",
            cs_dir=tmp_path,
            pdb_cache_dir=str(tmp_path),
            fetch_bmrb=False,
            blast_db="stub",
            blast_remote=False,
            min_pident=20.0,
            min_query_cov=0.4,
            ratio=1.5,
            sleep_between_uniprot=0.0,
            dry_run_rest=False,
        )

    blast_first = mock_blast.call_count

    mock_extract.side_effect = RuntimeError("_extract_sequence should not run when row is cached")
    mock_rcsb.side_effect = RuntimeError("RCS should not run when row is cached")
    mock_blast.side_effect = RuntimeError("blast should not run when row is cached")

    with patch(
        "scripts.domain_full_length_bifurcation.fetch_fasta_sequence",
        _fake_fasta,
    ):
        process_dataset(
            input_csv=inp,
            receptor_csv=rc,
            out_dir=tmp_path,
            basename="resume",
            cs_dir=tmp_path,
            pdb_cache_dir=str(tmp_path),
            fetch_bmrb=False,
            blast_db="stub",
            blast_remote=False,
            min_pident=20.0,
            min_query_cov=0.4,
            ratio=1.5,
            sleep_between_uniprot=0.0,
            dry_run_rest=False,
            force_recompute=False,
        )

    assert blast_first == mock_blast.call_count


@patch("scripts.domain_full_length_bifurcation.run_blastp")
@patch("scripts.domain_full_length_bifurcation.fetch_holo_polymer_context")
@patch("scripts.domain_full_length_bifurcation._extract_sequence_bmrb_then_pdb")
def test_process_dataset_force_recompute_reruns_compute(
    mock_extract,
    mock_rcsb,
    mock_blast,
    fake_csv,
    tmp_path,
):
    inp, rc = fake_csv
    mock_extract.side_effect = [
        ("M" * 120, "bmrb_aggregate", False),
        ("M" * 118, "bmrb_aggregate", False),
    ]
    ctx = MagicMock()
    ctx.ncbi_taxonomy_id = None
    ctx.uniprot_accessions = []
    ctx.preferred_uniprot_accession = None
    ctx.entity_id = None
    mock_rcsb.return_value = ctx
    mock_blast.return_value = [
        BlastHit(
            subject_accession="Q12923",
            pident=97.0,
            length=118,
            bitscore=230.0,
            qstart=1,
            qend=118,
            sstart=400,
            send=517,
            subject_title="x",
        ),
    ]

    def _fake_fasta(acc):
        return "M" * 500, 500

    with patch(
        "scripts.domain_full_length_bifurcation.fetch_fasta_sequence",
        _fake_fasta,
    ):
        process_dataset(
            input_csv=inp,
            receptor_csv=rc,
            out_dir=tmp_path,
            basename="force",
            cs_dir=tmp_path,
            pdb_cache_dir=str(tmp_path),
            fetch_bmrb=False,
            blast_db="stub",
            blast_remote=False,
            min_pident=20.0,
            min_query_cov=0.4,
            ratio=1.5,
            sleep_between_uniprot=0.0,
            dry_run_rest=False,
        )

    mock_blast.reset_mock()
    mock_extract.side_effect = [
        ("M" * 120, "bmrb_aggregate", False),
        ("M" * 118, "bmrb_aggregate", False),
    ]

    with patch(
        "scripts.domain_full_length_bifurcation.fetch_fasta_sequence",
        _fake_fasta,
    ):
        process_dataset(
            input_csv=inp,
            receptor_csv=rc,
            out_dir=tmp_path,
            basename="force",
            cs_dir=tmp_path,
            pdb_cache_dir=str(tmp_path),
            fetch_bmrb=False,
            blast_db="stub",
            blast_remote=False,
            min_pident=20.0,
            min_query_cov=0.4,
            ratio=1.5,
            sleep_between_uniprot=0.0,
            dry_run_rest=False,
            force_recompute=True,
        )

    assert mock_blast.call_count == 1
