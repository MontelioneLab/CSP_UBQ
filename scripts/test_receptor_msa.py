"""Unit tests for receptor five-way sequence MSA PNG helpers (offline)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from scripts.receptor_msa import (
    lift_partner_to_master,
    load_bifurcation_row,
    merge_reference_gaps,
    star_msa_strings,
    try_write_receptor_alignment_png,
)


def test_merge_gap_profile_roundtrip():
    center = "AC"
    c1 = "A-C"
    c2 = "AC"
    g1 = merge_reference_gaps(center, [c1, c2])
    assert "-" in g1
    assert g1.replace("-", "").upper() == center.upper()


def test_star_msa_equal_width():
    labeled = [
        ("a", "MKT"),
        ("b", "MK"),
        ("c", "KT"),
    ]
    rows = star_msa_strings(labeled)
    assert len(rows) == 3
    w = len(rows[0][1])
    for _lab, s in rows:
        assert len(s) == w


def test_lift_identity():
    center = "AC"
    c = "A-C"
    p = "A-C"
    m = merge_reference_gaps(center, [c])
    lifted = lift_partner_to_master(m, c, p, center)
    assert lifted == m


def test_load_bifurcation_domains_first(tmp_path: Path):
    dom = tmp_path / "tst_domains.csv"
    dom.write_text(
        "apo_bmrb,holo_bmrb,holo_pdb,bmrb_apo_seq,bmrb_holo_seq,uniprot_seq,uniprot_accession\n"
        "1,2,1abc,A,A,A,P1\n",
        encoding="utf-8",
    )
    row = load_bifurcation_row(tmp_path, "tst", "1", "2", "1abc")
    assert row is not None
    assert row["uniprot_accession"] == "P1"


@patch("scripts.receptor_msa.write_msa_png")
def test_try_write_png(mock_png, tmp_path: Path):
    dom = tmp_path / "tst2_domains.csv"
    dom.write_text(
        "apo_bmrb,holo_bmrb,holo_pdb,receptor_chain_id,bmrb_apo_seq,bmrb_holo_seq,uniprot_seq,uniprot_accession\n"
        "111,222,1xyz,A,ACDEFGHIK,MKTAYIAKQR,QWERTYSEQID,P00999\n",
        encoding="utf-8",
    )
    outp = tmp_path / "msa.png"
    ok, msg = try_write_receptor_alignment_png(
        data_dir=tmp_path,
        bifurcation_basename="tst2",
        apo_bmrb="111",
        holo_bmrb="222",
        holo_pdb="1xyz",
        receptor_chain="A",
        seq_apo_pdb="MKTAYIAKQRMMMM",
        seq_holo_pdb="MKTAYIAKQR",
        out_png=outp,
    )
    assert ok
    mock_png.assert_called_once()
    aligned = mock_png.call_args[0][0]
    assert len(aligned) >= 2
    w = len(aligned[0][1])
    for _l, s in aligned:
        assert len(s) == w

