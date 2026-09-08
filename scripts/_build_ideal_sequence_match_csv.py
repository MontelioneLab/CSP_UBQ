#!/usr/bin/env python3
"""Build ph0.5-filtered ideal-sequence-match cohort CSV from the user pair list."""

from __future__ import annotations

import csv
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent

# apo_bmrb, holo_bmrb, holo_pdb (from user table; order preserved, duplicates kept once)
USER_PAIRS = [
    ("34688", "4516", "1d5g"),
    ("6268", "5987", "1r8u"),
    ("51289", "5480", "1sy9"),
    ("5131", "6060", "1vj6"),
    ("6261", "6456", "1wa7"),
    ("6592", "7241", "2hug"),
    ("15503", "15504", "2jw1"),
    ("51289", "15624", "2jzi"),
    ("5555", "17344", "2l75"),
    ("51289", "17807", "2lgf"),
    ("17813", "17808", "2lgg"),
    ("17813", "17812", "2lgk"),
    ("18455", "18434", "2lsk"),
    ("51289", "17771", "2m0j"),
    ("51289", "17771", "2m0k"),
    ("51289", "19238", "2mes"),
    ("17769", "25230", "2mur"),
    ("18455", "25559", "2n1g"),
    ("5958", "25919", "2n9x"),
    ("51289", "30063", "5j8h"),
    ("51289", "34161", "5oeo"),
    ("17769", "30276", "5vf0"),
    ("18827", "34307", "6h8c"),
    ("28070", "34384", "6r5g"),
    ("30782", "30791", "7jyz"),
    ("30782", "30786", "7jq8"),
    ("30782", "30790", "7jyn"),
    ("52209", "25694", "2n55"),
    ("18824", "18825", "2m0u"),
    ("17128", "17127", "2l29"),
    ("15710", "11040", "2rol"),
    ("19757", "34224", "6fdt"),
    ("5318", "5310", "2pld"),
    ("17627", "4453", "1h8b"),
    ("25866", "25865", "2n8t"),
    ("50765", "5532", "1m4p"),
    ("4232", "5738", "1npq"),
    ("6558", "6559", "1ywi"),
    ("5659", "6823", "2b0f"),
    ("15670", "15671", "2k17"),
    ("5958", "15877", "2k6q"),
    ("16066", "16065", "2kc8"),
    ("15972", "16967", "2kyl"),
    ("17018", "17019", "2kzu"),
    ("17254", "17255", "2l4t"),
    ("15001", "17569", "2lbm"),
    ("16396", "17702", "2le8"),
    ("16835", "17879", "2li5"),
    ("10318", "17364", "2loz"),
    ("19428", "19429", "2mc6"),
    ("19913", "19914", "2mnz"),
    ("10235", "10236", "2roz"),
    ("11362", "11115", "2rr4"),
    ("7141", "11508", "2rsy"),
    ("4901", "36013", "5gow"),
    ("34068", "34067", "5mf9"),
    ("25544", "30500", "6e5n"),
    ("34247", "34248", "6g04"),
    ("30808", "30809", "7klr"),
    ("30926", "30950", "7s5j"),
    ("27579", "31080", "8sg2"),
    ("25229", "25230", "2mur"),
    ("15524", "15525", "2kwi"),
    ("11532", "18709", "2ly4"),
    ("18581", "18582", "2lvo"),
    ("6457", "18737", "2lz6"),
    ("6457", "25070", "2mre"),
    ("25883", "25829", "2n80"),
    ("25883", "25833", "2n83"),
    ("30092", "30093", "5jyv"),
    ("30599", "30605", "6oqj"),
    ("34724", "34727", "7zey"),
    ("16983", "34727", "7zey"),
    ("15670", "34986", "9qlm"),
    ("15125", "26041", "2ncz"),
    ("15125", "26042", "2nd0"),
    ("15125", "26043", "2nd1"),
    ("30782", "30367", "6bgg"),
    ("15125", "30373", "6bnh"),
    ("50148", "17270", "2l5e"),
    ("18250", "18238", "2lp8"),
]


def main() -> None:
    ph_path = _REPO / "data" / "CSP_UBQ_ph0.5_temp5C.csv"
    with ph_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        ph = {
            (
                r["apo_bmrb"].strip(),
                r["holo_bmrb"].strip(),
                r["holo_pdb"].strip().upper(),
            ): r
            for r in reader
        }

    matched = []
    missing = []
    seen = set()
    for apo, holo_bmrb, holo_pdb in USER_PAIRS:
        key = (apo, holo_bmrb, holo_pdb.upper())
        if key in seen:
            continue
        seen.add(key)
        if key in ph:
            matched.append(ph[key])
        else:
            missing.append(key)

    out = _REPO / "data" / "CSP_UBQ_ph0.5_temp5C_ideal_sequence_match.csv"
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(matched)

    miss_path = _REPO / "data" / "ideal_sequence_match_not_in_ph05.csv"
    with miss_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["apo_bmrb", "holo_bmrb", "holo_pdb"])
        writer.writerows(missing)

    print(f"wrote {len(matched)} rows -> {out}")
    print(f"not in ph0.5: {len(missing)} -> {miss_path}")


if __name__ == "__main__":
    main()
