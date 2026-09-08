#!/usr/bin/env python3
"""
Propagate apo/holo pair edits from CSP_UBQ_ph0.5_temp5C.csv into other data CSVs.

ph0.5 is the source of truth for pair edits (especially apo_bmrb / apo_pdb swaps).
Matching uses (holo_pdb, holo_bmrb) when that key is unique in ph0.5. Multi-apo
holos in ph0.5 are preserved as multi-row sets (no collapse).

Does not overwrite CSP_UBQ_ph0.5_temp5C.csv.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

PAIR_FIELDS = ("apo_bmrb", "apo_pdb", "holo_bmrb", "holo_pdb")
SKIP_FILES = frozenset({"CSP_UBQ_ph0.5_temp5C.csv"})


def _norm(v: str | None) -> str:
    return (v or "").strip()


def _key_hh(row: dict[str, str]) -> tuple[str, str]:
    return (_norm(row.get("holo_pdb")), _norm(row.get("holo_bmrb")))


def _key_ha(row: dict[str, str]) -> tuple[str, str]:
    return (_norm(row.get("holo_pdb")), _norm(row.get("apo_bmrb")))


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        fieldnames = list(reader.fieldnames)
        rows = [{k: (row.get(k) or "") for k in fieldnames} for row in reader]
    return fieldnames, rows


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _build_ph05_maps(
    ph_rows: list[dict[str, str]],
) -> tuple[dict[tuple[str, str], dict[str, str]], dict[tuple[str, str], list[dict[str, str]]]]:
    by_hh: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in ph_rows:
        by_hh[_key_hh(row)].append(row)

    unique_map: dict[tuple[str, str], dict[str, str]] = {}
    multi: dict[tuple[str, str], list[dict[str, str]]] = {}
    for key, group in by_hh.items():
        if len(group) == 1:
            unique_map[key] = group[0]
        else:
            multi[key] = group
    return unique_map, multi


def _copy_pair_fields(dst: dict[str, str], src: dict[str, str], fields: list[str]) -> list[str]:
    changed: list[str] = []
    for field in fields:
        if field not in dst:
            continue
        new_val = _norm(src.get(field))
        old_val = _norm(dst.get(field))
        if new_val != old_val:
            dst[field] = src.get(field, "")
            changed.append(f"{field}:{old_val}->{new_val}")
    return changed


def sync_file(
    path: Path,
    *,
    unique_map: dict[tuple[str, str], dict[str, str]],
    multi_map: dict[tuple[str, str], list[dict[str, str]]],
    ph_by_ha: dict[tuple[str, str], dict[str, str]],
    dry_run: bool,
) -> tuple[int, list[str], list[str]]:
    """Return (n_updated, change_lines, ambiguous_lines)."""
    fieldnames, rows = _read_csv(path)
    if "holo_pdb" not in fieldnames or "apo_bmrb" not in fieldnames:
        return 0, [], []

    has_holo_bmrb = "holo_bmrb" in fieldnames
    shared_fields = [f for f in PAIR_FIELDS if f in fieldnames]

    n_updated = 0
    changes: list[str] = []
    ambiguous: list[str] = []

    # Precompute which (holo, apo) already exist so unique-key swaps do not
    # create duplicate HA keys when both old and new apos are present.
    ha_counts: dict[tuple[str, str], int] = defaultdict(int)
    for r in rows:
        ha_counts[_key_ha(r)] += 1

    for idx, row in enumerate(rows):
        if not has_holo_bmrb:
            ha = _key_ha(row)
            if ha in ph_by_ha:
                changed = _copy_pair_fields(
                    row, ph_by_ha[ha], [f for f in shared_fields if f != "apo_bmrb"]
                )
                if changed:
                    n_updated += 1
                    changes.append(f"{path.name}:row{idx+2} ha-refresh " + ", ".join(changed))
            continue

        hh = _key_hh(row)
        ha = _key_ha(row)

        if hh in unique_map:
            src = unique_map[hh]
            desired_apo = _norm(src.get("apo_bmrb"))
            current_apo = _norm(row.get("apo_bmrb"))
            desired_ha = (_norm(src.get("holo_pdb")), desired_apo)

            if current_apo == desired_apo:
                changed = _copy_pair_fields(
                    row, src, [f for f in shared_fields if f != "apo_bmrb"]
                )
                if changed:
                    n_updated += 1
                    changes.append(
                        f"{path.name}:row{idx+2} holo={hh[0]} refresh " + ", ".join(changed)
                    )
                continue

            # Stale apo for this unique complex.
            if ha_counts.get(desired_ha, 0) > 0:
                # Desired pair already present — mark this stale row for removal.
                row["__delete__"] = "1"
                ha_counts[ha] -= 1
                n_updated += 1
                changes.append(
                    f"{path.name}:row{idx+2} holo={hh[0]} holo_bmrb={hh[1]} "
                    f"remove stale_apo={current_apo} (desired_apo={desired_apo} already present)"
                )
                continue

            old_ha = ha
            changed = _copy_pair_fields(row, src, shared_fields)
            if changed:
                ha_counts[old_ha] -= 1
                ha_counts[desired_ha] += 1
                n_updated += 1
                changes.append(
                    f"{path.name}:row{idx+2} holo={hh[0]} holo_bmrb={hh[1]} "
                    + ", ".join(changed)
                )
            continue

        if hh in multi_map:
            ph_apos = {_norm(r.get("apo_bmrb")) for r in multi_map[hh]}
            if ha in ph_by_ha:
                changed = _copy_pair_fields(
                    row, ph_by_ha[ha], [f for f in shared_fields if f != "apo_bmrb"]
                )
                if changed:
                    n_updated += 1
                    changes.append(
                        f"{path.name}:row{idx+2} multi-apo refresh holo={hh[0]} "
                        + ", ".join(changed)
                    )
            elif current_apo := _norm(row.get("apo_bmrb")):
                if current_apo not in ph_apos:
                    ambiguous.append(
                        f"{path.name}:row{idx+2} holo={hh[0]} holo_bmrb={hh[1]} "
                        f"apo={current_apo} not in ph0.5 multi-apo set {sorted(ph_apos)}"
                    )
            continue

        # Unrelated to ph0.5 hh keys: optional ha refresh if exact pair exists.
        if ha in ph_by_ha:
            changed = _copy_pair_fields(
                row, ph_by_ha[ha], [f for f in shared_fields if f != "apo_bmrb"]
            )
            if changed:
                n_updated += 1
                changes.append(f"{path.name}:row{idx+2} ha-refresh " + ", ".join(changed))

    rows = [r for r in rows if r.get("__delete__") != "1"]
    for r in rows:
        r.pop("__delete__", None)

    # Ensure every ph0.5 pair exists in CSP_UBQ.csv.
    if path.name == "CSP_UBQ.csv":
        present_ha = {_key_ha(r) for r in rows}
        for ha, src in ph_by_ha.items():
            if ha in present_ha:
                continue
            hh = _key_hh(src)
            converted = False
            if hh in unique_map:
                for row in rows:
                    if _key_hh(row) == hh and _key_ha(row) != ha:
                        if ha in {_key_ha(r) for r in rows}:
                            break
                        old_apo = _norm(row.get("apo_bmrb"))
                        _copy_pair_fields(row, src, shared_fields)
                        n_updated += 1
                        converted = True
                        present_ha.add(ha)
                        changes.append(
                            f"{path.name}: convert stale holo={hh[0]} apo {old_apo}->{ha[1]}"
                        )
                        break
            if not converted and ha not in {_key_ha(r) for r in rows}:
                new_row = {k: "" for k in fieldnames}
                for k in fieldnames:
                    if k in src:
                        new_row[k] = src[k]
                rows.append(new_row)
                n_updated += 1
                changes.append(f"{path.name}: append missing ph0.5 pair holo={ha[0]} apo={ha[1]}")

    if n_updated and not dry_run:
        _write_csv(path, fieldnames, rows)

    return n_updated, changes, ambiguous


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--ph05",
        type=Path,
        default=_REPO / "data" / "CSP_UBQ_ph0.5_temp5C.csv",
    )
    ap.add_argument(
        "--data-dir",
        type=Path,
        default=_REPO / "data",
    )
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    ph05 = args.ph05 if args.ph05.is_absolute() else _REPO / args.ph05
    data_dir = args.data_dir if args.data_dir.is_absolute() else _REPO / args.data_dir

    if not ph05.is_file():
        print(f"Missing ph0.5 CSV: {ph05}", file=sys.stderr)
        return 1

    _, ph_rows = _read_csv(ph05)
    unique_map, multi_map = _build_ph05_maps(ph_rows)
    ph_by_ha = {_key_ha(r): r for r in ph_rows}

    print(f"ph0.5 rows: {len(ph_rows)}")
    print(f"unique (holo_pdb, holo_bmrb) sync keys: {len(unique_map)}")
    print(f"multi-apo (holo_pdb, holo_bmrb) groups: {len(multi_map)}")
    for key, group in sorted(multi_map.items()):
        apos = sorted(_norm(r.get("apo_bmrb")) for r in group)
        print(f"  multi {key[0]}/{key[1]} apos={apos}")

    total_updated = 0
    all_changes: list[str] = []
    all_ambiguous: list[str] = []
    files_touched: list[str] = []

    csv_paths = sorted(data_dir.glob("*.csv"), key=lambda p: (p.name != "CSP_UBQ.csv", p.name))
    for path in csv_paths:
        if path.name in SKIP_FILES:
            continue
        n_updated, changes, ambiguous = sync_file(
            path,
            unique_map=unique_map,
            multi_map=multi_map,
            ph_by_ha=ph_by_ha,
            dry_run=args.dry_run,
        )
        if n_updated or ambiguous:
            files_touched.append(
                f"{path.name} (updated_rows={n_updated}, ambiguous={len(ambiguous)})"
            )
        total_updated += n_updated
        all_changes.extend(changes)
        all_ambiguous.extend(ambiguous)

    print("\n=== Audit ===")
    print(f"mode: {'dry-run' if args.dry_run else 'write'}")
    print(f"files with updates/ambiguities: {len(files_touched)}")
    for line in files_touched:
        print(f"  {line}")
    print(f"total row updates: {total_updated}")
    if all_changes:
        print("\nChanges:")
        for line in all_changes:
            print(f"  {line}")
    if all_ambiguous:
        print("\nAmbiguous (manual review):")
        for line in all_ambiguous:
            print(f"  {line}")
    else:
        print("\nAmbiguous: none")

    # Re-read CSP_UBQ after write (or from memory path).
    csp_path = data_dir / "CSP_UBQ.csv"
    _, csp_rows = _read_csv(csp_path)
    if args.dry_run:
        # Dry-run did not write; re-simulate gate on in-memory is harder — skip strict gate.
        print("\nDry-run complete (gate skipped).")
        return 0

    csp_ha = {_key_ha(r) for r in csp_rows}
    missing = [ha for ha in ph_by_ha if ha not in csp_ha]
    if missing:
        print(f"\nERROR: {len(missing)} ph0.5 pairs still missing from CSP_UBQ.csv:", file=sys.stderr)
        for ha in missing[:20]:
            print(f"  holo={ha[0]} apo={ha[1]}", file=sys.stderr)
        return 1

    remaining = 0
    csp_by_hh: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for r in csp_rows:
        csp_by_hh[_key_hh(r)].append(r)
    for hh, src in unique_map.items():
        desired = _norm(src.get("apo_bmrb"))
        group = csp_by_hh.get(hh, [])
        if not group:
            remaining += 1
            continue
        if not any(_norm(r.get("apo_bmrb")) == desired for r in group):
            remaining += 1

    print(f"\nUnique-key holos missing desired apo in CSP_UBQ.csv: {remaining}")
    if remaining:
        print("ERROR: sync incomplete for unique (holo_pdb, holo_bmrb) keys.", file=sys.stderr)
        return 1

    print("\nGate OK: all ph0.5 pairs present in CSP_UBQ.csv; unique-key apo swaps applied.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
