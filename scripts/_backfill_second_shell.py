#!/usr/bin/env python3
"""One-off backfill: remerge master_alignment.csv with second-shell columns."""
from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


def remereg_one(tgt: str):
    from scripts.merge_csv import merge_all_csv_files

    master = os.path.join(tgt, "master_alignment.csv")
    merge_all_csv_files(tgt, master)
    has_filter = os.path.exists(os.path.join(tgt, "second_shell_filter.csv"))
    return os.path.basename(tgt), has_filter


def main() -> int:
    outputs = _ROOT / "outputs"
    targets = sorted(
        str(p)
        for p in outputs.iterdir()
        if p.is_dir() and (p / "csp_table.csv").exists()
    )
    # Skip targets that already have a populated second-shell filter
    pending = []
    for t in targets:
        filt = Path(t) / "second_shell_filter.csv"
        if filt.exists() and filt.stat().st_size > 50:
            continue
        pending.append(t)
    print(f"Found {len(targets)} targets; pending={len(pending)}", flush=True)
    if not pending:
        return 0

    ok = fail = missing_filter = 0
    failures = []
    workers = min(4, os.cpu_count() or 2)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(remereg_one, t): t for t in pending}
        for i, fut in enumerate(as_completed(futs), 1):
            try:
                name, has_filter = fut.result()
                ok += 1
                if not has_filter:
                    missing_filter += 1
                if i % 10 == 0 or i == len(pending):
                    print(
                        f"[{i}/{len(pending)}] ok={ok} fail={fail} "
                        f"missing_filter={missing_filter} last={name}",
                        flush=True,
                    )
            except Exception as e:
                fail += 1
                failures.append((os.path.basename(futs[fut]), str(e)))
                print(f"FAIL {os.path.basename(futs[fut])}: {e}", flush=True)

    print(f"Done. ok={ok} fail={fail} missing_filter={missing_filter}", flush=True)
    for name, err in failures[:30]:
        print(f"  {name}: {err}")
    return 1 if fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
