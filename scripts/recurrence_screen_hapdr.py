#!/usr/bin/env python3
# ==========================================================
# author: Meethila Gade
# created on: 31-01-2026-15h-54m
# github: https://github.com/brain-discourse
# copyright: 2026
# description: Screen contiguity windows for recurrence across cells-collapse to segment-outputs one row per cell per chrom, per segment (segment = consecutive windows where:pass_basic is true;, pass_dir is true, top_run_dir != NA ;and top_run_dir stays the same across the segment)
# ==========================================================
import argparse
import csv
import os
import sys
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

try:
    import yaml  # type: ignore
except Exception:
    yaml = None


def eprint(*args):
    print(*args, file=sys.stderr)


def load_yaml(path: str) -> dict:
    if yaml is None:
        return {}
    with open(path, "r") as f:
        return yaml.safe_load(f) or {}


def read_tsv(path: str) -> List[dict]:
    with open(path, "r", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
    if not rows:
        raise ValueError(f"No rows in {path}")
    return rows


def write_tsv(path: str, rows: List[dict], fieldnames: List[str]) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def pass_window(r, min_usable_fraction, min_run_len, min_dir_mode_frac, min_dir_mode_count):
    n_total = int(r["n_sites_total"])
    if n_total <= 0:
        return False

    n_usable = int(r["n_sites_usable"])
    usable_frac = n_usable / n_total

    max_run = int(r["max_run_len_same_direction"])
    if usable_frac < min_usable_fraction or max_run < min_run_len:
        return False

    top_dir = r.get("top_run_dir", "NA")
    if top_dir == "NA":
        return False

    dir_frac = float(r.get("run_dir_mode_frac", "0.0"))
    dir_count = int(r.get("run_dir_mode_count", 0))
    if dir_frac < min_dir_mode_frac or dir_count < min_dir_mode_count:
        return False

    return True


def collapse_to_segments(
    lst: List[dict],
    min_usable_fraction: float,
    min_run_len: int,
    min_dir_mode_frac: float,
    min_dir_mode_count: int,
    min_segment_windows: int = 2,
) -> List[Dict[str, Any]]:
    """
    Collapse consecutive passing windows with same top_run_dir into segments.
    lst must be for ONE (cell_id, chrom) and sorted by window_index.
    """
    segments: List[Dict[str, Any]] = []
    cur: Optional[Dict[str, Any]] = None
    prev_idx: Optional[int] = None

    for r in lst:
        idx = int(r["window_index"])

        if prev_idx is not None and idx != prev_idx + 1:
            cur = None
        prev_idx = idx

        if not pass_window(r, min_usable_fraction, min_run_len, min_dir_mode_frac, min_dir_mode_count):
            cur = None
            continue

        d = r["top_run_dir"]
        start_pos = int(r["start_pos"])
        end_pos = int(r["end_pos"])
        mid_pos = int(r["window_mid_pos"]) if r.get("window_mid_pos") not in ("", None) else None
        df = float(r.get("run_dir_mode_frac", "0.0"))

        if cur is None or d != cur["segment_dir"]:
            cur = {
                "segment_dir": d,
                "segment_start_window_index": idx,
                "segment_end_window_index": idx,
                "segment_start_pos": start_pos,
                "segment_end_pos": end_pos,
                "segment_n_windows": 1,
                "segment_dir_frac_min": df,
                "segment_dir_frac_sum": df,
                "segment_mid_sum": mid_pos if mid_pos is not None else 0,
                "segment_mid_n": 1 if mid_pos is not None else 0,
            }
            segments.append(cur)
        else:
            cur["segment_end_window_index"] = idx
            cur["segment_end_pos"] = end_pos
            cur["segment_n_windows"] += 1
            cur["segment_dir_frac_min"] = min(cur["segment_dir_frac_min"], df)
            cur["segment_dir_frac_sum"] += df
            if mid_pos is not None:
                cur["segment_mid_sum"] += mid_pos
                cur["segment_mid_n"] += 1

    out: List[Dict[str, Any]] = []
    for s in segments:
        if s["segment_n_windows"] < min_segment_windows:
            continue
        s["segment_dir_frac_mean"] = s["segment_dir_frac_sum"] / s["segment_n_windows"]
        if s["segment_mid_n"] > 0:
            s["segment_mid_pos_mean"] = int(s["segment_mid_sum"] / s["segment_mid_n"])
        else:
            s["segment_mid_pos_mean"] = ""
            
        del s["segment_dir_frac_sum"], s["segment_mid_sum"], s["segment_mid_n"]
        out.append(s)

    return out


def main():
    ap = argparse.ArgumentParser(description="Collapse within-cell contiguity windows into direction-consistent segments")
    ap.add_argument("--contiguity_table", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--out_table", required=True)
    args = ap.parse_args()

    p = load_yaml(args.params)
    min_run_len = int(p.get("min_run_len", 25))
    min_usable_fraction = float(p.get("min_usable_fraction", 0.4))
    min_dir_mode_frac = float(p.get("min_dir_mode_frac", 0.75))
    min_dir_mode_count = int(p.get("min_dir_mode_count", 15))
    min_segment_windows = int(p.get("min_segment_windows", 2))

    rows = read_tsv(args.contiguity_table)

    grouped: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for r in rows:
        grouped[(r["cell_id"], r["chrom"])].append(r)

    out_rows: List[Dict[str, Any]] = []

    for (cid, chrom), lst in grouped.items():
        lst.sort(key=lambda x: int(x["window_index"]))
        segs = collapse_to_segments(
            lst,
            min_usable_fraction=min_usable_fraction,
            min_run_len=min_run_len,
            min_dir_mode_frac=min_dir_mode_frac,
            min_dir_mode_count=min_dir_mode_count,
            min_segment_windows=min_segment_windows,
        )
        for s in segs:
            out_rows.append({
                "cell_id": cid,
                "chrom": chrom,
                **s,
                "threshold_min_run_len": min_run_len,
                "threshold_min_usable_fraction": f"{min_usable_fraction:.3f}",
                "threshold_min_dir_mode_frac": f"{min_dir_mode_frac:.3f}",
                "threshold_min_dir_mode_count": min_dir_mode_count,
                "threshold_min_segment_windows": min_segment_windows,
            })

    out_rows.sort(key=lambda x: (x["cell_id"], x["chrom"], int(x["segment_start_window_index"])))

    fields = [
        "cell_id", "chrom",
        "segment_dir",
        "segment_start_window_index", "segment_end_window_index", "segment_n_windows",
        "segment_start_pos", "segment_end_pos", "segment_mid_pos_mean",
        "segment_dir_frac_mean", "segment_dir_frac_min",
        "threshold_min_run_len", "threshold_min_usable_fraction",
        "threshold_min_dir_mode_frac", "threshold_min_dir_mode_count",
        "threshold_min_segment_windows",
    ]
    write_tsv(args.out_table, out_rows, fields)
    eprint(f"Wrote: {args.out_table} (n={len(out_rows)} segments)")


if __name__ == "__main__":
    main()