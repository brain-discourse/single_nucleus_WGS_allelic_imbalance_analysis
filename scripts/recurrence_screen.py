#!/usr/bin/env python3
# ==========================================================
# author: Meethila Gade
# created on: 31-01-2026-15h-54m
# github: https://github.com/brain-discourse
# copyright: 2026
# description: Screen contiguity windows for recurrence across cells
# ==========================================================
import argparse
import csv
import os
import sys
from collections import Counter, defaultdict
from typing import Any, Dict, List, Tuple

try:
    import yaml  # type: ignore
except Exception:
    yaml = None

try:
    import matplotlib.pyplot as plt  # type: ignore
except Exception:
    plt = None


def eprint(*args):
    print(*args, file=sys.stderr)


def load_yaml(path: str) -> dict:
    if yaml is None:
        return {}
    with open(path, "r") as f:
        return yaml.safe_load(f) or {}


def read_contiguity_table(path: str) -> List[dict]:
    with open(path, "r", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = [row for row in r]
    if not rows:
        raise ValueError(f"No rows in {path}")
    return rows


def write_tsv(path: str, rows: List[dict], fieldnames: List[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def plot_cells_passing(windows: List[dict], out_dir: str) -> None:
    if plt is None or not windows:
        return
    xs = [int(w["window_mid_pos"]) for w in windows if w.get("window_mid_pos", "") != ""]
    ys = [int(w["n_cells_passing_basic"]) for w in windows if w.get("window_mid_pos", "") != ""]
    if not xs:
        return
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xs, ys)
    ax.set_title("Direction-consistent cells passing contiguity thresholds by position")
    ax.set_xlabel("Window midpoint (bp)")
    ax.set_ylabel("# cells passing(dir-consistent)")
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, "cells_passing_by_position.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def unique_mode(labels: List[str]) -> Tuple[str, bool]:
    """Return (mode_label, is_unique). If tie or empty => ("NA", False)."""
    if not labels:
        return ("NA", False)
    c = Counter(labels).most_common()
    if len(c) > 1 and c[0][1] == c[1][1]:
        return ("NA", False)
    return (c[0][0], True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contiguity_table", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--out_table", required=True)
    ap.add_argument("--out_plots", required=True)
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()

    p = load_yaml(args.params)
    min_run_len = int(p.get("min_run_len", 25))
    min_cells_recur = int(p.get("min_cells_recur", 2))
    min_usable_fraction = float(p.get("min_usable_fraction", 0.6))
    require_same_direction = bool(p.get("require_same_direction", False))
    min_dir_mode_frac = float(p.get("min_dir_mode_frac", 0.0))
    rows = read_contiguity_table(args.contiguity_table)

    grouped: Dict[Tuple[str, str, str, str], List[dict]] = defaultdict(list)
    for r in rows:
        key = (r["chrom"], r["start_pos"], r["end_pos"], r["window_index"])
        grouped[key].append(r)
        
    out_rows: List[Dict[str, Any]] = []

    for (chrom, start_pos, end_pos, win_idx), lst in grouped.items():
        mid_pos = int(lst[0]["window_mid_pos"]) if lst and lst[0].get("window_mid_pos") else ""

        passing_basic: List[str] = []
        passing_dir_pairs: List[Tuple[str, str]] = []

        for r in lst:
            n_total = int(r["n_sites_total"])
            n_usable = int(r["n_sites_usable"])
            max_run = int(r["max_run_len_same_direction"])
            if n_total <= 0:
                continue

            usable_frac = n_usable / n_total
            pass_basic = (usable_frac >= min_usable_fraction and max_run >= min_run_len)
            if not pass_basic:
                continue

            cid = r["cell_id"]
            passing_basic.append(cid)

            top_dir = r.get("top_run_dir", "NA")
            dir_frac = float(r.get("run_dir_mode_frac", "0.0"))
            dir_count = int(r.get("run_dir_mode_count", 0))
            
            pass_dir = (
                top_dir != "NA"
                and dir_frac >= min_dir_mode_frac
                and dir_count >= 10
            )

            if pass_dir:
                passing_dir_pairs.append((cid, top_dir))

        n_basic = len(set(passing_basic))

        if not require_same_direction:
            n_pass = n_basic
            if n_pass >= min_cells_recur:
                out_rows.append({
                    "chrom": chrom,
                    "start_pos": int(start_pos),
                    "end_pos": int(end_pos),
                    "window_index": int(win_idx),
                    "window_mid_pos": mid_pos,
                    "n_cells_passing_basic": n_basic,
                    "n_cells_passing_dir": n_basic,
                    "window_top_run_dir": "",
                    "cell_ids": ",".join(sorted(set(passing_basic))),
                    "require_same_direction": str(require_same_direction),
                    "threshold_min_run_len": min_run_len,
                    "threshold_min_usable_fraction": f"{min_usable_fraction:.3f}",
                    "threshold_min_cells_recur": min_cells_recur,
                    "threshold_min_dir_mode_frac": f"{min_dir_mode_frac:.3f}",
                })
            continue
    
        dir_labels = [d for (_, d) in passing_dir_pairs if d != "NA"]
        dir_mode, ok = unique_mode(dir_labels)
        if not ok:
            continue

        passing_dir = sorted({cid for (cid, d) in passing_dir_pairs if d == dir_mode})
        n_dir = len(passing_dir)

        if n_dir >= min_cells_recur:
            out_rows.append({
                "chrom": chrom,
                "start_pos": int(start_pos),
                "end_pos": int(end_pos),
                "window_index": int(win_idx),
                "window_mid_pos": mid_pos,
                "n_cells_passing_basic": n_basic,
                "n_cells_passing_dir": n_dir,
                "window_top_run_dir": dir_mode,
                "cell_ids": ",".join(passing_dir),
                "require_same_direction": str(require_same_direction),
                "threshold_min_run_len": min_run_len,
                "threshold_min_usable_fraction": f"{min_usable_fraction:.3f}",
                "threshold_min_cells_recur": min_cells_recur,
                "threshold_min_dir_mode_frac": f"{min_dir_mode_frac:.3f}",
            })

    out_rows.sort(key=lambda x: (x["chrom"], x["start_pos"], x["end_pos"]))

    fields = [
        "chrom", "start_pos", "end_pos", "window_index", "window_mid_pos",
        "n_cells_passing_basic", "n_cells_passing_dir", "window_top_run_dir",
        "cell_ids", "require_same_direction",
        "threshold_min_run_len", "threshold_min_usable_fraction", "threshold_min_cells_recur",
        "threshold_min_dir_mode_frac",
    ]
    write_tsv(args.out_table, out_rows, fields)
    eprint(f"Wrote: {args.out_table} (n={len(out_rows)} windows)")

    if (not args.no_plots) and plt is not None and out_rows:
        plot_cells_passing(out_rows, args.out_plots)


if __name__ == "__main__":
    main()