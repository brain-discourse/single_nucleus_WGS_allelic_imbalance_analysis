#!/usr/bin/env python3
# ==========================================================
# author: Meethila Gade
# created on: 31-01-2026-15h-53m
# github: https://github.com/brain-discourse
# copyright: 2026
# description: per cell run-length metrics across chr3 using a bulk-defined germline het SNP panel
# ==========================================================
import argparse
import csv
import gzip
import os
import sys
from collections import Counter
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

try:
    import yaml  # type: ignore
except Exception:
    yaml = None

try:
    import matplotlib.pyplot as plt  # type: ignore
except Exception:
    plt = None

BASES = ("A", "C", "G", "T")


@dataclass(frozen=True)
class Params:
    min_cell_dp: int
    min_major_allele_frac: float
    monoallelic_af_cut: float
    window_nsnp: int
    step_nsnp: int


def eprint(*args):
    print(*args, file=sys.stderr)


def open_maybe_gzip(path: str):
    if path and path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def load_yaml_params(path: str) -> Params:
    if yaml is None:
        raise RuntimeError("PyYAML not available. Install pyyaml or provide params differently.")
    with open(path, "r") as f:
        d = yaml.safe_load(f)
    return Params(
        min_cell_dp=int(d.get("min_cell_dp", 6)),
        min_major_allele_frac=float(d.get("min_major_allele_frac", 0.8)),
        monoallelic_af_cut=float(d.get("monoallelic_af_cut", 0.95)),
        window_nsnp=int(d.get("window_nsnp", 50)),
        step_nsnp=int(d.get("step_nsnp", 10)),
    )


def read_cell_ids(samples_tsv: str) -> List[str]:
    cell_ids = []
    with open(samples_tsv, "r", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("type") == "cell":
                cell_ids.append(row["id"])
    if not cell_ids:
        raise ValueError("No cell rows found in samples.tsv")
    return cell_ids


def load_sites_from_vcf(vcf_gz: str) -> List[Tuple[str, int, str, str]]:
    """
    Returns ordered list of (chrom, pos, ref, alt) in VCF order.
    Assumes already restricted to chr3 and filtered.
    """
    sites = []
    with open_maybe_gzip(vcf_gz) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            chrom = p[0]
            pos = int(p[1])
            ref = p[3].upper()
            alt = p[4].split(",")[0].upper()
            if ref in BASES and alt in BASES:
                sites.append((chrom, pos, ref, alt))
    if not sites:
        raise ValueError(f"No SNP sites read from {vcf_gz}")
    return sites


def load_phasing_map(phased_vcf_gz: str, sample_name: Optional[str] = None) -> Dict[Tuple[str, int], Tuple[Optional[int], Optional[str]]]:
    """
    Map site -> (ref_hap, ps)
    Where ref_hap is 1 or 2 if GT phased and determines which haplotype carries REF allele.
      - GT 0|1 => REF on hap1, ALT on hap2 => ref_hap=1
      - GT 1|0 => REF on hap2, ALT on hap1 => ref_hap=2
    ps is the phase set identifier if present (FORMAT:PS), else None.
    If unphased or missing, returns (None, None) for that site.
    """
    if not phased_vcf_gz or not os.path.exists(phased_vcf_gz):
        return {}

    ph_map: Dict[Tuple[str, int], Tuple[Optional[int], Optional[str]]] = {}
    with open_maybe_gzip(phased_vcf_gz) as f:
        sample_idx = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                if not samples:
                    return {}
                if sample_name:
                    if sample_name not in samples:
                        raise ValueError(f"Sample {sample_name} not found in phased VCF samples: {samples[:5]}...")
                    sample_idx = 9 + samples.index(sample_name)
                else:
                    sample_idx = 9
                continue
            p = line.rstrip("\n").split("\t")
            chrom = p[0]
            pos = int(p[1])
            fmt = p[8].split(":")
            if sample_idx is None:
                continue
            sample = p[sample_idx].split(":")
            fmt_map = {k: sample[i] if i < len(sample) else "" for i, k in enumerate(fmt)}
            gt = fmt_map.get("GT", "")
            ps = fmt_map.get("PS", None)
            if "|" not in gt:
                ph_map[(chrom, pos)] = (None, None)
                continue
            a, b = gt.split("|", 1)
            if a == "." or b == ".":
                ph_map[(chrom, pos)] = (None, None)
                continue
            if (a, b) == ("0", "1"):
                ph_map[(chrom, pos)] = (1, ps)
            elif (a, b) == ("1", "0"):
                ph_map[(chrom, pos)] = (2, ps)
            else:
                ph_map[(chrom, pos)] = (None, ps)
    return ph_map


def read_counts_for_cell(counts_path: str) -> Dict[Tuple[str, int], Tuple[int, int, int, int, int]]:
    """Return mapping (chrom,pos) -> (dp, A,C,G,T)"""
    if not os.path.exists(counts_path):
        raise FileNotFoundError(f"Missing counts file: {counts_path}")
    m = {}
    with open_maybe_gzip(counts_path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("chrom\tpos\tref\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 8:
                continue
            chrom, pos_s, ref_mp, dp_s, a, c, g, t = parts
            try:
                pos = int(pos_s)
                m[(chrom, pos)] = (int(dp_s), int(a), int(c), int(g), int(t))
            except ValueError:
                continue
    return m


def monoallelic_and_direction(dp: int, refc: int, altc: int, params: Params) -> Tuple[bool, Optional[str], float]:
    """
    Return (is_mono, direction_label, af)
    direction_label is REF_ONLY or ALT_ONLY if monoallelic else None.
    """
    if dp <= 0:
        return (False, None, float("nan"))
    af = altc / dp
    if dp < params.min_cell_dp:
        return (False, None, af)
    if af <= 1.0 - params.monoallelic_af_cut:
        return (True, "REF_ONLY", af)
    if af >= params.monoallelic_af_cut:
        return (True, "ALT_ONLY", af)
    return (False, None, af)


def direction_with_phasing(ref_only_or_alt_only: str, ref_hap: Optional[int]) -> str:
    """
    Convert REF_ONLY/ALT_ONLY into HAP1/HAP2 if ref_hap is known.
    """
    if ref_hap not in (1, 2):
        return ref_only_or_alt_only
    if ref_hap == 1:
        return "HAP1" if ref_only_or_alt_only == "REF_ONLY" else "HAP2"
    else:
        return "HAP2" if ref_only_or_alt_only == "REF_ONLY" else "HAP1"


def max_run_same_direction(states: List[str]) -> int:
    """
    Define runs as consecutive entries of the same monoallelic direction.
    BIALLELIC breaks runs.
    """
    best = 0
    cur_dir = None
    cur_len = 0
    for s in states:
        if s in ("BIALLELIC", "NA"):
            cur_dir = None
            cur_len = 0
            continue
        if cur_dir is None or s != cur_dir:
            cur_dir = s
            cur_len = 1
        else:
            cur_len += 1
        if cur_len > best:
            best = cur_len
    return best


def dir_site_states(states: List[str]) -> List[str]:
    """Return only monoallelic direction states (exclude BIALLELIC/NA)."""
    return [s for s in states if s not in ("BIALLELIC", "NA")]


def mode_and_stats(states: List[str]) -> Tuple[str, int, float, str]:
    """
    Among monoallelic direction states, compute:
      - run_dir_mode: direction label or NA/TIE
      - run_dir_mode_count: count for the mode (0 if NA)
      - run_dir_mode_frac: fraction among direction sites (0.0 if none)
      - top_run_dir: one of {HAP1,HAP2,REF_ONLY,ALT_ONLY} else NA
    """
    dirs = dir_site_states(states)
    n_dir = len(dirs)
    if n_dir == 0:
        return ("NA", 0, 0.0, "NA")

    c = Counter(dirs)
    top = c.most_common()
    if len(top) == 1:
        label, count = top[0]
        frac = count / n_dir
        return (label, count, frac, label)

    if top[0][1] == top[1][1]:
        max_count = top[0][1]
        frac = max_count / n_dir
        return ("TIE", max_count, frac, "NA")

    label, count = top[0]
    frac = count / n_dir
    return (label, count, frac, label)


def window_iter(n: int, w: int, step: int):
    i = 0
    while i + w <= n:
        yield i, i + w
        i += step


def write_tsv(path: str, rows: List[dict], fieldnames: List[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_track_gz(path: str, rows: List[dict], fieldnames: List[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def plot_track(cell_id: str, track_rows: List[dict], out_dir: str) -> None:
    if plt is None or not track_rows:
        return
    xs = [int(r["window_mid_pos"]) for r in track_rows]
    ys = [int(r["max_run_len_same_direction"]) for r in track_rows]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xs, ys)
    ax.set_title(f"Max monoallelic run length: {cell_id}")
    ax.set_xlabel("Window midpoint (bp)")
    ax.set_ylabel("Max run length (SNPs)")
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, f"max_run_len_vs_position_{cell_id}.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--hets_vcf", required=True)
    ap.add_argument("--phased_vcf", default="", help="Optional phased VCF.gz for same sites")
    ap.add_argument("--phased_sample", default="", help="Optional sample name in phased VCF")
    ap.add_argument("--counts_dir", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--out_table", required=True)
    ap.add_argument("--out_tracks_dir", required=True)
    ap.add_argument("--out_plots", required=True)
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()

    params = load_yaml_params(args.params)
    cell_ids = read_cell_ids(args.samples)
    sites = load_sites_from_vcf(args.hets_vcf)
    n_sites_in_vcf = len(sites)

    phased_map = {}
    phased_sample_used = ""
    phasing_mode = "unphased"

    if args.phased_vcf and os.path.exists(args.phased_vcf):
        phased_map = load_phasing_map(
            args.phased_vcf,
            sample_name=(args.phased_sample or None)
        )
        phasing_mode = "phased"
        phased_sample_used = args.phased_sample or "FIRST_SAMPLE"
        eprint(f"Loaded phasing for {len(phased_map)} sites from {args.phased_vcf} sample={phased_sample_used}")
    else:
        eprint("No phased VCF provided/found; using unphased REF_ONLY/ALT_ONLY directions.")
        
    eprint(f"Total SNP sites to analyze: {n_sites_in_vcf}")
    all_rows = []

    site_keys = [(c, p) for (c, p, _, _) in sites]
    site_pos = [p for (_, p, _, _) in sites]

    for cell_id in cell_ids:
        counts_path = os.path.join(args.counts_dir, f"{cell_id}.counts.tsv.gz")
        counts = read_counts_for_cell(counts_path)
        n_sites_in_counts_map = len(counts)
        n_sites_overlap = 0
        
        dp_list = []
        state_list = []
        usable_list = []
        mono_list = []

        for (chrom, pos, ref, alt) in sites:
            key = (chrom, pos)
            if key not in counts:
                dp_list.append(0)
                state_list.append("NA")
                usable_list.append(False)
                mono_list.append(False)
                continue
                
            n_sites_overlap += 1
            dp, a, c, g, t = counts[key]
            base_counts = {"A": a, "C": c, "G": g, "T": t}
            refc = base_counts.get(ref, 0)
            altc = base_counts.get(alt, 0)

            is_mono, _, af = monoallelic_and_direction(dp, refc, altc, params)
            usable = dp >= params.min_cell_dp
            if not usable:
                dp_list.append(dp)
                state_list.append("NA")
                usable_list.append(False)
                mono_list.append(False)
                continue

            dir_basic = None
            if af <= 1.0 - params.min_major_allele_frac:
                dir_basic = "REF_ONLY"
            elif af >= params.min_major_allele_frac:
                dir_basic = "ALT_ONLY"

            if dir_basic is not None:
                ref_hap, ps = phased_map.get(key, (None, None))
                dir_final = direction_with_phasing(dir_basic, ref_hap)
                dp_list.append(dp)
                state_list.append(dir_final)
                usable_list.append(True)
                mono_list.append(is_mono)
            else:
                dp_list.append(dp)
                state_list.append("BIALLELIC")
                usable_list.append(True)
                mono_list.append(False)

        track_rows = []
        win_index = 0
        for i0, i1 in window_iter(len(sites), params.window_nsnp, params.step_nsnp):
            start_pos = site_pos[i0]
            end_pos = site_pos[i1 - 1]
            mid_pos = int((start_pos + end_pos) / 2)

            usable_states = []
            n_usable = 0
            n_mono = 0

            for j in range(i0, i1):
                if usable_list[j]:
                    n_usable += 1
                    usable_states.append(state_list[j])
                    if mono_list[j]:
                        n_mono += 1

            frac_not_usable = 1.0 - (n_usable / (i1 - i0))
            frac_mono = (n_mono / n_usable) if n_usable > 0 else 0.0
            n_dir = len(dir_site_states(usable_states))
            
            if n_dir < 10:
                max_run = 0
                run_dir_mode, run_dir_mode_count, run_dir_mode_frac, top_run_dir = ("NA", 0, 0.0, "NA")
            else:
                max_run = max_run_same_direction(usable_states)
                run_dir_mode, run_dir_mode_count, run_dir_mode_frac, top_run_dir = mode_and_stats(usable_states)

            row = {
                "cell_id": cell_id,
                "chrom": sites[i0][0],
                "window_index": win_index,
                "run_dir_mode": run_dir_mode,
                "run_dir_mode_count": run_dir_mode_count,
                "run_dir_mode_frac": f"{run_dir_mode_frac:.6f}",
                "top_run_dir": top_run_dir,
                "start_pos": start_pos,
                "end_pos": end_pos,
                "window_mid_pos": mid_pos,
                "min_cell_dp": params.min_cell_dp,
                "monoallelic_af_cut": params.monoallelic_af_cut,
                "window_nsnp": params.window_nsnp,
                "step_nsnp": params.step_nsnp,
                "phasing_mode": phasing_mode,
                "phased_sample": phased_sample_used,
                "n_sites_in_vcf": n_sites_in_vcf,
                "n_sites_in_counts_map": n_sites_in_counts_map,
                "n_sites_overlap": n_sites_overlap,
                "n_sites_total": (i1 - i0),
                "n_sites_usable": n_usable,
                "n_monoallelic": n_mono,
                "frac_monoallelic": f"{frac_mono:.6f}",
                "frac_not_usable": f"{frac_not_usable:.6f}",
                "max_run_len_same_direction": max_run,
            }

            all_rows.append(row)
            track_rows.append(row)
            win_index += 1

        track_path = os.path.join(args.out_tracks_dir, f"{cell_id}.track.tsv.gz")
        track_fields = [
            "cell_id", "chrom", "window_index", "start_pos", "end_pos", "window_mid_pos",
            "min_cell_dp", "monoallelic_af_cut", "window_nsnp", "step_nsnp",
            "phasing_mode", "phased_sample",
            "n_sites_in_vcf", "n_sites_in_counts_map", "n_sites_overlap",
            "n_sites_total", "n_sites_usable", "n_monoallelic",
            "frac_monoallelic", "frac_not_usable", "max_run_len_same_direction",
            "run_dir_mode", "run_dir_mode_count", "run_dir_mode_frac", "top_run_dir",
        ]

        write_track_gz(track_path, track_rows, track_fields)

        if (not args.no_plots) and plt is not None:
            plot_track(cell_id, track_rows, args.out_plots)

        eprint(f"Wrote track: {track_path}")

    write_tsv(args.out_table, all_rows, track_fields)
    eprint(f"Wrote: {args.out_table}")


if __name__ == "__main__":
    main()