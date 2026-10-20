#!/usr/bin/env python3
# ==========================================================
# author: Meethila Gade
# created on: 31-01-2026-15h-51m
# github: https://github.com/brain-discourse
# copyright: 2026
# description: Compute cov binned monoallelic rates per cell at bulk-defined germline het SNPs
# ==========================================================
import argparse
import csv
import gzip
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

try:
    import yaml 
except Exception as e:
    yaml = None

try:
    import matplotlib.pyplot as plt 
except Exception:
    plt = None


BASES = ("A", "C", "G", "T")


@dataclass(frozen=True)
class Params:
    min_cell_dp: int
    monoallelic_af_cut: float
    bin_edges: List[int]


def eprint(*args):
    print(*args, file=sys.stderr)


def load_yaml_params(path: str) -> Params:
    if yaml is None:
        raise RuntimeError("PyYAML not available. Install pyyaml or replace YAML with CLI args.")
    with open(path, "r") as f:
        d = yaml.safe_load(f)
    min_cell_dp = int(d.get("min_cell_dp", 6))
    monoallelic_af_cut = float(d.get("monoallelic_af_cut", 0.95))
    bin_edges = d.get("bin_edges", [0, 5, 10, 15, 20, 30, 50, 1000])
    bin_edges = [int(x) for x in bin_edges]
    if len(bin_edges) < 2 or sorted(bin_edges) != bin_edges:
        raise ValueError("bin_edges must be a sorted list with >=2 entries")
    return Params(min_cell_dp=min_cell_dp, monoallelic_af_cut=monoallelic_af_cut, bin_edges=bin_edges)


def read_samples_tsv(samples_tsv: str) -> List[str]:
    cell_ids = []
    with open(samples_tsv, "r", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"id", "type"}
        if not required.issubset(set(r.fieldnames or [])):
            raise ValueError(f"{samples_tsv} missing required columns: {required}")
        for row in r:
            if row["type"] == "cell":
                cell_ids.append(row["id"])
    if not cell_ids:
        raise ValueError("No cell rows found in samples.tsv")
    return cell_ids


def open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def load_hets_refalt(hets_vcf_gz: str) -> Dict[Tuple[str, int], Tuple[str, str]]:
    """
    Load REF/ALT from VCF as site dictionary keyed by (CHROM, POS).
    Assumes biallelic SNPs; ALT may have multiple alleles -> uses first.
    """
    site = {}
    with open_maybe_gzip(hets_vcf_gz) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3].upper()
            alt = parts[4].split(",")[0].upper()
            if ref in BASES and alt in BASES:
                site[(chrom, pos)] = (ref, alt)
    if not site:
        raise ValueError(f"No SNP sites loaded from {hets_vcf_gz}")
    return site


def depth_bin(dp: int, edges: List[int]) -> str:
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        if lo <= dp < hi:
            return f"[{lo},{hi})"
    return f"[{edges[-2]},{edges[-1]})"


def is_monoallelic(dp: int, af: float, params: Params) -> bool:
    return (af <= 1.0 - params.monoallelic_af_cut) or (af >= params.monoallelic_af_cut)


def read_cell_counts(counts_path, hets, params) -> Tuple[List[Tuple[int, float]], int, int]:
    """
    Return list of (dp, af) for sites present in hets.
    AF computed as alt_count/dp using REF/ALT mapping from VCF.
    """
    out = []
    n_total = 0
    n_used = 0

    if not os.path.exists(counts_path):
        raise FileNotFoundError(f"Missing counts file: {counts_path}")

    with open_maybe_gzip(counts_path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("chrom\tpos"):
                continue

            chrom, pos_s, ref_mp, dp_s, a, c, g, t = line.rstrip("\n").split("\t")
            n_total += 1

            try:
                pos = int(pos_s)
            except ValueError:
                continue

            key = (chrom, pos)
            if key not in hets:
                continue

            ref, alt = hets[key]
            dp = int(dp_s)
            if dp < params.min_cell_dp:
                continue

            base_counts = {"A": int(a), "C": int(c), "G": int(g), "T": int(t)}
            alt_count = base_counts.get(alt, 0)
            n_used += 1
            af = alt_count / dp
            out.append((dp, af))

    return out, n_total, n_used


def write_tsv(path: str, rows: List[dict], fieldnames: List[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def plot_cell(cell_id: str, table_rows: List[dict], out_dir: str) -> None:
    if plt is None:
        eprint("[WARN] matplotlib not available; plots will not be generated.")
        return

    xs, ys = [], []
    for r in table_rows:
        b = r["depth_bin"]
        lo = int(b.split("[", 1)[1].split(",", 1)[0])
        hi = int(b.split(",", 1)[1].split(")", 1)[0])
        mid = (lo + hi) / 2.0
        xs.append(mid)
        ys.append(float(r["monoallelic_rate"]))
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xs, ys, marker="o")
    ax.set_title(f"Monoallelic rate vs depth: {cell_id}")
    ax.set_xlabel("Depth bin midpoint")
    ax.set_ylabel("Monoallelic rate")
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, f"monoallelic_rate_vs_depth_{cell_id}.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, help="config/samples.tsv")
    ap.add_argument("--hets_vcf", required=True, help="work/hets/chr3.hets.tight.vcf.gz")
    ap.add_argument("--counts_dir", required=True, help="work/counts/")
    ap.add_argument("--params", required=True, help="config/params.yaml")
    ap.add_argument("--out_table", required=True, help="results/tables/ado_by_depth.tsv")
    ap.add_argument("--out_dir", required=True, help="results/plots/ (optional plots)")
    ap.add_argument("--no_plots", action="store_true", help="Skip plots even if matplotlib available")
    args = ap.parse_args()

    params = load_yaml_params(args.params)
    cell_ids = read_samples_tsv(args.samples)
    hets = load_hets_refalt(args.hets_vcf)

    all_rows = []
    for cell_id in cell_ids:
        counts_path = os.path.join(args.counts_dir, f"{cell_id}.counts.tsv.gz")
        dp_af, n_total, n_used = read_cell_counts(counts_path, hets, params)

        agg = {}
        for dp, af in dp_af:
            b = depth_bin(dp, params.bin_edges)
            if b not in agg:
                agg[b] = [0, 0]
            agg[b][0] += 1
            if is_monoallelic(dp, af, params):
                agg[b][1] += 1

        bin_labels = []
        for i in range(len(params.bin_edges) - 1):
            lo, hi = params.bin_edges[i], params.bin_edges[i + 1]
            bin_labels.append(f"[{lo},{hi})")

        cell_rows = []
        for b in bin_labels:
            n_sites, n_mono = agg.get(b, [0, 0])
            rate = (n_mono / n_sites) if n_sites > 0 else "NA"
            row = {
                "cell_id": cell_id,
                "depth_bin": b,
                "n_sites": n_sites,
                "n_monoallelic": n_mono,
                "monoallelic_rate": rate,
            }
            all_rows.append(row)
            cell_rows.append(row)

        if (not args.no_plots) and plt is not None:
            plot_rows = [r for r in cell_rows if r["monoallelic_rate"] != "NA"]
            if plot_rows:
                plot_cell(cell_id, plot_rows, args.out_dir)

    fieldnames = ["cell_id", "depth_bin", "n_sites", "n_monoallelic", "monoallelic_rate"]
    write_tsv(args.out_table, all_rows, fieldnames)


if __name__ == "__main__":
    main()