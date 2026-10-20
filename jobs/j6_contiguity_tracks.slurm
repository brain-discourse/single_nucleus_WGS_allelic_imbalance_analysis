#!/bin/bash
# ==========================================================
# author: Meethila Gade
# created on: 31-01-2026-15h-43m
# github: https://github.com/brain-discourse
# copyright: 2026
# description: contiguity tracks per cell
# ==========================================================
#SBATCH -J A6_runs
#SBATCH -o logs/A6_runs.%j.out
#SBATCH -e logs/A6_runs.%j.err
#SBATCH -t 02-00:00:00
#SBATCH -c 2
#SBATCH --mem=18G

module load samtools

if [[ -n "${EXTRA_PATH:-}" ]]; then export PATH="${EXTRA_PATH}:${PATH}"; fi

# Required: `case` label 
: "${CASE:?Set CASE to a case label from the manifest, e.g. sbatch --export=ALL,CASE=CASE_A}"
# Project root: 
PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
DIR="${DIR:-$PROJECT_ROOT}"
MANIFEST="${MANIFEST:-$DIR/config/samples.tsv}"

echo "CASE=$CASE"

BULK_SAMPLE=$(awk -F'\t' -v c="$CASE" 'NR>1 && $2=="bulk" && $7==c{print $1; exit}' "$MANIFEST") 

VCF=$DIR/work/hets/${BULK_SAMPLE}.chr3.hets.tight.vcf.gz 
P_VCF=$DIR/work/phase/${BULK_SAMPLE}.chr3.hets.tight.phased.vcf.gz 

PARAMS=$DIR/config/params2-plot3.yaml 
OUT_TBL=$DIR/results/tables/${CASE}.contiguity_by_cell.tsv 
OUT_PLOTS_DIR=$DIR/results/plots/${CASE} 
OUT_WIN_DIR=$DIR/work/windows/${CASE} 

mkdir -p $DIR/results/tables "$OUT_WIN_DIR" 

CASE_SAMPLES=$DIR/work/qc/${CASE}.cells.tsv 
awk -F'\t' -v c="$CASE" 'NR==1{print; next} $2=="cell" && $7==c{print}' "$MANIFEST" > "$CASE_SAMPLES" 

python $DIR/scripts/contiguity_runs.py \
  --samples "$CASE_SAMPLES" \
  --hets_vcf "$VCF" \
  --phased_vcf "$P_VCF" \
  --phased_sample "$BULK_SAMPLE" \
  --counts_dir "$DIR/work/counts/$CASE" \
  --params "$PARAMS" \
  --out_table "$OUT_TBL" \
  --out_tracks_dir "$OUT_WIN_DIR" \
  --out_plots "$OUT_PLOTS_DIR" \
  --no_plots