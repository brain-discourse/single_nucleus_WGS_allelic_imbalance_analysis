#!/bin/bash
# ==========================================================
# author: Meethila Gade
# created on: 31-01-2026-15h-45m
# github: https://github.com/brain-discourse
# copyright: 2026
# description: recurrence across cell
# ==========================================================
#SBATCH -J J7_recur
#SBATCH -o logs/J7_recur.%j.out
#SBATCH -e logs/J7_recur.%j.err
#SBATCH -t 02-00:00:00
#SBATCH -c 1
#SBATCH --mem=18G
module load samtools

if [[ -n "${EXTRA_PATH:-}" ]]; then export PATH="${EXTRA_PATH}:${PATH}"; fi

# Required: `case` label 
: "${CASE:?Set CASE to a case label from the manifest, e.g. sbatch --export=ALL,CASE=CASE_A}"
# Project root
PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
DIR="${DIR:-$PROJECT_ROOT}"
MANIFEST="${MANIFEST:-$DIR/config/samples.tsv}"

BULK_SAMPLE=$(awk -F'\t' -v c="$CASE" 'NR>1 && $2=="bulk" && $7==c{print $1; exit}' "$MANIFEST")

VCF=$DIR/work/hets/${BULK_SAMPLE}.chr3.hets.tight.vcf.gz
PARAMS=$DIR/config/params2-plot3.yaml
CONT_TBL=$DIR/results/tables/${CASE}.contiguity_by_cell.tsv

OUT_TBL=$DIR/results/tables/${CASE}.recurrence_windows-hap-3.3.tsv
OUT_REC_DIR=$DIR/work/recurrence/${CASE}

mkdir -p $DIR/results/tables "$OUT_REC_DIR"

CASE_SAMPLES=$DIR/work/qc/${CASE}.cells.tsv
awk -F'\t' -v c="$CASE" 'NR==1{print; next} $2=="cell" && $7==c{print}' "$MANIFEST" > "$CASE_SAMPLES"

python $DIR/scripts/recurrence_screen_hapdr.py \
  --contiguity_table "$CONT_TBL" \
  --params "$PARAMS" \
  --out_table "$OUT_TBL"

# plot 3.3 - changed secment windows from 2-->3 , min_dir_mode_frac: 0.75-->0.8
