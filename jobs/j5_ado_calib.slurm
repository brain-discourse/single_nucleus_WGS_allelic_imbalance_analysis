#!/bin/bash
# ==========================================================
# author: Meethila Gade
# created on: 31-01-2026-15h-40m
# github: https://github.com/brain-discourse
# copyright: 2026
# description: covg binned monoallelic rate
# ==========================================================
#SBATCH -J A5_ado
#SBATCH -o logs/A5_ado.%j.out
#SBATCH -e logs/A5_ado.%j.err
#SBATCH -t 02-00:00:00
#SBATCH -c 2
#SBATCH --mem=18G

module load samtools

if [[ -n "${EXTRA_PATH:-}" ]]; then export PATH="${EXTRA_PATH}:${PATH}"; fi

: "${CASE:?Set CASE to a case label from the manifest, e.g. sbatch --export=ALL,CASE=CASE_A}"
# Project root:
PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
DIR="${DIR:-$PROJECT_ROOT}"
MANIFEST="${MANIFEST:-$DIR/config/samples.tsv}"

echo "CASE=$CASE" 

BULK_SAMPLE=$(awk -F'\t' -v c="$CASE" 'NR>1 && $2=="bulk" && $7==c{print $1; exit}' "$MANIFEST") 

VCF=$DIR/work/hets/${BULK_SAMPLE}.chr3.hets.tight.vcf.gz 
PARAMS=$DIR/config/params.yaml 
OUT_TBL=$DIR/results/tables/${CASE}.ado_by_depth.tsv 
OUT_PLOTS_DIR=$DIR/results/plots/${CASE} 

mkdir -p $DIR/results/tables $DIR/work/qc 

# Create case-specific samples list (cells only) 
CASE_SAMPLES=$DIR/work/qc/${CASE}.cells.tsv 
awk -F'\t' -v c="$CASE" 'NR==1{print; next} $2=="cell" && $7==c{print}' "$MANIFEST" > "$CASE_SAMPLES" 

python $DIR/scripts/ado_calibration.py \
  --samples "$CASE_SAMPLES" \
  --hets_vcf "$VCF" \
  --counts_dir "$DIR/work/counts/$CASE" \
  --params "$PARAMS" \
  --out_table "$OUT_TBL" \
  --out_dir "$OUT_PLOTS_DIR" \
  --no_plots