#!/bin/bash
# ==========================================================
# author: Meethila Gade
# created on: 31-01-2026-15h-37m
# github: https://github.com/brain-discourse
# copyright: 2026
# description: per cell allele counts at het snps
# ==========================================================
#SBATCH -J A4_counts
#SBATCH -o logs/A4_counts.%j.out
#SBATCH -e logs/A4_counts.%j.err
#SBATCH -t 02-00:00:00
#SBATCH -c 2
#SBATCH --mem=8G

module load samtools

if [[ -n "${EXTRA_PATH:-}" ]]; then export PATH="${EXTRA_PATH}:${PATH}"; fi

# Required: CASE label 
: "${CASE:?Set CASE to a case label from the manifest, e.g. sbatch --export=ALL,CASE=CASE_A}"
# Project root: 
PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
DIR="${DIR:-$PROJECT_ROOT}"
MANIFEST="${MANIFEST:-$DIR/config/samples.tsv}"

echo "CASE=$CASE" 

REF="${REF_FASTA:?Set REF_FASTA to the GRCh38 FASTA used for alignment (with .fai)}"

BULK_SAMPLE=$(awk -F'\t' -v c="$CASE" 'NR>1 && $2=="bulk" && $7==c{print $1; exit}' "$MANIFEST") 

BED=$DIR/work/hets/${BULK_SAMPLE}.chr3.hets.tight.bed 
VCF=$DIR/work/hets/${BULK_SAMPLE}.chr3.hets.tight.vcf.gz 

OUTDIR=$DIR/work/counts/$CASE 
mkdir -p "$OUTDIR" 

# One job per case:
while IFS=$'\t' read -r CELL_ID _ CELL_BAM _ <&3; do
  echo "CELL=$CELL_ID"

  [[ -s "${CELL_BAM}.bai" ]] || samtools index "$CELL_BAM" 

  OUT=$OUTDIR/${CELL_ID}.counts.tsv.gz 

  samtools mpileup -f "$REF" -l "$BED" -aa -q 20 -Q 20 "$CELL_BAM" \
    | python3 "$DIR/scripts/allele_counts.py" \
    | bgzip -c > "$OUT" 

  tabix -f -s 1 -b 2 -e 2 "$OUT" 2>/dev/null || true
done 3< <(awk -F'\t' -v c="$CASE" 'NR>1 && $2=="cell" && $7==c' "$MANIFEST")
