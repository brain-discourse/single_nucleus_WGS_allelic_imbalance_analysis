#!/bin/bash
# ==========================================================
# author: Meethila Gade
# created on: 31-01-2026-15h-32m
# github: https://github.com/brain-discourse
# copyright: 2026
# description: stat summary-continue after phasing
# ==========================================================
#SBATCH -J A3_phase
#SBATCH -o logs/A3_phase.%j.out
#SBATCH -e logs/A3_phase.%j.err
#SBATCH -t 02-00:00:00
#SBATCH -c 2
#SBATCH --mem=80G
module load python/3.12.4
module load samtools
module load gatk


if [[ -n "${EXTRA_PATH:-}" ]]; then export PATH="${EXTRA_PATH}:${PATH}"; fi

# Project root:
PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
DIR="${DIR:-$PROJECT_ROOT}"
MANIFEST="${MANIFEST:-$DIR/config/samples.tsv}"
REF="${REF_FASTA:-}" 

mkdir -p "$DIR/work/phase" "$DIR/logs"

# One job per case
: "${CASE:?Set CASE to a case label from the manifest, e.g. sbatch --export=ALL,CASE=CASE_A}"
ROW=$(awk -F'\t' -v c="$CASE" 'NR>1 && $2=="bulk" && $7==c{print; exit}' "$MANIFEST")
[[ -n "$ROW" ]] || { echo "[ERROR] No bulk row for CASE=$CASE in $MANIFEST" >&2; exit 2; }

BULK_SAMPLE=$(echo "$ROW" | cut -f1) 
BULK_BAM=$(echo "$ROW" | cut -f3) 

IN_VCF=$DIR/work/hets/${BULK_SAMPLE}.chr3.hets.tight.vcf.gz 
OUT_VCF=$DIR/work/phase/${BULK_SAMPLE}.chr3.hets.tight.phased.vcf.gz 
OUT_STATS=$DIR/work/phase/${BULK_SAMPLE}.chr3.hets.tight.phased.stats.tsv 

if [[ ! -s "${BULK_BAM}.bai" ]]; then 
    samtools index "$BULK_BAM" 
fi 

# generate stats 
whatshap stats --tsv "$OUT_STATS" "$OUT_VCF"