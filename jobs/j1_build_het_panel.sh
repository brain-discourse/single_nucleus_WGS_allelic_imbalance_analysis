#!/bin/bash
# ==========================================================
# author: Meethila Gade
# created on: 31-01-2026-15h-25m
# github: https://github.com/brain-discourse
# copyright: 2026
# description: build clean het SNP panel
# ==========================================================
#SBATCH -J A1_hetpanel
#SBATCH -o logs/A1_hetpanel.%j.out
#SBATCH -e logs/A1_hetpanel.%j.err
#SBATCH -t 02-00:00:00
#SBATCH -c 2
#SBATCH --mem=8G
module load samtools
module load gatk

# Project root:
PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
DIR="${DIR:-$PROJECT_ROOT}"
JOINT_VCF_CHR3="${JOINT_VCF_CHR3:?Set JOINT_VCF_CHR3 to the joint-genotyped chr3 SNP VCF (bgzipped, indexed)}"
BULK_SAMPLE="${BULK_SAMPLE:?Set BULK_SAMPLE to the bulk sample name as it appears in the joint VCF}"

mkdir -p $DIR/work/hets $DIR/results/tables

OUTVCF=$DIR/work/hets/${BULK_SAMPLE}.chr3.hets.vcf.gz
OUTBED=$DIR/work/hets/${BULK_SAMPLE}.chr3.hets.bed

OUTVCF_strict=$DIR/work/hets/${BULK_SAMPLE}.chr3.hets.tight.vcf.gz
OUTBED_strict=$DIR/work/hets/${BULK_SAMPLE}.chr3.hets.tight.bed

OUTSUM=$DIR/results/tables/${BULK_SAMPLE}.chr3.hets_density.txt
OUTSUM_strict=$DIR/results/tables/${BULK_SAMPLE}.chr3.hets.tight_density.txt


# broad filter to high-confidence biallelic SNP hets
bcftools view -r chr3 -m2 -M2 -v snps "$JOINT_VCF_CHR3" \
| bcftools view -s "$BULK_SAMPLE" \
| bcftools filter -i '(FILTER="." || FILTER="PASS") && GT="het" && FORMAT/DP>=20 && FORMAT/GQ>=30 && FORMAT/AD[0:0]>=5 && FORMAT/AD[0:1]>=5' \
| bgzip -c > "$OUTVCF"

tabix -p vcf "$OUTVCF"

# strict filter to high-confidence biallelic SNP hets
bcftools view "$OUTVCF" \
| bcftools filter -i '(FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]))>=0.30 && (FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]))<=0.70' \
| bgzip -c > "$OUTVCF_strict"

tabix -p vcf "$OUTVCF_strict"

# make bed of snp sites
bcftools query -f'%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' "$OUTVCF" > "$OUTBED"
bcftools query -f'%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' "$OUTVCF_strict" > "$OUTBED_strict"

# density summary
awk 'BEGIN{OFS="\t"}{c++; if(prev>0){d=$3-prev; sum+=d} prev=$3}
     END{print "nsnp",c; print "mean_bp_spacing",(c>1?sum/(c-1):"NA")}' "$OUTBED" > "$OUTSUM"

awk 'BEGIN{OFS="\t"}{c++; if(prev>0){d=$3-prev; sum+=d} prev=$3}
     END{print "nsnp",c; print "mean_bp_spacing",(c>1?sum/(c-1):"NA")}' "$OUTBED_strict" > "$OUTSUM_strict"

echo "DONE" > $DIR/work/hets/${BULK_SAMPLE}.DONE.ok