# Single-nucleus WGS allelic imbalance analysis

Code for evaluating whether apparent monoallelic signals in whole-genome-amplified single-nucleus DNA reflect stochastic allelic sampling or extend across contiguous, consistently oriented haplotype regions. The workflow was developed for an exploratory analysis of nuclei with high _PIK3CA_ p.E545K variant-allele signals in focal cortical malformations, as described in [Chapter 4 of the dissertation](#citation).

The analysis uses donor-specific germline heterozygous SNPs from matched bulk whole-genome sequencing (WGS) to characterize depth-dependent allelic dropout, single-nucleus allelic balance, haplotype-directional contiguity, and candidate recurrent imbalance to assess evidence of regional, haplotype-consistent imbalance compatible with a loss-of-heterozygosity (LOH) mechanism in the high-variant signal group.

## 1. Scientific question and scope

Whole-genome amplification can generate uneven locus coverage and stochastic loss of representation of one allele. Consequently, a single-nucleus variant-allele fraction near 0 or 1 cannot be assumed to represent a biological copy-number or LOH event. The analysis asks whether nuclei showing unusually strong _PIK3CA_ variant-allele signals also demonstrate coherent allelic imbalance over neighboring phased heterozygous SNPs, beyond isolated sites affected by depth-dependent sampling.

The supplied implementation:

- Restricts SNP-panel construction and downstream analysis to chromosome 3(`chr3`; also embedded in file naming conventions).
- Generates donor-specific heterozygous SNP panels from a pre-existing joint-genotyped bulk WGS VCF and phases them using a matched bulk BAM.
- Counts single-nucleus reference and alternate alleles at panel SNPs, summarizes depth-stratified monoallelic observations, and evaluates consistent directional haplotype signals across sliding windows.
- Collapses qualifying windows into per-nucleus candidate segments. A separate script supports window-level recurrence screening.

## 2. Analysis methods

| Stage                                    | Implementation                                                               | Main operation                                                                                                                                                                                                                                                                                              |
| ---------------------------------------- | ---------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Heterozygous SNP selection**           | `jobs/j1_build_het_panel.sh`                                                 | `bcftools` selects chromosome 3 biallelic SNPs for one bulk sample with `PASS`/`.` filter status, heterozygous genotype, DP ≥ 20, GQ ≥ 30, and ≥ 5 reads for each allele. A stricter panel additionally requires bulk alternate-allele fraction of 0.30–0.70. Produces VCF, BED, and SNP-density summaries. |
| **Read-backed phasing**                  | `jobs/j2_phase_het_panel_whatshap.sh`; `jobs/j3_phase_het_panel_whatshap.sh` | WhatsHap phases the strict panel with the matched bulk BAM and a reference FASTA. `j2` compresses the phased VCF and calculates phase statistics; `j3` recalculates statistics from an existing phased VCF.                                                                                                 |
| **Single-nucleus allele counting**       | `jobs/j4_cell_allele_counts.slurm`; `scripts/allele_counts.py`               | Uses `samtools mpileup -aa -q 20 -Q 20` at heterozygous panel sites and derives per-base counts (`chrom`, `pos`, `ref`, `dp`, `A`, `C`, `G`, `T`).                                                                                                                                                          |
| **Allelic dropout calibration**          | `jobs/j5_ado_calib.slurm`; `scripts/ado_calibration.py`                      | At sites passing `min_cell_dp`, calculates alternate-allele fraction and depth-binned monoallelic rates using `monoallelic_af_cut`. This is strictly used as a calibration.                                                                                                                                 |
| **Haplotype contiguity**                 | `jobs/j6_contiguity_tracks.slurm`; `scripts/contiguity_runs.py`              | Assigns directional allelic support at eligible sites, then summarizes sliding windows (`window_nsnp`, `step_nsnp`) by longest directional run, modal haplotype direction, directional counts/fractions, and monoallelic fractions.                                                                         |
| **Candidate segment construction**       | `jobs/j7_recurrence_screen.slurm`; `scripts/recurrence_screen_hapdr.py`      | Applies usable-site, directional-run, directional-fraction, and directional-count criteria to windows and merges adjacent passing windows with matching direction into per-nucleus segments.                                                                                                                |
| **Optional cross-nucleus window screen** | `scripts/recurrence_screen.py`                                               | Counts nuclei meeting directional-window criteria with a shared unique modal direction and applies `min_cells_recur`.                                                                                                                                                                                       |

### Interpreting the key thresholds

| Setting                        | Level                    | Interpretation                                                                                                               |
| ------------------------------ | ------------------------ | ---------------------------------------------------------------------------------------------------------------------------- |
| `monoallelic_af_cut = 0.90`    | SNP site                 | A site is classified as monoallelic when alternate-allele fraction is ≥ 0.90 or ≤ 0.10.                                      |
| `min_major_allele_frac = 0.70` | SNP site                 | In `params2-plot3.yaml`, this is the minimum major-allele fraction for assigning a directional label in contiguity analysis. |
| `min_dir_mode_frac = 0.80`     | Window/segment screening | The modal haplotype direction must account for ≥ 80% of directional sites in a qualifying window.                            |

## 4. Software

The workflow requires the following command-line tools:

- `samtools`, `bcftools`, `bgzip`, and `tabix` (htslib)
- `whatshap` for read-backed phasing
- Python 3 and PyYAML for configuration; matplotlib for optional plots

## 5. Required inputs

| Input                                | Expected format                                                                          | Stages     |
| ------------------------------------ | ---------------------------------------------------------------------------------------- | ---------- |
| Joint-genotyped chromosome 3 SNP VCF | Indexed, bgzipped VCF with bulk sample genotype and `GT`, `AD`, `DP`, `GQ` FORMAT values | `j1`       |
| Matched bulk WGS alignment           | Coordinate-sorted BAM and `.bai`                                                         | `j2`       |
| Single-nucleus WGS alignments        | Coordinate-sorted BAMs and indexes; `j4` may create missing indexes                      | `j4`       |
| Reference genome                     | GRCh38 FASTA and `.fai`, matching the alignment/reference build                          | `j2`, `j4` |
| Sample manifest                      | `config/samples.tsv`, tab-separated                                                      | `j2`–`j7`  |

Sample manifest columns include:

| Column | Name   | Description                                                                         |
| ------ | ------ | ----------------------------------------------------------------------------------- |
| 1      | `id`   | Bulk sample name in the VCF for `bulk` rows; per-nucleus identifier for `cell` rows |
| 2      | `type` | `bulk` or `cell`                                                                    |
| 3      | `bam`  | BAM path                                                                            |
| 4      | `bai`  | BAM index path                                                                      |
| 5      | `vcf`  | path to vcf file                                                                    |
| 6      | `tbi`  | path to sample index                                                                |
| 7      | `case` | Donor/case identifier relating cell rows to their bulk sample                       |

## 6. Configuration and execution

From the repository root:

```bash
cp config/samples.example.tsv config/samples.tsv
export PROJECT_ROOT="$(pwd)"
export REF_FASTA=/path/to/reference/GRCh38.fasta
```

Edit `config/samples.tsv` locally, then set the relevant environment variables:

- `PROJECT_ROOT` or `DIR`: Repository root
- `MANIFEST`: Defaults to `$DIR/config/samples.tsv`
- `REF_FASTA`: Indexed GRCh38 reference FASTA
- `JOINT_VCF_CHR3`, `BULK_SAMPLE`: Indexed bulk joint VCF and the bulk sample to select
- `CASE`: Case label

### Workflow

```text
j1  Build donor-specific chr3 heterozygous SNP panel
---> j2  Phase strict SNP panel with matched bulk WGS BAM (j3: optional statistics)
 ----> j4  Count reference/alternate alleles in single nuclei
  -------> j5  ADO by depth
  -------> j6  Haplotype-directional windows
   ------------------------------------------------------> j7 Qualifying segments per nucleus
```

## 7. Data availability and privacy

This code release excludes individual-level sequencing reads, BAMs, variant calls, real sample manifests and identifiers, sample-specific outputs, and SLURM logs.

## 8. Citation

Gade, M. (2026). _Somatic Mosaicism in Focal Cortical Malformations: From Tissue-Level Variant Discovery to Single-Cell Genotype-Transcriptome Resolution_ [Doctoral dissertation, University of North Carolina at Chapel Hill]. Carolina Digital Repository. https://doi.org/10.17615/39vr-n649

Relevant section: Chapter 4, “Single-nucleus WGS allelic-balance assessment.”
