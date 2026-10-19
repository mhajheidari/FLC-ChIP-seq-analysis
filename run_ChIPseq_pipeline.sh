#!/bin/bash
# shellcheck disable=SC1091,SC2154
# (config.sh is created by the user at runtime from config.sh.example and
# is gitignored, so shellcheck cannot follow it (SC1091) or see the
# variables it defines (SC2154) -- RAW_DIRS, GENOME_FASTA, OUTDIR, names,
# ids, pairs, MOCK_SAMPLES, etc. throughout this file all come from there.)
set -e
set -u
set -o pipefail

# ==============================================
# ChIP-seq analysis pipeline
# FastQC -> fastp -> Bowtie2 -> samtools -> MACS3 -> MultiQC
#
# All paths and sample information live in config.sh (not tracked in
# git -- see config.sh.example for the template). This script itself
# should not need editing to run on a different machine or dataset.
# ==============================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config.sh"

if [ ! -f "${CONFIG_FILE}" ]; then
  echo "ERROR: config.sh not found next to this script." >&2
  echo "  Copy the template and edit it first:" >&2
  echo "    cp config.sh.example config.sh" >&2
  echo "    nano config.sh" >&2
  exit 1
fi
# shellcheck source=config.sh
source "${CONFIG_FILE}"

# ==============================================
# CHECK REQUIRED TOOLS ARE ON PATH
# ==============================================
REQUIRED_TOOLS=(fastqc fastp bowtie2 bowtie2-build samtools macs3 multiqc)
missing=0
for tool in "${REQUIRED_TOOLS[@]}"; do
  if ! command -v "${tool}" &> /dev/null; then
    echo "ERROR: required tool '${tool}' not found on PATH." >&2
    missing=1
  fi
done
if [ "${missing}" -ne 0 ]; then
  echo "Install the missing tool(s) (e.g. via conda/mamba) and re-run." >&2
  exit 1
fi

# ==============================================
# CREATE DIRECTORIES
# ==============================================
mkdir -p "${OUTDIR}/fastq_merged"
mkdir -p "${OUTDIR}/fastqc_raw"
mkdir -p "${OUTDIR}/fastp"
mkdir -p "${OUTDIR}/fastqc_trimmed"
mkdir -p "${OUTDIR}/bam"
mkdir -p "${OUTDIR}/logs"
mkdir -p "${OUTDIR}/macs3"
mkdir -p "${OUTDIR}/Bowtie2Index"
mkdir -p "${OUTDIR}/multiqc"

# ==============================================
# 1. MERGE TECHNICAL REPLICATES ACROSS RUN DIRECTORIES
# ==============================================
echo "=== Merging technical replicates ==="
for i in "${!names[@]}"; do
  name=${names[$i]}
  id=${ids[$i]}

  # Find all R1 and R2 files for this sample across every configured run directory
  r1_files=$(find "${RAW_DIRS[@]}" -maxdepth 1 -name "*_${id}_*_R1_001.fastq.gz" | sort)
  r2_files=$(find "${RAW_DIRS[@]}" -maxdepth 1 -name "*_${id}_*_R2_001.fastq.gz" | sort)

  if [ -z "${r1_files}" ] || [ -z "${r2_files}" ]; then
    echo "  WARNING: no FASTQ files found for ${name} (id ${id}) in RAW_DIRS -- skipping" >&2
    continue
  fi

  # shellcheck disable=SC2086
  # (r1_files/r2_files intentionally holds a space-separated list of
  # multiple filenames; word-splitting here is required so cat receives
  # each as a separate argument.)
  cat ${r1_files} > "${OUTDIR}/fastq_merged/${name}_R1.fastq.gz"
  # shellcheck disable=SC2086
  cat ${r2_files} > "${OUTDIR}/fastq_merged/${name}_R2.fastq.gz"

  echo "  ${name}: merged $(echo "${r1_files}" | wc -w) R1 file(s)"
done

# ==============================================
# 2. FASTQC ON RAW (MERGED) READS
# ==============================================
echo "=== Running FastQC on raw merged reads ==="
fastqc \
  --threads "${THREADS}" \
  --outdir "${OUTDIR}/fastqc_raw" \
  "${OUTDIR}"/fastq_merged/*_R1.fastq.gz "${OUTDIR}"/fastq_merged/*_R2.fastq.gz

# ==============================================
# 3. ADAPTER TRIMMING / QUALITY FILTERING WITH FASTP
# ==============================================
echo "=== Trimming and quality-filtering with fastp ==="
for name in "${names[@]}"; do
  echo "  fastp: ${name}"
  fastp \
    --thread "${THREADS}" \
    --in1 "${OUTDIR}/fastq_merged/${name}_R1.fastq.gz" \
    --in2 "${OUTDIR}/fastq_merged/${name}_R2.fastq.gz" \
    --out1 "${OUTDIR}/fastp/${name}_R1.trimmed.fastq.gz" \
    --out2 "${OUTDIR}/fastp/${name}_R2.trimmed.fastq.gz" \
    --json "${OUTDIR}/fastp/${name}.fastp.json" \
    --html "${OUTDIR}/fastp/${name}.fastp.html" \
    2> "${OUTDIR}/logs/${name}_fastp.log"
done

# ==============================================
# 4. FASTQC ON TRIMMED READS (POST-FASTP QC)
# ==============================================
echo "=== Running FastQC on trimmed reads ==="
fastqc \
  --threads "${THREADS}" \
  --outdir "${OUTDIR}/fastqc_trimmed" \
  "${OUTDIR}"/fastp/*_R1.trimmed.fastq.gz "${OUTDIR}"/fastp/*_R2.trimmed.fastq.gz

# ==============================================
# 5. BUILD BOWTIE2 INDEX
# ==============================================
INDEX_PREFIX="${OUTDIR}/Bowtie2Index/TAIR10"
if [ ! -f "${INDEX_PREFIX}.1.bt2" ]; then
  echo "=== Building Bowtie2 index ==="
  bowtie2-build "${GENOME_FASTA}" "${INDEX_PREFIX}"
else
  echo "=== Bowtie2 index already exists, skipping build ==="
fi

# ==============================================
# 6. MAP TRIMMED READS WITH BOWTIE2 AND PROCESS WITH SAMTOOLS
# ==============================================
echo "=== Mapping with Bowtie2 ==="
for name in "${names[@]}"; do
  echo "  Processing ${name}..."
  bowtie2 \
    -p "${THREADS}" \
    --very-sensitive \
    -x "${INDEX_PREFIX}" \
    -1 "${OUTDIR}/fastp/${name}_R1.trimmed.fastq.gz" \
    -2 "${OUTDIR}/fastp/${name}_R2.trimmed.fastq.gz" \
    --no-unal \
    --no-discordant \
    --no-mixed \
    -S "${OUTDIR}/bam/${name}.sam" \
    2> "${OUTDIR}/logs/${name}_bowtie2.log"

  # Convert SAM to sorted BAM and index
  samtools view -bS -@ "${THREADS}" "${OUTDIR}/bam/${name}.sam" | \
    samtools sort -@ "${THREADS}" -o "${OUTDIR}/bam/${name}.sorted.bam"
  samtools index "${OUTDIR}/bam/${name}.sorted.bam"

  # QC: flagstat
  samtools flagstat "${OUTDIR}/bam/${name}.sorted.bam" > "${OUTDIR}/logs/${name}_flagstat.txt"

  # Remove large SAM file
  rm "${OUTDIR}/bam/${name}.sam"
done

# ==============================================
# 7. PEAK CALLING WITH MACS3
# ==============================================
echo "=== Calling peaks with MACS3 ==="
for pair in "${pairs[@]}"; do
  chip=${pair%%:*}
  input=${pair##*:}

  echo "  ${chip} vs ${input}"

  mkdir -p "${OUTDIR}/macs3/${chip}"

  macs3 callpeak \
    -t "${OUTDIR}/bam/${chip}.sorted.bam" \
    -c "${OUTDIR}/bam/${input}.sorted.bam" \
    -f BAMPE \
    -g "${GENOME_SIZE}" \
    -n "${chip}" \
    --outdir "${OUTDIR}/macs3/${chip}" \
    -q 0.05
done

# ==============================================
# 8. POOL INPUT REPLICATES AS BACKGROUND FOR MOCK IP PEAK CALLING
# ==============================================
# C1/C2 (mock IP) have no dedicated matched input of their own, so the
# input replicates from the ChIP:input pairs above are pooled to
# provide a background of comparable depth/composition.
echo "=== Pooling input replicates as background for mock IP peak calling ==="

mkdir -p "${OUTDIR}/bam_pooled"
POOLED_INPUT="${OUTDIR}/bam_pooled/input_pooled.sorted.bam"

if [ ! -f "${POOLED_INPUT}" ]; then
  input_bams=()
  for pair in "${pairs[@]}"; do
    input_bams+=("${OUTDIR}/bam/${pair##*:}.sorted.bam")
  done

  samtools merge -@ "${THREADS}" -f \
    "${OUTDIR}/bam_pooled/input_pooled.bam" \
    "${input_bams[@]}"
  samtools sort -@ "${THREADS}" -o "${POOLED_INPUT}" "${OUTDIR}/bam_pooled/input_pooled.bam"
  samtools index "${POOLED_INPUT}"
  rm "${OUTDIR}/bam_pooled/input_pooled.bam"
  echo "  Pooled input created: ${POOLED_INPUT}"
else
  echo "  Pooled input already exists, skipping merge"
fi

# ==============================================
# 9. PEAK CALLING FOR MOCK IP SAMPLES VS POOLED INPUT
# ==============================================
echo "=== Calling peaks for mock IP samples vs pooled input ==="
for mock in "${MOCK_SAMPLES[@]}"; do
  echo "  ${mock} vs pooled input"

  mkdir -p "${OUTDIR}/macs3/${mock}"

  macs3 callpeak \
    -t "${OUTDIR}/bam/${mock}.sorted.bam" \
    -c "${POOLED_INPUT}" \
    -f BAMPE \
    -g "${GENOME_SIZE}" \
    -n "${mock}" \
    --outdir "${OUTDIR}/macs3/${mock}" \
    -q 0.05
done

# ==============================================
# 10. AGGREGATE QC WITH MULTIQC
# ==============================================
echo "=== Aggregating QC reports with MultiQC ==="
multiqc "${OUTDIR}" \
  --outdir "${OUTDIR}/multiqc" \
  --force

echo "=== Pipeline finished ==="
echo "MultiQC summary report: ${OUTDIR}/multiqc/multiqc_report.html"
echo "Next: run find_and_annotate_reproducible_peaks.R to compute reproducible,"
echo "mock-filtered peak sets and annotate them, then proceed to motif analysis."
