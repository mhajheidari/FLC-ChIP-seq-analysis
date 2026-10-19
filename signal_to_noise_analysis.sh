#!/bin/bash
# shellcheck disable=SC1091,SC2154
set -euo pipefail

# ============================================================
# Signal-to-noise metrics for FLC ChIP-seq
#   FRiP, deepTools fingerprint/JS distance, strand
#   cross-correlation (optional), and fold enrichment at the
#   reproducible, mock-IP-filtered peak set.
# ============================================================
# Reads its paths and sample sheet from config.sh (see
# config.sh.example), the same configuration file used by
# run_ChIPseq_pipeline.sh and run_motif_enrichment.sh.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config.sh"
if [ ! -f "${CONFIG_FILE}" ]; then
  echo "ERROR: ${CONFIG_FILE} not found. Copy config.sh.example to config.sh and edit it first." >&2
  exit 1
fi
# shellcheck source=/dev/null
source "${CONFIG_FILE}"

# ------------------------------------------------------------
# Tool checks
# ------------------------------------------------------------
REQUIRED_TOOLS=(samtools bedtools bc plotFingerprint awk sort)
for tool in "${REQUIRED_TOOLS[@]}"; do
  if ! command -v "${tool}" >/dev/null 2>&1; then
    echo "ERROR: required tool '${tool}' not found on PATH." >&2
    exit 1
  fi
done

# ------------------------------------------------------------
# Paths and sample sheet (derived from config.sh)
# ------------------------------------------------------------
SNR_DIR="${OUTDIR}/signal_to_noise"
mkdir -p "${SNR_DIR}"

PEAKS="${OUTDIR}/reproducible_peaks/P_reproducible_peaks_mockfiltered.bed"
BAMDIR="${OUTDIR}/bam"

# "pairs" (e.g. A1:Ai1 A2:Ai2 ...) and "MOCK_SAMPLES" come from config.sh.
CHIP_SAMPLES=()
INPUT_SAMPLES=()
for pair in "${pairs[@]}"; do
  CHIP_SAMPLES+=("${pair%%:*}")
  INPUT_SAMPLES+=("${pair##*:}")
done
# MOCK_SAMPLES is already defined in config.sh.

JSD_REF="${BAMDIR}/${INPUT_SAMPLES[0]}.sorted.bam"

# ------------------------------------------------------------
# 0. Sanity checks
# ------------------------------------------------------------
if [ ! -f "${PEAKS}" ]; then
  echo "ERROR: peak file not found: ${PEAKS}" >&2
  echo "       (expected output of find_and_annotate_reproducible_peaks.R, P1-P4 ChIP replicates)" >&2
  exit 1
fi

N_PEAKS=$(wc -l < "${PEAKS}")
echo "Peak file: ${PEAKS}"
echo "Number of reproducible, mock-IP-filtered peaks: ${N_PEAKS}"
echo

# ------------------------------------------------------------
# 1. FRiP
# ------------------------------------------------------------
echo "=== Computing FRiP ==="
FRIP_TABLE="${SNR_DIR}/FRiP_values.tsv"
echo -e "Sample\tGroup\tReads_in_peaks\tTotal_reads\tFRiP" > "${FRIP_TABLE}"

compute_frip () {
  local sample=$1
  local group=$2
  local bam="${BAMDIR}/${sample}.sorted.bam"

  if [ ! -f "${bam}" ]; then
    echo "  [skip] ${bam} not found"
    return
  fi

  local inpeaks
  inpeaks=$(samtools view -c -F 0x904 -f 0x2 -L "${PEAKS}" "${bam}")
  local total
  total=$(samtools view -c -F 0x904 -f 0x2 "${bam}")

  local frip
  frip=$(echo "scale=6; ${inpeaks}/${total}" | bc)

  echo -e "${sample}\t${group}\t${inpeaks}\t${total}\t${frip}" >> "${FRIP_TABLE}"
  echo "  ${sample} (${group}): FRiP = ${frip}"
}

for s in "${CHIP_SAMPLES[@]}";  do compute_frip "$s" "ChIP";  done
for s in "${INPUT_SAMPLES[@]}"; do compute_frip "$s" "Input"; done
for s in "${MOCK_SAMPLES[@]}";  do compute_frip "$s" "MockIP"; done

# ------------------------------------------------------------
# 2. Fingerprint plot (deepTools)
# ------------------------------------------------------------
echo
echo "=== Generating fingerprint plot ==="

FP_ARGS=()
LABELS=()
for s in "${CHIP_SAMPLES[@]}";  do FP_ARGS+=("${BAMDIR}/${s}.sorted.bam"); LABELS+=("ChIP_${s}");  done
for s in "${INPUT_SAMPLES[@]}"; do FP_ARGS+=("${BAMDIR}/${s}.sorted.bam"); LABELS+=("Input_${s}"); done
for s in "${MOCK_SAMPLES[@]}";  do FP_ARGS+=("${BAMDIR}/${s}.sorted.bam"); LABELS+=("MockIP_${s}"); done

plotFingerprint \
  -b "${FP_ARGS[@]}" \
  --labels "${LABELS[@]}" \
  --extendReads \
  --skipZeros \
  --JSDsample "${JSD_REF}" \
  --outQualityMetrics "${SNR_DIR}/fingerprint_metrics.txt" \
  --plotFile "${SNR_DIR}/fingerprint.pdf" \
  --outRawCounts "${SNR_DIR}/fingerprint_raw_counts.txt" \
  --numberOfProcessors "${THREADS:-6}"

plotFingerprint \
  -b "${FP_ARGS[@]}" \
  --labels "${LABELS[@]}" \
  --extendReads --skipZeros \
  --JSDsample "${JSD_REF}" \
  --plotFile "${SNR_DIR}/fingerprint.png" \
  --numberOfProcessors "${THREADS:-6}"

# ------------------------------------------------------------
# 3. Strand cross-correlation (optional, paired-end)
# ------------------------------------------------------------
echo
echo "=== Strand cross-correlation (phantompeakqualtools) ==="

SPP=$(command -v run_spp.R 2>/dev/null || true)
if [ -n "${SPP}" ]; then
  for s in "${CHIP_SAMPLES[@]}"; do
    bam="${BAMDIR}/${s}.sorted.bam"
    Rscript "${SPP}" \
      -c="${bam}" \
      -p="${THREADS:-6}" \
      -savp="${SNR_DIR}/spp_${s}.pdf" \
      -out="${SNR_DIR}/spp_${s}.txt" \
      -odir="${SNR_DIR}"
  done
else
  echo "  [skip] run_spp.R not found in PATH (phantompeakqualtools not installed)."
fi

# ------------------------------------------------------------
# 4. Fold enrichment restricted to the reproducible peak set
# ------------------------------------------------------------
echo
echo "=== Fold enrichment at the ${N_PEAKS} reproducible peaks ==="

FE_TABLE="${SNR_DIR}/fold_enrichment_median.tsv"
echo -e "Sample\tN_peaks_overlapping\tMedian_fold_enrichment" > "${FE_TABLE}"

for s in "${CHIP_SAMPLES[@]}"; do
    np="${OUTDIR}/macs3/${s}/${s}_peaks.narrowPeak"
    if [ ! -f "${np}" ]; then
        echo "  [skip] ${np} not found"
        continue
    fi

    intersected=$(mktemp)
    bedtools intersect -a "${np}" -b "${PEAKS}" -wa -u > "${intersected}"

    n=$(wc -l < "${intersected}")
    if [ "${n}" -eq 0 ]; then
        echo "  ${s}: no overlapping peaks"
        rm "${intersected}"
        continue
    fi

    # Portable median: sort numerically with the external "sort" command
    # (no reliance on gawk-only features such as asort()).
    median_fe=$(cut -f7 "${intersected}" | sort -n | awk '
        { a[NR] = $1; c = NR }
        END {
            if (c % 2 == 1) { print a[(c + 1) / 2] }
            else            { print (a[c / 2] + a[c / 2 + 1]) / 2 }
        }')
    echo -e "${s}\t${n}\t${median_fe}" >> "${FE_TABLE}"
    echo "  ${s}: median fold enrichment = ${median_fe} (${n} peaks overlapping the reproducible set)"

    rm "${intersected}"
done

echo
echo "=== Done. Outputs in ${SNR_DIR} ==="
