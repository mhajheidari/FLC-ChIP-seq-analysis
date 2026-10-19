#!/usr/bin/env Rscript

# ============================================================
# Summary figure: FRiP + fold enrichment + JS distance
# ============================================================
# Reads OUTDIR, CHIP_REPLICATES and MOCK_SAMPLES from config.R
# (see config.R.example), the same configuration file used by
# find_and_annotate_reproducible_peaks.R. Run this script from
# the repository directory that contains config.R, or set the
# CONFIG_FILE environment variable to its path.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

config_file <- Sys.getenv("CONFIG_FILE", unset = "config.R")
if (!file.exists(config_file)) {
  stop(sprintf(
    "config.R not found at '%s'. Copy config.R.example to config.R (next to find_and_annotate_reproducible_peaks.R), edit it, and run this script from that directory, or set CONFIG_FILE=/path/to/config.R.",
    config_file
  ))
}
source(config_file)

# config.R defines CHIP_REPLICATES (e.g. c("A1","A2","A3","A4")) and
# MOCK_SAMPLES (e.g. c("C1","C2")). The matched input controls follow the
# same "<ChIP name>" -> "<letter>i<number>" convention used throughout the
# pipeline (A1 -> Ai1, etc.).
INPUT_SAMPLES <- sub("^([A-Za-z]+)", "\\1i", CHIP_REPLICATES)

SNR_DIR <- file.path(OUTDIR, "signal_to_noise")
if (!dir.exists(SNR_DIR)) {
  stop(sprintf("signal_to_noise directory not found: %s (run signal_to_noise_analysis.sh first)", SNR_DIR))
}
setwd(SNR_DIR)

# Reproducible, mock-IP-filtered peak set (for the panel subtitle; the
# peak count is read from the actual file rather than hardcoded).
PEAKS_FILE <- file.path(OUTDIR, "reproducible_peaks", "P_reproducible_peaks_mockfiltered.bed")
n_peaks <- if (file.exists(PEAKS_FILE)) length(readLines(PEAKS_FILE)) else NA

# ------------------------------------------------------------
# 1. FRiP
# ------------------------------------------------------------
frip <- read.delim("FRiP_values.tsv", stringsAsFactors = FALSE)
frip$Group <- factor(frip$Group, levels = c("ChIP", "Input", "MockIP"))
frip$FRiP_pct <- frip$FRiP * 100

# Preserve sample order: ChIP, Input, MockIP (from config.R)
frip$Sample <- factor(frip$Sample,
                      levels = c(CHIP_REPLICATES, INPUT_SAMPLES, MOCK_SAMPLES))

frip_subtitle <- if (is.na(n_peaks)) {
  "Peak set: reproducible, mock-IP-filtered"
} else {
  sprintf("Peak set: reproducible, mock-IP-filtered (n = %d)", n_peaks)
}

p_frip <- ggplot(frip, aes(x = Sample, y = FRiP_pct, fill = Group)) +
  geom_col(colour = "black", width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  annotate("text", x = 0.6, y = 1.15,
           label = "ENCODE TF threshold (1%)",
           hjust = 0, colour = "red", size = 3) +
  scale_fill_manual(values = c(ChIP = "#2166AC",
                               Input = "#B2182B",
                               MockIP = "#4DAF4A")) +
  labs(x = NULL,
       y = "FRiP (%)",
       title = "A. Fraction of reads in peaks (FRiP)",
       subtitle = frip_subtitle) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

# ------------------------------------------------------------
# 2. Fold enrichment (reproducible peak set)
# ------------------------------------------------------------
fe <- read.delim("fold_enrichment_median.tsv", stringsAsFactors = FALSE)
fe$Sample <- factor(fe$Sample, levels = CHIP_REPLICATES)

fe_title <- if (is.na(n_peaks)) {
  "B. Fold enrichment at the reproducible peaks"
} else {
  sprintf("B. Fold enrichment at the %d reproducible peaks", n_peaks)
}

p_fe <- ggplot(fe, aes(x = Sample, y = Median_fold_enrichment)) +
  geom_col(fill = "#2166AC", colour = "black", width = 0.7) +
  geom_text(aes(label = N_peaks_overlapping),
            vjust = -0.5, size = 3, colour = "grey30") +
  labs(x = NULL,
       y = "Median fold enrichment (MACS3)",
       title = fe_title,
       subtitle = "Number above each bar = peaks overlapping the reproducible set") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ------------------------------------------------------------
# 3. JS distance from fingerprint metrics (optional)
# ------------------------------------------------------------
# deepTools writes a table with a header line starting with "#"
# Columns: Sample, Reads, Fingerprint, JS Distance (of fingerprint), ...
# We only need Sample and JS distance.

fp <- tryCatch({
  raw <- readLines("fingerprint_metrics.txt")
  # Find the header line
  hdr_idx <- grep("^Sample", raw)
  if (length(hdr_idx) == 0) hdr_idx <- grep("^#Sample", raw)
  if (length(hdr_idx) > 0) {
    # Read from header, skipping the leading "#"
    hdr <- sub("^#", "", raw[hdr_idx])
    lines <- c(hdr, raw[(hdr_idx + 1):length(raw)])
    read.delim(text = paste(lines, collapse = "\n"),
               stringsAsFactors = FALSE)
  } else {
    NULL
  }
}, error = function(e) NULL)

if (!is.null(fp) && "Sample" %in% colnames(fp)) {
  # Column 4 is usually JS Distance; find by name if possible
  js_col <- grep("JS", colnames(fp), value = TRUE)[1]
  if (!is.na(js_col)) {
    fp$Group <- ifelse(grepl("^ChIP", fp$Sample), "ChIP",
                ifelse(grepl("^Input", fp$Sample), "Input",
                ifelse(grepl("^MockIP", fp$Sample), "MockIP", "Other")))
    fp$Group <- factor(fp$Group, levels = c("ChIP", "Input", "MockIP"))

    p_js <- ggplot(fp, aes(x = Sample, y = .data[[js_col]], fill = Group)) +
      geom_col(colour = "black", width = 0.7) +
      scale_fill_manual(values = c(ChIP = "#2166AC",
                                   Input = "#B2182B",
                                   MockIP = "#4DAF4A")) +
      labs(x = NULL,
           y = "JS distance vs input",
           title = "C. Fingerprint JS distance",
           subtitle = "Higher = more deviation from input distribution") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  } else {
    p_js <- NULL
  }
} else {
  p_js <- NULL
}

# ------------------------------------------------------------
# 4. Combine panels
# ------------------------------------------------------------
if (!is.null(p_js)) {
  combined <- p_frip / p_fe / p_js +
    plot_annotation(tag_levels = NULL,
                    title = "Signal-to-noise metrics for FLC ChIP-seq")
  ggsave("signal_to_noise_summary.pdf", combined, width = 9, height = 11)
  ggsave("signal_to_noise_summary.png", combined, width = 9, height = 11, dpi = 300)
} else {
  combined <- p_frip / p_fe +
    plot_annotation(title = "Signal-to-noise metrics for FLC ChIP-seq")
  ggsave("signal_to_noise_summary.pdf", combined, width = 9, height = 8)
  ggsave("signal_to_noise_summary.png", combined, width = 9, height = 8, dpi = 300)
}

# Save tables as CSV for supplementary
write.csv(frip, "FRiP_values.csv", row.names = FALSE)
write.csv(fe,   "fold_enrichment_median.csv", row.names = FALSE)

cat("Summary figure written: signal_to_noise_summary.pdf / .png\n")
cat("Tables: FRiP_values.csv, fold_enrichment_median.csv\n")
