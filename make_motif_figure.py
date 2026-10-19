# python make_motif_figure.py \
#  --streme streme_out/streme.html \
#  --tomtom tomtom_out/tomtom.html \
#  --centrimo centrimo_out/centrimo.html \
#  --out figure1
#!/usr/bin/env python3
"""
make_motif_figure.py

Generates Figure 1 (Panel A: sequence-logo comparison; Panel B: CentriMo
central-enrichment curves) directly from the JSON data embedded in the
MEME-Suite HTML output files (streme.html, tomtom.html, centrimo.html).

No image is fabricated: every number and coordinate plotted here is parsed
straight out of the "var data = {...}" block that STREME/TOMTOM/CentriMo
write into their own HTML reports. Re-running this script on your own
machine, against your own output files, reproduces the figure exactly.

Requirements:
    pip install matplotlib numpy pandas logomaker pillow

Usage:
    python make_motif_figure.py \
        --streme streme.html \
        --tomtom tomtom.html \
        --centrimo centrimo.html \
        --out figure1

    (Adjust --query-id / --target-id / --denovo-db / --denovo-id /
     --published-db / --published-id below if your motif IDs differ.)
"""

import argparse
import json
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import logomaker


# ----------------------------------------------------------------------
# Step 1: extract the embedded "var data = {...}" JSON block from a
# MEME-Suite HTML report (works for streme.html, tomtom.html, centrimo.html,
# meme.html, etc. -- they all embed their results the same way).
# ----------------------------------------------------------------------
def load_meme_suite_json(html_path):
    with open(html_path, encoding="utf-8", errors="replace") as f:
        content = f.read()

    start = content.find("var data")
    if start == -1:
        raise ValueError(f"Could not find 'var data' block in {html_path}")
    start = content.find("{", start)

    # Walk the braces to find the matching closing brace (the JSON block
    # can be large and contain nested objects/arrays).
    depth = 0
    end = None
    for i in range(start, len(content)):
        c = content[i]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    if end is None:
        raise ValueError(f"Unbalanced braces while parsing {html_path}")

    return json.loads(content[start:end])


# ----------------------------------------------------------------------
# Step 2: helpers
# ----------------------------------------------------------------------
def revcomp_pwm(pwm, alphabet):
    """Reverse-complement a position-probability matrix (list of [A,C,G,T]
    rows) given the alphabet's symbol order."""
    comp_map = {"A": "T", "C": "G", "G": "C", "T": "A"}
    comp_idx = [alphabet.index(comp_map[a]) for a in alphabet]
    pwm = np.asarray(pwm)
    return pwm[::-1][:, comp_idx]


def pwm_to_info_df(pwm, alphabet):
    df = pd.DataFrame(np.asarray(pwm), columns=alphabet)
    return logomaker.transform_matrix(df, from_type="probability", to_type="information")


def centrimo_curve(motif_entry, seqlen, n_sequences):
    """Return (x, y): x = position relative to peak center (bp),
    y = fraction of sequences with a motif site starting at that position."""
    sites = np.asarray(motif_entry["sites"], dtype=float)
    motif_len = motif_entry["len"]
    n_positions = seqlen - motif_len + 1
    starts = np.arange(n_positions)
    x = starts + motif_len / 2.0 - seqlen / 2.0
    y = sites / n_sequences
    return x, y


def smooth(y, window=7):
    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode="same")


def find_centrimo_motif(centrimo_data, db, motif_id):
    for m in centrimo_data["motifs"]:
        if m["db"] == db and m["id"] == motif_id:
            return m
    raise ValueError(f"CentriMo motif not found: db={db}, id={motif_id}")


# ----------------------------------------------------------------------
# Step 3: main
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--streme", required=True, help="Path to streme.html")
    ap.add_argument("--tomtom", required=True, help="Path to tomtom.html")
    ap.add_argument("--centrimo", required=True, help="Path to centrimo.html")
    ap.add_argument("--out", default="figure1", help="Output file prefix (no extension)")

    # IDs identifying which motifs to plot -- edit these if your top hit
    # differs. Defaults match a de novo STREME motif "1-..." matched via
    # TOMTOM to a published motif at index 4 of the target list, and their
    # CentriMo entries at db=0 (STREME output) and db=3 (published-motif db).
    ap.add_argument("--tomtom-query-index", type=int, default=0,
                     help="Index into tomtom['queries'] for the de novo motif")
    ap.add_argument("--tomtom-target-index", type=int, default=4,
                     help="Index into tomtom['targets'] for the published motif")
    ap.add_argument("--centrimo-denovo-db", type=int, default=0)
    ap.add_argument("--centrimo-denovo-id", default=None,
                     help="Motif id string in centrimo.html; default = auto from TOMTOM query id")
    ap.add_argument("--centrimo-published-db", type=int, default=3)
    ap.add_argument("--centrimo-published-id", default=None,
                     help="Motif id string in centrimo.html; default = auto from TOMTOM target id")

    args = ap.parse_args()

    streme_data = load_meme_suite_json(args.streme)
    tomtom_data = load_meme_suite_json(args.tomtom)
    centrimo_data = load_meme_suite_json(args.centrimo)

    alphabet = [s["symbol"] for s in tomtom_data["alphabet"]["symbols"][:4]]

    # --- Panel A inputs: de novo motif PWM and published motif PWM ---
    query = tomtom_data["queries"][args.tomtom_query_index]
    target = tomtom_data["targets"][args.tomtom_target_index]

    match = None
    for mm in tomtom_data["all_matches"][args.tomtom_query_index]["matches"]:
        if mm["idx"] == args.tomtom_target_index:
            match = mm
            break
    if match is None:
        sys.exit("Could not find a TOMTOM match between the requested query and target indices.")

    pwm_denovo = np.asarray(query["pwm"])
    pwm_published = np.asarray(target["pwm"])
    if match.get("rc"):
        pwm_published_plot = revcomp_pwm(pwm_published, alphabet)
        orientation_note = "reverse complement"
    else:
        pwm_published_plot = pwm_published
        orientation_note = "forward strand"

    info_denovo = pwm_to_info_df(pwm_denovo, alphabet)
    info_published = pwm_to_info_df(pwm_published_plot, alphabet)

    # --- Panel B inputs: CentriMo site-position histograms ---
    denovo_id = args.centrimo_denovo_id or query["id"]
    published_id = args.centrimo_published_id or target["id"]

    m_denovo = find_centrimo_motif(centrimo_data, args.centrimo_denovo_db, denovo_id)
    m_published = find_centrimo_motif(centrimo_data, args.centrimo_published_db, published_id)

    seqlen = centrimo_data["seqlen"]
    n_sequences = len(centrimo_data["sequences"])

    x_denovo, y_denovo = centrimo_curve(m_denovo, seqlen, n_sequences)
    x_published, y_published = centrimo_curve(m_published, seqlen, n_sequences)

    y_denovo_s = smooth(y_denovo)
    y_published_s = smooth(y_published)

    # --- Build the figure ---
    fig = plt.figure(figsize=(7.5, 7.2))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1.1, 1.1, 1.6], hspace=0.65)

    axA1 = fig.add_subplot(gs[0])
    logomaker.Logo(info_denovo, ax=axA1, color_scheme="classic")
    axA1.set_title("De novo motif (STREME)", fontsize=10, loc="left")
    axA1.set_ylabel("bits", fontsize=8)
    axA1.set_xticks([])
    axA1.text(0.01, 0.90, "A", transform=axA1.transAxes, fontsize=13,
               fontweight="bold", va="top", ha="left")

    axA2 = fig.add_subplot(gs[1])
    logomaker.Logo(info_published, ax=axA2, color_scheme="classic")
    axA2.set_title(f"Published motif ({target['id']}/{target.get('alt', '')}; {orientation_note})",
                    fontsize=9, loc="left")
    axA2.set_ylabel("bits", fontsize=8)
    axA2.set_xlabel("Position (bp)", fontsize=8)
    axA2.text(0.99, 1.35, f"TOMTOM E = {match['ev']}", transform=axA2.transAxes,
               fontsize=9, ha="right", va="top", style="italic")

    axB = fig.add_subplot(gs[2])
    axB.plot(x_denovo, y_denovo_s, color="#B22222", lw=1.6, label="De novo motif")
    axB.plot(x_published, y_published_s, color="#1f5fa8", lw=1.6, label="Published motif")
    axB.axvline(0, color="grey", lw=0.8, linestyle="--")
    axB.set_xlabel("Distance from peak center (bp)", fontsize=9)
    axB.set_ylabel("Probability of site", fontsize=9)
    axB.set_title("CentriMo central enrichment", fontsize=10, loc="left")
    axB.legend(fontsize=8, frameon=False, loc="upper right")
    axB.text(0.01, 0.95, "B", transform=axB.transAxes, fontsize=13,
              fontweight="bold", va="top")
    axB.set_xlim(-seqlen / 2, seqlen / 2)

    fig.savefig(f"{args.out}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{args.out}.pdf", bbox_inches="tight")
    print(f"Saved {args.out}.png and {args.out}.pdf")
    print(f"De novo motif: {query['id']}  (n={len(x_denovo)} positions, seqlen={seqlen})")
    print(f"Published motif: {target['id']}  match E-value={match['ev']}, rc={match.get('rc')}")


if __name__ == "__main__":
    main()
