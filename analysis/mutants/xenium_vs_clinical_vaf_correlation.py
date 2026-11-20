"""
Fig. 2f, ED2h, ED2i: Xenium JAK2 V617F VAF vs. clinical NGS (MSK IMPACT) VAF: correlation, and probe
true-positive / false-discovery (ROC/FDR) analysis.

Xenium VAF = JAK2_ALT_probe_count / (JAK2_ALT_probe_count + JAK2_WT_probe_count),
computed from raw per-sample JAK2_ALT_1848:T / JAK2_WT_1848 transcript-spot counts
(every detected molecule for that probe; matches the ALT/WT counts and VAF column
in Supplementary Table 1).

"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats
from pathlib import Path

# Global aesthetics
# pdf.fonttype=42 embeds TrueType fonts in the PDF (text stays editable in Illustrator)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["font.family"]  = "Arial"
matplotlib.rcParams["font.size"]    = 12

from mfx import paths

OUT_DIR = paths.results("vaf")

# Hard-coded data
# (sample_id, xenium_jak2_vaf, clinical_ngs_vaf, group)
# xenium_vaf = JAK2_ALT_raw / (JAK2_ALT_raw + JAK2_WT_raw)
# clinical_vaf = JAK2 V617F VAF from MSK IMPACT / clinical NGS
# group: MFUT | Pre-Tx | Remission | Relapse | CTMS

DATA_MF8 = [
    # MF8 (JAK2 V617F + TET2 L1830Hfs*4 + TP53 R273C)
    ("MF8-T2",  0.255551,  0.281,  "Pre-Tx"),
    ("MF8-T3",  0.006320,  0.000,  "Relapse"),
    ("MF8-T4",  0.217575,  0.000,  "Relapse"),
    ("MF8-T5",  0.633157,  0.070,  "Relapse"),
]

# Main cohort (MF8 excluded from regression)
DATA = [
    # Transplant cohort
    # MF1  (JAK2 V617F + EZH2 N673K + NF1 L1064P + NUP98 P649R)
    ("MF1-T1",   0.607843,  0.949,  "Pre-Tx"),
    ("MF1-T2",   0.058824,  0.000,  "Remission"),
    ("MF1-T3",   0.050847,  0.000,  "Remission"),
    ("MF1-T4",   0.014113,  0.000,  "Remission"),
    # MF3  (JAK2 V617F + DNMT3A R882H + BCOR + EZH2)
    ("MF3-T1",   0.440353,  0.424,  "Pre-Tx"),
    ("MF3-T2",   0.001380,  0.000,  "Remission"),
    ("MF3-T3",   0.010538,  0.000,  "Remission"),
    ("MF3-T4",   0.015198,  0.130,  "Remission"),
    # MF4  (JAK2 V617F + ASXL1 Y974* + SRSF2 + IDH2 R140W + CBL + PTEN)
    ("MF4-T1",   0.304779,  0.428,  "Pre-Tx"),
    ("MF4-T2",   0.028871,  0.000,  "Relapse"),
    ("MF4-T3",   0.085003,  0.085,  "Relapse"),
    ("MF4-T4",   0.282068,  0.410,  "Relapse"),
    # MF5  (ASXL1/NRAS; non-JAK2-driven)
    ("MF5-T1",   0.086476,  0.139,  "Pre-Tx"),
    ("MF5-T2",   0.021339,  0.000,  "Remission"),
    ("MF5-T3",   0.014788,  0.000,  "Remission"),
    ("MF5-T4",   0.005510,  0.000,  "Remission"),
    ("MF5-T5",   0.018382,  0.000,  "Remission"),
    # MF6  (JAK2 V617F + ASXL1 G646Wfs*12 + ASXL1 R693* + TET2)
    ("MF6-T2",   0.633929,  0.868,  "Pre-Tx"),
    ("MF6-T3",   0.030177,  0.000,  "Relapse"),
    ("MF6-T4",   0.120490,  0.153,  "Relapse"),
    ("MF6-T5",   0.395114,  0.649,  "Relapse"),
    ("MF6-T6",   0.008102,  0.000,  "Remission"),
    # MF7  (JAK2 V617F + ASXL1 + NRAS G12D + TET2 + TET3)
    ("MF7-T2",   0.756991,  0.869,  "Pre-Tx"),
    ("MF7-T5",   0.016616,  0.000,  "Remission"),
    # MF10 (JAK2 V617F + ASXL1 E877* + SF3B1 E592K)
    ("MF10-T1",  0.474149,  0.434,  "Pre-Tx"),
    ("MF10-T2",  0.033364,  0.000,  "Relapse"),
    # MF11 (JAK2 V617F + TET2 E808Kfs*5)
    ("MF11-T1",  0.480410,  0.570,  "Pre-Tx"),
    ("MF11-T2",  0.013944,  0.000,  "Remission"),
    # MF12 (JAK2 V617F + TET2 I1105Rfs*24 + LZTR1 + NOTCH2 + PHF6)
    ("MF12-T1",  0.402847,  0.510,  "Pre-Tx"),
    ("MF12-T2",  0.010394,  0.000,  "Relapse"),
    # MF13 (JAK2 V617F + ASXL1 R693* + SRSF2 P95H + IDH2 R140Q + NF1 + TP53)
    ("MF13-T1",  0.479246,  0.188,  "Pre-Tx"),
    # MF16 (JAK2 V617F)
    ("MF16-T1",  0.680751,  0.948,  "Pre-Tx"),
    ("MF16-T2",  0.012295,  0.0021, "Relapse"),
    # MFUT (untreated MPN controls)
    ("MFUT-1",   0.483090,  0.555,  "MFUT"),   # JAK2 55.5%, SRSF2 P95H
    ("MFUT-2",   0.737085,  0.695,  "MFUT"),   # JAK2 69.5%, IDH2 R140Q 50%
    ("MFUT-3",   0.643145,  0.540,  "MFUT"),   # JAK2 54.0%
    ("MFUT-4",   0.177251,  0.414,  "MFUT"),   # JAK2 41.4%, KRAS L19F, CUX1 L385*
    ("MFUT-5",   0.435166,  0.516,  "MFUT"),   # JAK2 51.6%, TET2 T245Pfs*5
    ("MFUT-6",   0.893058,  0.913,  "MFUT"),   # JAK2 91.3%, ASXL1 E635Rfs*15
    ("MFUT-7",   0.365873,  0.568,  "MFUT"),   # JAK2 56.8%
    ("MFUT-8",   0.395783,  0.588,  "MFUT"),   # JAK2 58.8%, CEBPA, DNMT3A R882H, TP53 M246V
    ("MFUT-9",   0.522471,  0.480,  "MFUT"),   # JAK2 48.0%, NRAS G12R, TET2 Q1034*
    ("MFUT-10",  0.336516,  0.310,  "MFUT"),   # JAK2 31.0% (PB), Chr20q-, BCORL1 del
    # CTMS (JAKi-treated MPN)
    ("CTMS4-T1", 0.538350,  0.918,  "CTMS"),   # JAK2 91.8% (pre-trial), ARID2
    ("CTMS7-T2", 0.584362,  0.848,  "CTMS"),   # JAK2 84.8%
    ("CTMS8-T3", 0.502195,  0.739,  "CTMS"),   # JAK2 73.9% (pre-trial), MCL1, NXF1
    ("CTMS9-T1", 0.034091,  0.049,  "CTMS"),   # JAK2 4.9%, TET2 R1516*
    ("CTMS-11",  0.679216,  0.772,  "CTMS"),   # JAK2 77.2%
]

# Group styling
GROUP_ORDER  = ["MFUT", "Pre-Tx", "Remission", "Relapse", "CTMS"]
GROUP_LABELS = {
    "MFUT":      "Untreated",
    "Pre-Tx":    "Pre-transplant",
    "Remission": "Remission",
    "Relapse":   "Relapse",
    "CTMS":      "JAKi-treated",
}
GROUP_COLORS = {
    "MFUT":      "#d62728",
    "Pre-Tx":    "#ff7f0e",
    "Remission": "#1f77b4",
    "Relapse":   "#e377c2",
    "CTMS":      "#2ca02c",
}
GROUP_MARKERS = {
    "MFUT":      "o",
    "Pre-Tx":    "s",
    "Remission": "^",
    "Relapse":   "D",
    "CTMS":      "P",
}

# Unpack main cohort
sample_ids    = [d[0] for d in DATA]
xenium_vafs   = np.array([d[1] for d in DATA])
clinical_vafs = np.array([d[2] for d in DATA])
groups        = [d[3] for d in DATA]

mf8_xen  = np.array([d[1] for d in DATA_MF8])
mf8_clin = np.array([d[2] for d in DATA_MF8])

# Linear regression (main cohort only, MF8 excluded)
slope, intercept, r_value, p_value, _ = stats.linregress(clinical_vafs, xenium_vafs)
r2 = r_value ** 2

# Same fit with a random intercept per patient (several sections per patient)
import pandas as pd
import statsmodels.formula.api as smf
vaf = pd.DataFrame({"xenium": xenium_vafs, "clinical": clinical_vafs, "patient": [s.split("-")[0] for s in sample_ids]})
mixed = smf.mixedlm("xenium ~ clinical", vaf, groups=vaf.patient).fit(reml=True)

print(f"Main cohort (MF8 excluded): n = {len(DATA)}")
print(f"R² = {r2:.4f}   slope = {slope:.4f}   p = {p_value:.3g}")
print(f"Mixed model (patient random intercept): slope = {mixed.params['clinical']:.4f}   p = {mixed.pvalues['clinical']:.3g}")
pd.DataFrame([{"n": len(DATA), "r2": r2, "slope": slope, "intercept": intercept, "p": p_value,
               "mixed_slope": mixed.params["clinical"], "mixed_p": mixed.pvalues["clinical"]}]
             ).to_csv(OUT_DIR / "fig2f_regression.csv", index=False)

# Fig. 2f: Xenium vs. clinical NGS JAK2 V617F VAF correlation (MF8 excluded)
fig, ax = plt.subplots(figsize=(6, 6))

for sid, xv, cv, grp in zip(sample_ids, xenium_vafs, clinical_vafs, groups):
    ax.scatter(cv, xv, color=GROUP_COLORS[grp], marker=GROUP_MARKERS[grp],
               s=110, edgecolors="none", zorder=3)

x_line = np.linspace(0, 1, 300)
ax.plot(x_line, slope * x_line + intercept, color="#808080", lw=1.5, zorder=2)
ax.plot([0, 1], [0, 1], "--", color="black", lw=1.0, zorder=1)

annot = f"R² = {r2:.3f}\np = {p_value:.0e}\nslope = {slope:.3f}\nn = {len(DATA)}"
ax.text(0.60, 0.04, annot, transform=ax.transAxes, ha="left", va="bottom", fontsize=15)

ax.set_xlim(-0.02, 1.03)
ax.set_ylim(-0.02, 1.03)
ax.set_aspect("equal", adjustable="box")  # equal axis lengths so y=x renders at a true 45 deg
ax.set_xlabel("Clinical NGS JAK2 V617F VAF", fontsize=15)
ax.set_ylabel("Xenium JAK2 V617F VAF\n[ALT / (ALT + WT) transcripts]", fontsize=15)
ax.text(-0.16, 1.05, "f", transform=ax.transAxes, fontsize=22, fontweight="bold")

legend_elements = [
    Line2D([0],[0], marker=GROUP_MARKERS[g], color="w",
           markerfacecolor=GROUP_COLORS[g], markeredgecolor="none", markersize=11, label=GROUP_LABELS[g])
    for g in GROUP_ORDER
] + [
    Line2D([0],[0], color="#808080", lw=1.5, label="Linear fit"),
    Line2D([0],[0], color="black", lw=1.0, ls="--", label="y = x"),
]
leg = ax.legend(handles=legend_elements, fontsize=12, loc="center left", bbox_to_anchor=(1.03, 0.5),
                 frameon=True, edgecolor="black", title="Group", title_fontsize=13)
leg.get_title().set_fontweight("bold")
ax.tick_params(labelsize=13)
ax.spines[["top","right"]].set_visible(False)

plt.tight_layout(pad=0.4)
fig.savefig(OUT_DIR / "xenium_vs_clinical_vaf_correlation.pdf", dpi=300, bbox_inches="tight")
fig.savefig(OUT_DIR / "xenium_vs_clinical_vaf_correlation.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved Fig. 2f to {OUT_DIR}")

# ED2h, ED2i: ROC and FDR analysis of the JAK2 ALT probe (all samples, incl. MF8)
#
# Ground-truth positive (JAK2+): clinical_vaf > 0 (JAK2 V617F confirmed)
# Ground-truth negative (JAK2-): clinical_vaf == 0 (remission / non-JAK2 patient)
all_xen  = np.concatenate([xenium_vafs, mf8_xen])
all_clin = np.concatenate([clinical_vafs, mf8_clin])

is_pos = all_clin > 0
is_neg = all_clin == 0

thresholds = np.linspace(0, 1, 500)
tpr_list, fpr_list, fdr_list = [], [], []
for t in thresholds:
    tp = np.sum(is_pos & (all_xen >= t))
    fp = np.sum(is_neg & (all_xen >= t))
    fn = np.sum(is_pos & (all_xen < t))
    tn = np.sum(is_neg & (all_xen < t))
    tpr_list.append(tp / (tp + fn) if (tp + fn) > 0 else 0)
    fpr_list.append(fp / (fp + tn) if (fp + tn) > 0 else 0)
    fdr_list.append(fp / (tp + fp) if (tp + fp) > 0 else 0)

tpr = np.array(tpr_list)
fpr = np.array(fpr_list)
fdr = np.array(fdr_list)

auc = float(np.trapezoid(tpr, fpr) * -1)  # fpr is in descending order, so negate

j_scores = tpr - fpr
opt_idx  = np.argmax(j_scores)
opt_t    = thresholds[opt_idx]
opt_tpr  = tpr[opt_idx]
opt_fpr  = fpr[opt_idx]
opt_fdr  = fdr[opt_idx]

print(f"\nProbe TP/FP analysis (all samples, n={len(all_xen)})")
print(f"  JAK2+ (clinical_vaf > 0): {is_pos.sum()}")
print(f"  JAK2− (clinical_vaf = 0): {is_neg.sum()}")
print(f"  ROC AUC: {auc:.3f}")
print(f"  Optimal threshold (Youden): {opt_t:.3f}  TPR={opt_tpr:.3f}  FPR={opt_fpr:.3f}  FDR={opt_fdr:.3f}")

fig3, axes = plt.subplots(1, 2, figsize=(11, 4.8))

ax_roc = axes[0]
ax_roc.plot(fpr, tpr, color="#1f77b4", lw=2.5)
ax_roc.scatter([opt_fpr], [opt_tpr], color="#d62728", s=120, zorder=5)
ax_roc.plot([0,1],[0,1],"--", color="#aaaaaa", lw=1.2)
ax_roc.text(0.97, 0.22, f"AUC = {auc:.3f}", transform=ax_roc.transAxes,
            ha="right", va="top", fontsize=14)
ax_roc.text(0.97, 0.13,
            f"Optimal Threshold = {opt_t:.2f}\n(TPR={opt_tpr:.2f}, FPR={opt_fpr:.2f})",
            transform=ax_roc.transAxes, ha="right", va="top", fontsize=12)
ax_roc.set_xlabel("False positive rate (1-specificity)", fontsize=14)
ax_roc.set_ylabel("True positive rate (sensitivity)", fontsize=14)
ax_roc.text(-0.16, 1.05, "d", transform=ax_roc.transAxes, fontsize=22, fontweight="bold")
ax_roc.set_xlim(-0.02, 1.02); ax_roc.set_ylim(-0.02, 1.02)
ax_roc.tick_params(labelsize=12)
ax_roc.spines[["top","right"]].set_visible(False)

ax_fdr = axes[1]
ax_fdr.plot(thresholds, fdr, color="#ff7f0e", lw=2.5)
ax_fdr.axvline(opt_t, color="#d62728", ls="--", lw=1.5)
ax_fdr.axhline(opt_fdr, color="#d62728", ls="--", lw=1.5)
ax_fdr.axhline(0.05, color="#9467bd", ls=":", lw=1.5)
ax_fdr.text(opt_t + 0.02, 0.97, f"Optimal VAF Threshold = {opt_t:.0%}",
            transform=ax_fdr.get_xaxis_transform(), ha="left", va="top",
            fontsize=12, fontweight="bold")

pos_xv = all_xen[is_pos]
neg_xv = all_xen[is_neg]
ax_fdr.scatter(pos_xv, np.full_like(pos_xv, -0.020), marker="|",
               color="#ff7f0e", s=80, lw=2, zorder=4)
ax_fdr.scatter(neg_xv, np.full_like(neg_xv, -0.030), marker="|",
               color="#1f77b4", s=80, lw=2, zorder=4)

ax_fdr.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
ax_fdr.set_xlabel("VAF threshold tested", fontsize=14)
ax_fdr.set_ylabel("False Discovery Rate\n(FP / calls above threshold)", fontsize=14)
ax_fdr.text(-0.16, 1.05, "e", transform=ax_fdr.transAxes, fontsize=22, fontweight="bold")
ax_fdr.set_xlim(-0.02, 1.02); ax_fdr.set_ylim(-0.03, 1.05)

legend_elements_e = [
    Line2D([0],[0], color="#d62728", lw=1.5, ls="--", label=f"FDR = {opt_fdr:.3f}"),
    Line2D([0],[0], color="#9467bd", lw=1.5, ls=":", label="FDR = 0.05"),
    mpatches.Patch(color="#ff7f0e", label=f"JAK2+ (n={is_pos.sum()})"),
    mpatches.Patch(color="#1f77b4", label=f"JAK2− (n={is_neg.sum()})"),
]
ax_fdr.legend(handles=legend_elements_e, fontsize=11, loc="upper right",
              frameon=True, edgecolor="black")
ax_fdr.tick_params(labelsize=12)
ax_fdr.spines[["top","right"]].set_visible(False)

plt.tight_layout(pad=0.4)
fig3.savefig(OUT_DIR / "jak2_probe_tp_fp_analysis.pdf", dpi=300, bbox_inches="tight")
fig3.savefig(OUT_DIR / "jak2_probe_tp_fp_analysis.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved ED2h, ED2i to {OUT_DIR}")
