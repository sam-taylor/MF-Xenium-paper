"""Cell-type taxonomy of the cohort object: names at the three tiers, broad labels and display labels.

Columns of the object: ct_tier1_final, ct_tier2_final, ct_tier3_final, ct_broad_final (names below),
is_ct_final_resolved (False for Indeterminate and T/NK mixed).
"""
import pandas as pd

INDETERMINATE = "Indeterminate"
TNK_MIXED = "T_NK_Mixed"
UNRESOLVED = [INDETERMINATE, TNK_MIXED]

# tier 3, display label, tier 2, tier 2 label, tier 1, broad label, panel group (heatmap blocks)
_ROWS = [
    ("Megakaryocyte_ITGB3_TGFB1", "MK", "Megakaryocyte_ITGB3_TGFB1", "MK", "Myeloid", "Megakaryocyte", "Megakaryocyte"),
    ("Progenitor_Proliferating_TFRC_MYC", "Prog Prolif", "Progenitor_Proliferating_TFRC_MYC", "Prog Prolif", "Erythroid", "Proliferating progenitor", "Erythroid"),
    ("Erythroblast_Proliferating_MKI67_SLC4A1", "Ery Prolif", "Erythroblast_Proliferating_MKI67_SLC4A1", "Ery Prolif", "Erythroid", "Erythroblast", "Erythroid"),
    ("Erythroblast_Mid_SLC4A1_S100A9", "Ery Mid", "Erythroblast_Mid_SLC4A1_S100A9", "Ery Mid", "Erythroid", "Erythroblast", "Erythroid"),
    ("Erythroblast_Late_SLC4A1_ALAS2", "Ery Late", "Erythroblast_Late_SLC4A1_ALAS2", "Ery Late", "Erythroid", "Erythroblast", "Erythroid"),
    ("Myeloid_Progenitor_Early_MPO_CTSG", "MyPro Early", "Myeloid_Progenitor_Early_MPO_CTSG", "MyPro Early", "Myeloid", "Myeloid progenitor", "Myeloid Progenitor"),
    ("Granulocyte_Progenitor_MPO_CEACAM8", "Gran Prog", "Myeloid_Progenitor_Early_MPO_CTSG", "MyPro Early", "Myeloid", "Myeloid progenitor", "Myeloid Progenitor"),
    ("Monocyte_VCAN_FCN1", "Mono VCAN", "Monocyte_VCAN_FCN1", "Mono VCAN", "Myeloid", "Monocyte", "Mature Myeloid"),
    ("Neutrophil_Mature_ITGAX_MME", "Neut Mature", "Neutrophil_Mature_ITGAX_MME", "Neut Mature", "Myeloid", "Neutrophil", "Mature Myeloid"),
    ("Neutrophil_Immature_CEACAM8_ITGAM", "Neut Immature", "Neutrophil_Immature_CEACAM8_ITGAM", "Neut Immature", "Myeloid", "Neutrophil", "Mature Myeloid"),
    ("Macrophage_CD163_MPEG1", "Mac CD163", "Macrophage_CD163_MPEG1", "Mac CD163", "Myeloid", "Macrophage", "Mature Myeloid"),
    ("Mast_Cell_GATA2_CPA3", "Mast", "Mast_Cell_GATA2_CPA3", "Mast", "Myeloid", "Mast cell", "Mature Myeloid"),
    ("pDC_IRF8_MPEG1", "pDC", "pDC_IRF8_MPEG1", "pDC", "Myeloid", "pDC", "Mature Myeloid"),
    ("NK_Cell_GNLY_CD56bright", "NK bright", "NK_Cell", "NK", "Lymphoid", "T/NK cell", "Lymphoid"),
    ("NK_Cell_GNLY_CD56dim", "NK dim", "NK_Cell", "NK", "Lymphoid", "T/NK cell", "Lymphoid"),
    ("T_Cell_CD8_Naive", "T CD8 Naive", "T_Cell_CD8", "T CD8", "Lymphoid", "T/NK cell", "Lymphoid"),
    ("T_Cell_CD8_Memory", "T CD8 Mem", "T_Cell_CD8", "T CD8", "Lymphoid", "T/NK cell", "Lymphoid"),
    ("T_Cell_CD4_Naive", "T CD4 Naive", "T_Cell_CD4", "T CD4", "Lymphoid", "T/NK cell", "Lymphoid"),
    ("T_Cell_CD4_Memory", "T CD4 Mem", "T_Cell_CD4", "T CD4", "Lymphoid", "T/NK cell", "Lymphoid"),
    ("T_Cell_Reg", "Treg", "T_Cell_CD4", "T CD4", "Lymphoid", "T/NK cell", "Lymphoid"),
    ("T_Cell_ISG", "T ISG", "T_Cell_Other", "T other", "Lymphoid", "T/NK cell", "Lymphoid"),
    ("T_Cell_Prolif", "T Prolif", "T_Cell_Other", "T other", "Lymphoid", "T/NK cell", "Lymphoid"),
    ("B_Cell_IGHM_MZB1", "B Cell", "B_Cell_IGHM_MZB1", "B Cell", "Lymphoid", "B cell", "Lymphoid"),
    ("Plasma_Cell_IGKC_IGHG1", "Plasma", "Plasma_Cell_IGKC_IGHG1", "Plasma", "Lymphoid", "Plasma cell", "Lymphoid"),
    ("Stromal_CXCL12_LUM", "Stromal LUM", "Stromal_CXCL12_LUM", "Stromal LUM", "Stromal", "Stromal cell", "Stromal / Vascular"),
    ("Endothelial_PECAM1_PLVAP", "Endothelial", "Endothelial_PECAM1_PLVAP", "Endothelial", "Endothelial", "Endothelial", "Stromal / Vascular"),
    ("Pericyte_IGFBP7_SPARCL1", "Pericyte", "Pericyte_IGFBP7_SPARCL1", "Pericyte", "Stromal", "Pericyte", "Stromal / Vascular"),
    ("Adipocyte_PLIN1_CD36", "Adipocyte", "Adipocyte_PLIN1_CD36", "Adipocyte", "Stromal", "Adipocyte", "Stromal / Vascular"),
    ("Osteoclast_MMP9_ITGB3", "Osteoclast", "Osteoclast_MMP9_ITGB3", "Osteoclast", "Myeloid", "Osteoclast", "Stromal / Vascular"),
    (TNK_MIXED, "T/NK mixed", TNK_MIXED, "T/NK mixed", "Lymphoid", "T/NK cell", "Unresolved"),
    (INDETERMINATE, "Indeterminate", INDETERMINATE, "Indeterminate", INDETERMINATE, INDETERMINATE, "Unresolved"),
]
TYPES = pd.DataFrame(_ROWS, columns=["tier3", "label3", "tier2", "label2", "tier1", "broad", "panel_group"])
LABEL3 = dict(zip(TYPES.tier3, TYPES.label3))
LABEL2 = dict(zip(TYPES.tier2, TYPES.label2))


def label3(tier3: str) -> str:
    return LABEL3.get(tier3, tier3)


def label2(tier2: str) -> str:
    return LABEL2.get(tier2, tier2)


# lineages that can carry the clone: megakaryocytic, erythroid and myeloid cells (not lymphoid, not stromal)
CLONE_ELIGIBLE_TIER3 = list(TYPES.tier3[TYPES.panel_group.isin(["Megakaryocyte", "Erythroid", "Myeloid Progenitor", "Mature Myeloid"])])
