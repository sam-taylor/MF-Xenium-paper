"""Niche groups of the cohort object: 20 k-means niches (obs niche_clr50_km_k20_final) joined by Ward linkage into 8
groups (obs niche_group8_clr50_km_k20_final_label)."""

GROUPS_K20 = {
    "Perivascular MK":            [13, 14],
    "Adipo-vascular stroma":      [11],
    "Stromal-granulocytic":       [1, 2, 6],
    "Granulopoietic":             [0, 7, 9, 19],
    "Erythro-myeloid progenitor": [15, 18],
    "Erythropoietic":             [3, 4, 8, 17],
    "Late erythroid-B":           [12, 16],
    "Low-density / unresolved":   [5, 10],
}
GROUP_ORDER = list(GROUPS_K20)
GROUP_OF_K20 = {k: g for g, ks in GROUPS_K20.items() for k in ks}
K20_ORDER = [k for ks in GROUPS_K20.values() for k in ks]
