"""
Standard AAL90/AAL116 atlas labels in canonical Tzourio-Mazoyer 2002 order.

The AAL atlas indexes regions 1..116. Indices 1-90 are cortical (45 paired
L/R regions), 91-108 are cerebellar, 109-116 are vermal.

Source: Tzourio-Mazoyer N et al. (2002). Automated anatomical labeling of
activations in SPM using a macroscopic anatomical parcellation of the MNI
MRI single-subject brain. NeuroImage 15(1):273-289.

This is the standard ordering distributed with all AAL-parcellated DTI
connectomes including the BNU/HCP collection.
"""

# Index → (region_name, hemisphere) where hemisphere is L, R, or M (midline)
AAL116 = [
    (1,   "Precentral",          "L"),
    (2,   "Precentral",          "R"),
    (3,   "Frontal_Sup",         "L"),
    (4,   "Frontal_Sup",         "R"),
    (5,   "Frontal_Sup_Orb",     "L"),
    (6,   "Frontal_Sup_Orb",     "R"),
    (7,   "Frontal_Mid",         "L"),
    (8,   "Frontal_Mid",         "R"),
    (9,   "Frontal_Mid_Orb",     "L"),
    (10,  "Frontal_Mid_Orb",     "R"),
    (11,  "Frontal_Inf_Oper",    "L"),
    (12,  "Frontal_Inf_Oper",    "R"),
    (13,  "Frontal_Inf_Tri",     "L"),
    (14,  "Frontal_Inf_Tri",     "R"),
    (15,  "Frontal_Inf_Orb",     "L"),
    (16,  "Frontal_Inf_Orb",     "R"),
    (17,  "Rolandic_Oper",       "L"),
    (18,  "Rolandic_Oper",       "R"),
    (19,  "Supp_Motor_Area",     "L"),
    (20,  "Supp_Motor_Area",     "R"),
    (21,  "Olfactory",           "L"),
    (22,  "Olfactory",           "R"),
    (23,  "Frontal_Sup_Medial",  "L"),
    (24,  "Frontal_Sup_Medial",  "R"),
    (25,  "Frontal_Med_Orb",     "L"),
    (26,  "Frontal_Med_Orb",     "R"),
    (27,  "Rectus",              "L"),
    (28,  "Rectus",              "R"),
    (29,  "Insula",              "L"),
    (30,  "Insula",              "R"),
    (31,  "Cingulum_Ant",        "L"),
    (32,  "Cingulum_Ant",        "R"),
    (33,  "Cingulum_Mid",        "L"),
    (34,  "Cingulum_Mid",        "R"),
    (35,  "Cingulum_Post",       "L"),
    (36,  "Cingulum_Post",       "R"),
    (37,  "Hippocampus",         "L"),
    (38,  "Hippocampus",         "R"),
    (39,  "ParaHippocampal",     "L"),
    (40,  "ParaHippocampal",     "R"),
    (41,  "Amygdala",            "L"),
    (42,  "Amygdala",            "R"),
    (43,  "Calcarine",           "L"),
    (44,  "Calcarine",           "R"),
    (45,  "Cuneus",              "L"),
    (46,  "Cuneus",              "R"),
    (47,  "Lingual",             "L"),
    (48,  "Lingual",             "R"),
    (49,  "Occipital_Sup",       "L"),
    (50,  "Occipital_Sup",       "R"),
    (51,  "Occipital_Mid",       "L"),
    (52,  "Occipital_Mid",       "R"),
    (53,  "Occipital_Inf",       "L"),
    (54,  "Occipital_Inf",       "R"),
    (55,  "Fusiform",            "L"),
    (56,  "Fusiform",            "R"),
    (57,  "Postcentral",         "L"),
    (58,  "Postcentral",         "R"),
    (59,  "Parietal_Sup",        "L"),
    (60,  "Parietal_Sup",        "R"),
    (61,  "Parietal_Inf",        "L"),
    (62,  "Parietal_Inf",        "R"),
    (63,  "SupraMarginal",       "L"),
    (64,  "SupraMarginal",       "R"),
    (65,  "Angular",             "L"),
    (66,  "Angular",             "R"),
    (67,  "Precuneus",           "L"),
    (68,  "Precuneus",           "R"),
    (69,  "Paracentral_Lobule",  "L"),
    (70,  "Paracentral_Lobule",  "R"),
    (71,  "Caudate",             "L"),
    (72,  "Caudate",             "R"),
    (73,  "Putamen",             "L"),
    (74,  "Putamen",             "R"),
    (75,  "Pallidum",            "L"),
    (76,  "Pallidum",            "R"),
    (77,  "Thalamus",            "L"),
    (78,  "Thalamus",            "R"),
    (79,  "Heschl",              "L"),
    (80,  "Heschl",              "R"),
    (81,  "Temporal_Sup",        "L"),
    (82,  "Temporal_Sup",        "R"),
    (83,  "Temporal_Pole_Sup",   "L"),
    (84,  "Temporal_Pole_Sup",   "R"),
    (85,  "Temporal_Mid",        "L"),
    (86,  "Temporal_Mid",        "R"),
    (87,  "Temporal_Pole_Mid",   "L"),
    (88,  "Temporal_Pole_Mid",   "R"),
    (89,  "Temporal_Inf",        "L"),
    (90,  "Temporal_Inf",        "R"),
    (91,  "Cerebelum_Crus1",     "L"),
    (92,  "Cerebelum_Crus1",     "R"),
    (93,  "Cerebelum_Crus2",     "L"),
    (94,  "Cerebelum_Crus2",     "R"),
    (95,  "Cerebelum_3",         "L"),
    (96,  "Cerebelum_3",         "R"),
    (97,  "Cerebelum_4_5",       "L"),
    (98,  "Cerebelum_4_5",       "R"),
    (99,  "Cerebelum_6",         "L"),
    (100, "Cerebelum_6",         "R"),
    (101, "Cerebelum_7b",        "L"),
    (102, "Cerebelum_7b",        "R"),
    (103, "Cerebelum_8",         "L"),
    (104, "Cerebelum_8",         "R"),
    (105, "Cerebelum_9",         "L"),
    (106, "Cerebelum_9",         "R"),
    (107, "Cerebelum_10",        "L"),
    (108, "Cerebelum_10",        "R"),
    (109, "Vermis_1_2",          "M"),
    (110, "Vermis_3",            "M"),
    (111, "Vermis_4_5",          "M"),
    (112, "Vermis_6",            "M"),
    (113, "Vermis_7",            "M"),
    (114, "Vermis_8",            "M"),
    (115, "Vermis_9",            "M"),
    (116, "Vermis_10",           "M"),
]

# Lobe mapping for coarsening
REGION_TO_LOBE = {
    "Precentral":          "Motor",
    "Frontal_Sup":         "Frontal",
    "Frontal_Sup_Orb":     "Frontal",
    "Frontal_Mid":         "Frontal",
    "Frontal_Mid_Orb":     "Frontal",
    "Frontal_Inf_Oper":    "Frontal",
    "Frontal_Inf_Tri":     "Frontal",
    "Frontal_Inf_Orb":     "Frontal",
    "Rolandic_Oper":       "Motor",
    "Supp_Motor_Area":     "Motor",
    "Olfactory":           "Frontal",
    "Frontal_Sup_Medial":  "Frontal",
    "Frontal_Med_Orb":     "Frontal",
    "Rectus":              "Frontal",
    "Insula":              "Insula",
    "Cingulum_Ant":        "Cingulate",
    "Cingulum_Mid":        "Cingulate",
    "Cingulum_Post":       "Cingulate",
    "Hippocampus":         "Hippocampus",
    "ParaHippocampal":     "Temporal",
    "Amygdala":            "Subcortical",
    "Calcarine":           "Occipital",
    "Cuneus":              "Occipital",
    "Lingual":             "Occipital",
    "Occipital_Sup":       "Occipital",
    "Occipital_Mid":       "Occipital",
    "Occipital_Inf":       "Occipital",
    "Fusiform":            "Temporal",
    "Postcentral":         "Somatosensory",
    "Parietal_Sup":        "Parietal",
    "Parietal_Inf":        "Parietal",
    "SupraMarginal":       "Parietal",
    "Angular":             "Parietal",
    "Precuneus":           "Parietal",
    "Paracentral_Lobule":  "Motor",
    "Caudate":             "Subcortical",
    "Putamen":             "Subcortical",
    "Pallidum":            "Subcortical",
    "Thalamus":            "Thalamus",
    "Heschl":              "Auditory",
    "Temporal_Sup":        "Temporal",
    "Temporal_Pole_Sup":   "Temporal",
    "Temporal_Mid":        "Temporal",
    "Temporal_Pole_Mid":   "Temporal",
    "Temporal_Inf":        "Temporal",
}

# All cerebellar/vermal entries → Cerebellum
for idx, name, hemi in AAL116:
    if name not in REGION_TO_LOBE:
        if "Cerebelum" in name or "Vermis" in name:
            REGION_TO_LOBE[name] = "Cerebellum"
        else:
            REGION_TO_LOBE[name] = "Other"


def index_to_lobe_hemi(idx_1based):
    """Returns (lobe, hemi) for a 1-based AAL index, or None."""
    if idx_1based < 1 or idx_1based > 116:
        return None
    _, name, hemi = AAL116[idx_1based - 1]
    lobe = REGION_TO_LOBE.get(name, "Other")
    return (lobe, hemi)


def coarsen_aal_graph(G_full):
    """Take a graph whose nodes are 1-based AAL indices, return a coarse
    graph whose nodes are 'L-Lobe', 'R-Lobe', 'M-Cerebellum' etc."""
    import networkx as nx
    from collections import defaultdict
    coarse_map = {}
    for n in G_full.nodes():
        lh = index_to_lobe_hemi(n)
        if lh is None: continue
        lobe, hemi = lh
        coarse_map[n] = f"{hemi}-{lobe}"
    Gc = nx.Graph()
    weights = defaultdict(int)
    for u, v in G_full.edges():
        cu = coarse_map.get(u); cv = coarse_map.get(v)
        if cu is None or cv is None: continue
        if cu != cv:
            weights[tuple(sorted([cu, cv]))] += 1
    for (a, b), w in weights.items():
        if w >= 2:    # threshold for AAL — sparse data
            Gc.add_edge(a, b)
    # Drop isolated nodes
    Gc.remove_nodes_from([n for n in list(Gc.nodes()) if Gc.degree(n) == 0])
    return Gc


if __name__ == "__main__":
    print(f"AAL116 atlas with {len(AAL116)} regions defined.")
    print(f"Region-to-lobe coverage: {len(REGION_TO_LOBE)} entries.")
    # Test offset issue: the BNU/HCP files use 0-based indexing.
    # We'll handle both in coarsen.
    lobes = set()
    for idx, name, hemi in AAL116:
        lobe = REGION_TO_LOBE.get(name, "Other")
        lobes.add(f"{hemi}-{lobe}")
    print(f"Distinct coarse nodes: {len(lobes)}")
    for l in sorted(lobes):
        print(f"  {l}")
