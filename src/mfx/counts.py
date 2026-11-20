"""Integer transcript counts of the cohort object."""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem, sparse_dataset

from . import cohort

# allele-specific HLA probes mark host/donor cells and are therefore confounded with clone identity
HLA_ALLELE_PROBES = {"HLA_DPA1*02:01:01:01_b1+b2", "HLA_DPB1*02:01:02:01_b1+b2", "HLA_DPB1*13:01:01:01_b1+b2"}


def raw_gex(rows=None):
    """Raw integer counts of the gene-expression panel as CSR (cells x genes) and the gene names.

    Mutation probes, negative controls and HLA allele probes are dropped. rows: optional sorted array of
    row numbers to read."""
    g = cohort.open_group()
    raw_var = list(read_elem(g["raw"]["var"]).index)
    gex = set(cohort.var_names())
    keep = [i for i, name in enumerate(raw_var) if name in gex and name not in HLA_ALLELE_PROBES]
    X = sparse_dataset(g["raw"]["X"])
    X = X[:] if rows is None else X[np.asarray(rows)]
    return sp.csr_matrix(X)[:, keep], [raw_var[i] for i in keep]


def layer_counts(rows=None, genes=None):
    """layers['counts'] (the 465 embedding features) as CSR, optionally for rows and a gene subset."""
    g = cohort.open_group()
    X = sparse_dataset(g["layers"]["counts"])
    X = X[:] if rows is None else X[np.asarray(rows)]
    X = sp.csr_matrix(X)
    names = cohort.var_names()
    if genes is not None:
        X = X[:, [names.get_loc(x) for x in genes]]
        names = pd.Index(genes)
    return X, names


def log_norm(rows=None, genes=None):
    """X of the cohort object (counts scaled to 50 per cell, log1p) as CSR for rows and a gene subset."""
    g = cohort.open_group()
    X = sparse_dataset(g["X"])
    X = X[:] if rows is None else X[np.asarray(rows)]
    X = sp.csr_matrix(X)
    names = cohort.var_names()
    if genes is not None:
        X = X[:, [names.get_loc(x) for x in genes]]
        names = pd.Index(genes)
    return X, names
