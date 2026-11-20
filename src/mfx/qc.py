"""Per-section outlier flag."""
import numpy as np
import statsmodels.formula.api as smf
from scipy.stats import median_abs_deviation


def _mad_outlier(x, lower, upper):
    med, mad = np.median(x), median_abs_deviation(x)
    return (x < med - lower * mad) | (x > med + upper * mad)


def outlier_flag(obs):
    """True where >= 2 of 4 MAD criteria hold (log counts, log genes, top-20 share, count-gene residual)."""
    fit = smf.ols("log1p_n_genes_by_counts ~ log1p_total_counts", obs).fit()
    resid = obs["log1p_n_genes_by_counts"] - fit.predict()
    n = (_mad_outlier(obs["log1p_total_counts"], 3, 5).astype(int) + _mad_outlier(obs["log1p_n_genes_by_counts"], 3, 5).astype(int)
         + _mad_outlier(obs["pct_counts_in_top_20_genes"], 5, 5).astype(int) + _mad_outlier(resid, 5, 5).astype(int))
    return n >= 2
