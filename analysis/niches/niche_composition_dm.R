#!/usr/bin/env Rscript
# Dirichlet-multinomial model (brms) of niche-group composition per section, counts ~ disease group + (1 | patient),
# contrasts of the expected proportions (Fig. 4d): posterior mean (percentage points), 95% interval, P(opposite sign),
# and the spread of each niche group's expected proportion across the four disease groups (posterior SD, percentage points).
# Input:  <work>/results/niche_composition/{niche_group_counts,niche_group_columns}.csv
# Output: same folder: dm_fit.rds, dm_contrasts.csv, dm_stage_variation.csv, dm_diagnostics.csv
suppressPackageStartupMessages({
  library(brms)
  library(data.table)
  library(posterior)
})
source(file.path(Sys.getenv("MFX_ROOT", unset = "."), "src/mfx.R"))

out <- mfx_results("niche_composition")
options(mc.cores = 4)
GROUPS <- c("Pre-Tx", "MFUT", "Remission", "Relapse")

data <- as.data.frame(fread(file.path(out, "niche_group_counts.csv")))
cols <- fread(file.path(out, "niche_group_columns.csv"))
resp <- cols$response
data$counts <- as.matrix(data[, resp])
data$patient_id <- factor(data$patient_id)
data$group <- factor(data$group, levels = GROUPS)

pri <- prior(gamma(2, 0.1), class = phi)
for (dpar in paste0("mu", resp[-1])) {
  pri <- c(pri,
           prior_string("normal(0, 1.5)", class = "b", dpar = dpar),
           prior_string("student_t(3, 0, 1)", class = "sd", dpar = dpar))
}
fit <- brm(
  bf(counts | trials(total) ~ group + (1 | patient_id)),
  data = data, family = dirichlet_multinomial(), prior = pri,
  chains = 4, iter = 4000, warmup = 1000, seed = 20261109,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  file = file.path(out, "dm_fit"), file_refit = "on_change", refresh = 200
)

nd <- data.frame(group = factor(GROUPS, levels = GROUPS), total = 1e6L, patient_id = NA)
ep <- posterior_epred(fit, newdata = nd, re_formula = NA) / 1e6  # draws x group x niche
dimnames(ep)[[2]] <- GROUPS
prop <- function(g) ep[, g, ]
contr <- list("Pre-Tx vs MFUT" = prop("Pre-Tx") - prop("MFUT"),
              "Relapse vs Pre-Tx" = prop("Relapse") - prop("Pre-Tx"),
              "Remission vs Pre-Tx" = prop("Remission") - prop("Pre-Tx"))
res <- rbindlist(lapply(names(contr), function(cn) {
  x <- contr[[cn]]
  data.table(niche = cols$niche, group_contrast = cn,
             diff_pts = 100 * colMeans(x),
             ci_lo_pts = 100 * apply(x, 2, quantile, 0.025),
             ci_hi_pts = 100 * apply(x, 2, quantile, 0.975),
             p_opposite_sign = pmin(colMeans(x < 0), colMeans(x > 0)))
}))
res[, credible_95 := (ci_lo_pts > 0) | (ci_hi_pts < 0)]
base <- data.table(niche = cols$niche, group_contrast = "Pre-Tx level", diff_pts = 100 * colMeans(prop("Pre-Tx")))
fwrite(rbind(res, base, fill = TRUE), file.path(out, "dm_contrasts.csv"))

sd_pts <- 100 * apply(ep, c(1, 3), sd)
fwrite(data.table(niche = cols$niche, sd_pts = apply(sd_pts, 2, median),
                  sd_pts_lo = apply(sd_pts, 2, quantile, 0.025), sd_pts_hi = apply(sd_pts, 2, quantile, 0.975)),
       file.path(out, "dm_stage_variation.csv"))

d <- summarise_draws(as_draws_df(fit))
np <- nuts_params(fit)
fwrite(data.table(max_rhat = max(d$rhat, na.rm = TRUE), min_ess_bulk = min(d$ess_bulk, na.rm = TRUE),
                  divergences = sum(np$Parameter == "divergent__" & np$Value == 1)),
       file.path(out, "dm_diagnostics.csv"))
print(res[order(group_contrast, diff_pts)], digits = 3)
