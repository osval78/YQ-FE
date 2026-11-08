###############################################################################
# 12_yq_exact.R -- fully efficient (exact) phenotype-quantile compression YQ-FE
# Verifies that the two-stage YQ-FE reproduces the one-stage GBLUP exactly
# (same REML estimates and same predictions) in three algebraically equivalent
# forms, and compares it with the original (approximate) YQ and with the reduced
# versions YQ-FE(k) that keep only the k leading within-group components.
#   Data: first replicate of simulated scenario S1 (n = 6,400, p = 3,000) and the
#   first real-data partition (n = 6,434, p = 5,527); YQ groups, Q = 1,000
#   (exactness also checked for Q = 250 and 3,000 with the marker form).
# Outputs: outputs/results/yq_fe_checks.csv, yq_fe_reduced.csv
###############################################################################
source("00_functions.R")
resd <- file.path("..", "outputs", "results")
D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds"))
pt <- readRDS(file.path(resd, "real_partition_0.rds"))
ds <- make_sim_rep(1, 0)
DATA <- list(Simulation_S1 = list(Ztr = ds$Z[ds$tr, ], ytr = ds$y[ds$tr], Zte = ds$Z[ds$te, ], yte = ds$y[ds$te]),
             Real_wheat = list(Ztr = D$Z[pt$tr, ], ytr = D$y[pt$tr], Zte = D$Z[pt$te, ], yte = D$y[pt$te]))
rm(D, ds); gc(verbose = FALSE)
chk <- list(); red <- list()
for (nm in names(DATA)) {
  X <- DATA[[nm]]; message(nm)
  t0 <- tnow(); ff <- fit_wgblup(X$Ztr, X$ytr, rep(1, length(X$ytr))); pf <- predict_wgblup(ff, X$Zte); tf <- tnow() - t0
  rec <- function(method, Q, f, pr, t, nobs) chk[[length(chk) + 1]] <<- data.frame(data = nm, method = method, Q = Q,
    delta = f$delta, mu = if (is.null(f$mu)) f$mu_reml else f$mu, r_test = cor(pr, X$yte), slope = cov(pr, X$yte) / var(pr),
    max_abs_diff_vs_full = max(abs(pr - pf)), time = t, n_obs_stage2 = nobs)
  rec("One-stage GBLUP (n x n)", length(X$ytr), list(delta = ff$delta, mu = ff$mu_reml), pf, tf, length(X$ytr))
  for (Q in c(250, 1000, 3000)) {
    lab <- make_groups("YQ", X$ytr, Q)
    t0 <- tnow(); S1 <- yq_fe_stage1(X$Ztr, X$ytr, lab); t1 <- tnow() - t0
    t0 <- tnow(); fa <- fit_yq_fe(S1); ta <- tnow() - t0
    rec("YQ-FE, marker form", Q, fa, predict_lin(fa, X$Zte), t1 + ta, NA)
    t0 <- tnow(); S <- summarise_groups(X$Ztr, X$ytr, lab); fy <- fit_wgblup(S$Zbar, S$ybar, S$m)
    rec("YQ (original, approximate)", Q, list(delta = fy$delta, mu = fy$mu_reml), predict_wgblup(fy, X$Zte), tnow() - t0, Q)
    if (Q == 1000) {
      t0 <- tnow(); fb <- fit_yq_fe_kernel(X$Ztr, X$ytr, lab); tb <- tnow() - t0
      rec("YQ-FE, two-stage data form (means + within contrasts)", Q, fb, predict_lin(fb, X$Zte), tb, length(X$ytr))
      t0 <- tnow(); fc <- fit_yq_fe_compact(S1); tc <- tnow() - t0
      rec("YQ-FE, compact form (means + r_W components + residual)", Q, fc, predict_lin(fc, X$Zte), t1 + tc, fc$nobs)
      for (k in c(0, 250, 500, 1000, 2000, fc$rW)) {
        t0 <- tnow(); fk <- fit_yq_fe_compact(S1, k = k, exact = (k == fc$rW)); tk <- tnow() - t0
        pk <- predict_lin(fk, X$Zte)
        red[[length(red) + 1]] <- data.frame(data = nm, Q = Q, k = fk$k, rW = fc$rW, exact = (k == fc$rW),
          n_obs_stage2 = Q + fk$k, delta = fk$delta, r_test = cor(pk, X$yte), pct_full = 100 * cor(pk, X$yte) / cor(pf, X$yte),
          slope = cov(pk, X$yte) / var(pk), mse_ratio = mean((pk - X$yte)^2) / mean((pf - X$yte)^2),
          max_abs_diff_vs_full = max(abs(pk - pf)), time = t1 + tk)
        message(sprintf("   k = %d done", k))
      }
    }
    message(sprintf("  Q = %d done", Q)); gc(verbose = FALSE)
  }
}
chk <- do.call(rbind, chk); red <- do.call(rbind, red)
print(chk, digits = 6); print(red, digits = 4)
write.csv(chk, file.path(resd, "yq_fe_checks.csv"), row.names = FALSE)
write.csv(red, file.path(resd, "yq_fe_reduced.csv"), row.names = FALSE)
