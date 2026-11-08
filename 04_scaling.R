###############################################################################
# 04_scaling.R -- computing time and memory as a function of n (p = 3,000)
# Full GBLUP (n <= 8,000) versus the YQ, SQC and Nystrom pipelines (Q = 1,000)
# Timing only: random standardized markers, phenotype with 50 causal markers.
# Revision: adds the exact marker-space GBLUP (p x p cross-products), the Min_GRM
# selection rule (testing set = the 500 individuals that are predicted; the time
# includes the n x 500 cross-relationships) and a second experiment in which p
# varies (n = 8,000). Resumable: if a row of scale_results.csv lacks t_mingrm,
# only the Min_GRM pipeline is timed for that n.
# v6: adds the Nystrom pipeline to the experiment in which p varies.
# Resumable: only missing columns are timed.
# v8: adds the fully efficient phenotype-quantile compression (YQ-FE, Q = 1,000:
# Stage 1 group means and within-group SSCP, Stage 2 marker-space REML) to both
# experiments (columns t_yqfe).
###############################################################################
source("00_functions.R")
P <- 3000; Q <- 1000
out <- file.path("..", "outputs", "results", "scale_results.csv")
old <- if (file.exists(out) && "t_exact" %in% names(read.csv(out))) read.csv(out) else NULL
if (!is.null(old)) old <- old[, setdiff(names(old), c("t_gva", "t_yqgv"))]
time_mingrm <- function(Z, y, Znew, Q) {
  t0 <- tnow(); ix <- select_mingrm(Z, Znew, Q); f <- fit_wgblup(Z[ix, ], y[ix], rep(1, Q))
  pr <- predict_wgblup(f, Znew); tnow() - t0
}
if (!is.null(old) && !"t_mingrm" %in% names(old)) old$t_mingrm <- NA
if (!is.null(old) && !"t_yqfe" %in% names(old)) old$t_yqfe <- NA
time_yqfe <- function(Z, y, Znew, Q) {
  t0 <- tnow(); S1 <- yq_fe_stage1(Z, y, make_groups("YQ", y, Q)); f <- fit_yq_fe(S1); pr <- predict_lin(f, Znew)
  tnow() - t0
}
rows <- if (is.null(old)) list() else split(old, seq_len(nrow(old)))
for (n in c(2000, 4000, 6000, 8000, 12000, 20000)) {
  set.seed(n)
  if (!is.null(old) && n %in% old$n) {
    k <- which(old$n == n)
    if (!is.na(old$t_mingrm[k]) && !is.na(old$t_yqfe[k])) next
    Z <- matrix(rnorm(n * P), n, P); y <- rowSums(Z[, 1:50]) * 0.3 + rnorm(n); Znew <- Z[1:500, ]
    if (is.na(old$t_mingrm[k])) rows[[k]]$t_mingrm <- time_mingrm(Z, y, Znew, Q)
    if (is.na(old$t_yqfe[k])) rows[[k]]$t_yqfe <- time_yqfe(Z, y, Znew, Q)
    print(rows[[k]]); write.csv(do.call(rbind, rows), out, row.names = FALSE)
    rm(Z); gc(verbose = FALSE); next
  }
  Z <- matrix(rnorm(n * P), n, P)
  y <- rowSums(Z[, 1:50]) * 0.3 + rnorm(n)
  Znew <- Z[1:500, ]
  t_full <- NA
  if (n <= 6000) {   # full n x n GBLUP timed up to n = 6,000 (about 2 min); larger n extrapolated (Figure 5A)
    t0 <- tnow(); f <- fit_wgblup(Z, y, rep(1, n)); pr <- predict_wgblup(f, Znew)
    t_full <- tnow() - t0; rm(f); gc(verbose = FALSE)
  }
  t0 <- tnow(); lab <- make_groups("YQ", y, Q); S <- summarise_groups(Z, y, lab)
  t_yq_sum <- tnow() - t0
  t0 <- tnow(); f <- fit_wgblup(S$Zbar, S$ybar, S$m); pr <- predict_wgblup(f, Znew)
  t_yq_fit <- tnow() - t0
  t0 <- tnow(); PC <- top_pcs(Z, 30, 1); t_pcs <- tnow() - t0
  t0 <- tnow(); lab <- make_groups("SQC", y, Q, PC); t_grp <- tnow() - t0
  t0 <- tnow(); S <- summarise_groups(Z, y, lab); t_sum <- tnow() - t0
  t0 <- tnow(); f <- fit_wgblup(S$Zbar, S$ybar, S$m); pr <- predict_wgblup(f, Znew)
  t_fit <- tnow() - t0
  t0 <- tnow(); fm <- fit_marker_gblup(Z, y); pr <- predict_marker_gblup(fm, Znew)
  t_exact <- tnow() - t0; rm(fm); gc(verbose = FALSE)
  t_mingrm <- time_mingrm(Z, y, Znew, Q)
  t_yqfe <- time_yqfe(Z, y, Znew, Q)
  t0 <- tnow(); ny <- fit_nystrom(Z, y, Q, 1); pr <- predict_nystrom(ny, Znew)
  t_nys <- tnow() - t0; rm(ny); gc(verbose = FALSE)
  rows[[length(rows) + 1]] <- data.frame(n = n, p = P, Q = Q, t_full = t_full,
    t_yq = t_yq_sum + t_yq_fit, t_pcs = t_pcs, t_groups = t_grp, t_summary = t_sum,
    t_fit = t_fit, t_sqc = t_pcs + t_grp + t_sum + t_fit, t_nys = t_nys,
    t_exact = t_exact, t_mingrm = t_mingrm, t_yqfe = t_yqfe, mem_exact_GB = 8 * P * P * 3 / 1e9,
    mem_full_GB = 8 * n * n * 3 / 1e9, mem_nys_GB = 8 * n * Q * 3 / 1e9,
    mem_yq_GB = 8 * (Q * P + 2 * Q * Q) / 1e9, mem_Z_GB = 8 * n * P / 1e9)
  print(rows[[length(rows)]]); write.csv(do.call(rbind, rows), out, row.names = FALSE)
  rm(Z, S, f, PC); gc(verbose = FALSE)
}

## ---- second experiment: n = 8,000, p varies (exact marker-space GBLUP vs YQ) ----
out <- file.path("..", "outputs", "results", "scale_p_results.csv"); n <- 8000
old <- if (file.exists(out)) read.csv(out) else NULL
if (!is.null(old) && !"t_nys" %in% names(old)) old$t_nys <- NA
if (!is.null(old) && !"t_yqfe" %in% names(old)) old$t_yqfe <- NA
rows <- if (is.null(old)) list() else split(old, seq_len(nrow(old)))
for (P in c(1000, 2000, 4000, 6000, 8000)) {
  if (!is.null(old) && P %in% old$p) {
    k <- which(old$p == P); if (!is.na(old$t_nys[k]) && !is.na(old$t_yqfe[k])) next
    set.seed(P); Z <- matrix(rnorm(n * P), n, P); y <- rowSums(Z[, 1:50]) * 0.3 + rnorm(n); Znew <- Z[1:500, ]
    if (is.na(old$t_nys[k])) {
      t0 <- tnow(); ny <- fit_nystrom(Z, y, Q, 1); pr <- predict_nystrom(ny, Znew); rows[[k]]$t_nys <- tnow() - t0
      rm(ny); gc(verbose = FALSE) }
    if (is.na(old$t_yqfe[k])) rows[[k]]$t_yqfe <- time_yqfe(Z, y, Znew, Q)
    print(rows[[k]]); write.csv(do.call(rbind, rows), out, row.names = FALSE); rm(Z); gc(verbose = FALSE); next
  }
  set.seed(P); Z <- matrix(rnorm(n * P), n, P); y <- rowSums(Z[, 1:50]) * 0.3 + rnorm(n); Znew <- Z[1:500, ]
  t0 <- tnow(); fm <- fit_marker_gblup(Z, y); pr <- predict_marker_gblup(fm, Znew); t_exact <- tnow() - t0
  rm(fm); gc(verbose = FALSE)
  t0 <- tnow(); S <- summarise_groups(Z, y, make_groups("YQ", y, Q)); f <- fit_wgblup(S$Zbar, S$ybar, S$m)
  pr <- predict_wgblup(f, Znew); t_yq <- tnow() - t0
  t0 <- tnow(); ny <- fit_nystrom(Z, y, Q, 1); pr <- predict_nystrom(ny, Znew); t_nys <- tnow() - t0; rm(ny)
  rows[[length(rows) + 1]] <- data.frame(n = n, p = P, Q = Q, t_exact = t_exact, t_yq = t_yq,
    mem_exact_GB = 8 * P * P * 3 / 1e9, mem_full_GB = 8 * n * n * 3 / 1e9, mem_yq_GB = 8 * (Q * P + 2 * Q * Q) / 1e9,
    t_nys = t_nys, t_yqfe = time_yqfe(Z, y, Znew, Q))
  print(rows[[length(rows)]]); write.csv(do.call(rbind, rows), out, row.names = FALSE)
  rm(Z, S, f); gc(verbose = FALSE)
}
