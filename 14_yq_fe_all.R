###############################################################################
# 14_yq_fe_all.R -- fully efficient phenotype-quantile compression (YQ-FE,
# Section 2.10) on every simulated replicate (15) and real-data partition (5),
# with the same data, partitions and phenotype-quantile groups as YQ in
# 02/03 (Q = 250, 1000, 2000, 3000). Stage 1 = yq_fe_stage1() (group means and
# pooled within-group SSCP), Stage 2 = fit_yq_fe() (marker form, lme4 REML).
# The exact marker-space GBLUP (identical to the one-stage GBLUP) is refitted in
# each replicate to measure the largest absolute difference between the YQ-FE
# predictions and the one-stage predictions.
# Single core (OPENBLAS_NUM_THREADS=1). Resumable.
# Usage: Rscript 14_yq_fe_all.R [sim|real]
# Outputs: outputs/results/yq_fe_sim_results.csv, yq_fe_real_results.csv
###############################################################################
source("00_functions.R")
resd <- file.path("..", "outputs", "results")
part <- commandArgs(TRUE)[1]; if (is.na(part)) part <- "sim"
QS <- c(250, 1000, 2000, 3000)
eval_fe <- function(Ztr, ytr, Zte, yte, gte) {
  fm <- fit_marker_gblup(Ztr, ytr); pm <- predict_marker_gblup(fm, Zte)
  rows <- list()
  for (Q in QS) {
    t0 <- tnow(); lab <- make_groups("YQ", ytr, Q); S1 <- yq_fe_stage1(Ztr, ytr, lab); t1 <- tnow() - t0
    t0 <- tnow(); f <- fit_yq_fe(S1); pr <- predict_lin(f, Zte); t2 <- tnow() - t0
    rows[[length(rows) + 1]] <- data.frame(method = "YQFE", Q = Q, shrink = "REML", t(metrics(pr, yte, gte)),
      time = t1 + t2, time_stage1 = t1, time_stage2 = t2, groups = Q, delta = f$delta, delta_onestage = fm$delta,
      max_abs_diff_pred_vs_onestage = max(abs(pr - pm)), r_onestage = cor(pm, yte))
    rm(S1, f); gc(verbose = FALSE)
  }
  do.call(rbind, rows)
}
if (part == "sim") {
  out <- file.path(resd, "yq_fe_sim_results.csv")
  rows <- if (file.exists(out)) read.csv(out) else NULL
  done <- if (is.null(rows)) character(0) else unique(paste(rows$scenario, rows$rep))
  for (si in seq_along(SIM$SCEN)) for (rep in seq_len(SIM$NREP) - 1) {
    sc <- SIM$SCEN[si]; if (paste(sc, rep) %in% done) next
    message(sprintf("%s rep %d", sc, rep + 1)); t0 <- tnow()
    d <- make_sim_rep(si, rep)
    res <- eval_fe(d$Z[d$tr, ], d$y[d$tr], d$Z[d$te, ], d$y[d$te], d$g[d$te])
    res$scenario <- sc; res$rep <- rep; res$n_train <- length(d$tr); res$p <- SIM$P
    rows <- rbind(rows, res); write.csv(rows, out, row.names = FALSE)
    rm(d); gc(verbose = FALSE); message(sprintf("  %.0f s", tnow() - t0))
  }
}
if (part == "real") {
  D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds"))
  out <- file.path(resd, "yq_fe_real_results.csv")
  rows <- if (file.exists(out)) read.csv(out) else NULL
  done <- if (is.null(rows)) integer(0) else unique(rows$rep)
  for (rep in 0:4) {
    if (rep %in% done) next
    pt <- readRDS(file.path(resd, sprintf("real_partition_%d.rds", rep)))
    message(sprintf("real partition %d", rep + 1)); t0 <- tnow()
    res <- eval_fe(D$Z[pt$tr, ], D$y[pt$tr], D$Z[pt$te, ], D$y[pt$te], NULL)
    res$scenario <- "real_wheat"; res$rep <- rep; res$n_train <- length(pt$tr); res$p <- ncol(D$Z)
    rows <- rbind(rows, res); write.csv(rows, out, row.names = FALSE)
    message(sprintf("  %.0f s", tnow() - t0)); gc(verbose = FALSE)
  }
}
