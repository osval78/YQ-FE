###############################################################################
# 17_yq_fen_all.R -- Nystrom-accelerated fully efficient phenotype-quantile
# method (YQ-FEN) on every simulated replicate (15) and real-data partition (5),
# with the same data, partitions and phenotype-quantile groups as YQ and YQ-FE
# (Q = 250, 1000, 2000, 3000). Stage 1 = yq_fe_stage1() (identical to YQ-FE);
# Stage 2 = fit_yq_fen(): Nystrom approximation of the YQ-FE cross-product
# matrix (k = 500 landmark markers by default), moment estimate of delta after
# projecting out the h = 200 leading Nystrom directions, and Nystrom-
# preconditioned conjugate gradients on the YQ-FE mixed-model equations. At Q = 1000 the script also (i) re-times YQ-FE on the
# same machine, (ii) runs YQ-FEN with k = 0 (no Nystrom: plain moments, plain CG), 250 and 1000,
# and (iii) runs the CG solver with the REML delta of YQ-FE (check: it must
# reproduce YQ-FE). The exact one-stage GBLUP (marker space) is the reference.
# Single core (OPENBLAS_NUM_THREADS=1). Resumable.
# Usage: Rscript 17_yq_fen_all.R [sim|real]
# Outputs: outputs/results/yq_fen_sim_results.csv, yq_fen_real_results.csv
###############################################################################
source("00_functions.R")
resd <- file.path("..", "outputs", "results")
part <- commandArgs(TRUE)[1]; if (is.na(part)) part <- "sim"
QS <- c(250, 1000, 2000, 3000)
eval_fen <- function(Ztr, ytr, Zte, yte, gte) {
  t0 <- tnow(); fm <- fit_marker_gblup(Ztr, ytr); pm <- predict_marker_gblup(fm, Zte); t_ex <- tnow() - t0
  rows <- list()
  rec <- function(variant, Q, k, pr, t1, t2, f, extra = list())
    rows[[length(rows) + 1]] <<- data.frame(method = "YQFEN", variant = variant, Q = Q, k = k, shrink = "MOM",
      t(metrics(pr, yte, gte)), time = t1 + t2, time_stage1 = t1, time_stage2 = t2,
      time_moments = f$time_moments, time_nystrom = f$time_nystrom, time_cg = f$time_cg, iter = f$iter,
      delta = f$delta, delta_onestage = fm$delta, max_abs_diff_pred_vs_onestage = max(abs(pr - pm)),
      r_onestage = cor(pm, yte), time_onestage_exact = t_ex,
      time_yqfe = if (is.null(extra$tfe)) NA else extra$tfe)
  for (Q in QS) {
    t0 <- tnow(); lab <- make_groups("YQ", ytr, Q); S1 <- yq_fe_stage1(Ztr, ytr, lab); t1 <- tnow() - t0
    t0 <- tnow(); f <- fit_yq_fen(S1, k = 500); pr <- predict_lin(f, Zte); t2 <- tnow() - t0
    ex <- list()
    if (Q == 1000) { t0 <- tnow(); fe <- fit_yq_fe(S1); pfe <- predict_lin(fe, Zte); ex$tfe <- t1 + tnow() - t0 }
    rec("k500", Q, 500, pr, t1, t2, f, ex)
    if (Q == 1000) {
      for (k in c(0, 250, 1000)) {
        t0 <- tnow(); fk <- fit_yq_fen(S1, k = k); pk <- predict_lin(fk, Zte); rec(paste0("k", k), Q, k, pk, t1, tnow() - t0, fk) }
      t0 <- tnow(); fr <- fit_yq_fen(S1, k = 500, delta = fe$delta); prr <- predict_lin(fr, Zte)
      rec("k500_REMLdelta", Q, 500, prr, t1, tnow() - t0, fr)
      rows[[length(rows)]]$max_abs_diff_vs_yqfe <- max(abs(prr - pfe))
    }
    rm(S1, f); gc(verbose = FALSE)
  }
  out <- do.call(rbind, lapply(rows, function(r) { if (is.null(r$max_abs_diff_vs_yqfe)) r$max_abs_diff_vs_yqfe <- NA; r }))
  out
}
if (part == "sim") {
  out <- file.path(resd, "yq_fen_sim_results.csv")
  rows <- if (file.exists(out)) read.csv(out) else NULL
  done <- if (is.null(rows)) character(0) else unique(paste(rows$scenario, rows$rep))
  for (si in seq_along(SIM$SCEN)) for (rep in seq_len(SIM$NREP) - 1) {
    sc <- SIM$SCEN[si]; if (paste(sc, rep) %in% done) next
    message(sprintf("%s rep %d", sc, rep + 1)); t0 <- tnow()
    d <- make_sim_rep(si, rep)
    res <- eval_fen(d$Z[d$tr, ], d$y[d$tr], d$Z[d$te, ], d$y[d$te], d$g[d$te])
    res$scenario <- sc; res$rep <- rep; res$n_train <- length(d$tr); res$p <- SIM$P
    rows <- rbind(rows, res); write.csv(rows, out, row.names = FALSE)
    rm(d); gc(verbose = FALSE); message(sprintf("  %.0f s", tnow() - t0))
  }
}
if (part == "real") {
  D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds"))
  out <- file.path(resd, "yq_fen_real_results.csv")
  rows <- if (file.exists(out)) read.csv(out) else NULL
  done <- if (is.null(rows)) integer(0) else unique(rows$rep)
  for (rep in 0:4) {
    if (rep %in% done) next
    pt <- readRDS(file.path(resd, sprintf("real_partition_%d.rds", rep)))
    message(sprintf("real partition %d", rep + 1)); t0 <- tnow()
    res <- eval_fen(D$Z[pt$tr, ], D$y[pt$tr], D$Z[pt$te, ], D$y[pt$te], NULL)
    res$scenario <- "real_wheat"; res$rep <- rep; res$n_train <- length(pt$tr); res$p <- ncol(D$Z)
    rows <- rbind(rows, res); write.csv(rows, out, row.names = FALSE)
    message(sprintf("  %.0f s", tnow() - t0)); gc(verbose = FALSE)
  }
}
