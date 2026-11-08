###############################################################################
# 20_revision_eval.R -- revision analyses at Q = 1,000 (v12):
#  * agreement with the one-stage GBLUP (largest and RMS difference, overlap of
#    the top 10% selected) for YQ-FEN, YQ, YQ-CAL, YQ-CVCAL, NYS, SUB;
#  * equal-budget calibration: YQ, NYS and the one-stage GBLUP all fitted to the
#    same 90% and recalibrated on the same 10% (CAL); YQ with cross-fitted
#    calibration (CVCAL, 5 folds, all data used for the fit);
#  * YQ-FEN: moment estimate of delta (identifiability, truncation), solver error
#    (REML delta) and sensitivity to k, h and the landmark seed (original data);
#  * four additional simulated scenarios (h2 = 0.2, h2 = 0.8, full-sib families,
#    new families in the testing set) and leave-one-data-set-out validation.
# Usage: Rscript 20_revision_eval.R [orig_sim|orig_real|new_sim|lodo]   (resumable)
# Outputs: outputs/results/rev_<part>.csv and rev_sens_<part>.csv
###############################################################################
source("00_functions.R")
resd <- file.path("..", "outputs", "results")
part <- commandArgs(TRUE)[1]; if (is.na(part)) part <- "orig_sim"
Q <- 1000
eval_rev <- function(Ztr, ytr, Zte, yte, gte, seed, do_fe = FALSE, sens = FALSE) {
  rows <- list(); sens_rows <- list()
  t0 <- tnow(); fm <- fit_marker_gblup(Ztr, ytr); pm <- predict_marker_gblup(fm, Zte); tex <- tnow() - t0
  rec <- function(method, pr, t, delta = NA, extra = list()) rows[[length(rows) + 1]] <<- data.frame(
    method = method, t(metrics(pr, yte, gte)), t(agree(pr, pm)), time = t, delta = delta,
    truncated = if (is.null(extra$tr)) NA else extra$tr, rel_det = if (is.null(extra$rd)) NA else extra$rd,
    iter = if (is.null(extra$it)) NA else extra$it)
  rec("EXACT", pm, tex, fm$delta)
  t0 <- tnow(); lab <- make_groups("YQ", ytr, Q); S1 <- yq_fe_stage1(Ztr, ytr, lab); ts1 <- tnow() - t0
  if (do_fe) { t0 <- tnow(); f <- fit_yq_fe(S1); rec("YQFE", predict_lin(f, Zte), ts1 + tnow() - t0, f$delta) }
  t0 <- tnow(); f <- fit_yq_fen(S1, k = 500, h = 200); rec("YQFEN", predict_lin(f, Zte), ts1 + tnow() - t0, f$delta,
    list(tr = f$truncated, rd = f$rel_det, it = f$iter))
  t0 <- tnow(); f <- fit_yq_fen(S1, k = 500, delta = fm$delta); rec("YQFEN_REMLdelta", predict_lin(f, Zte), ts1 + tnow() - t0, f$delta, list(it = f$iter))
  if (sens) {
    Tt <- yq_fe_totals(S1)
    for (sd in 1:5) for (k in c(250, 500, 1000, 2000)) {
      if (sd > 1 && k != 500) next
      t0 <- tnow(); L <- nystrom_markers(Tt$Tc, k, sd); tn <- tnow() - t0
      for (h in c(0, 50, 100, 200, 400)) {
        if (sd > 1 && h != 200) next
        hh <- min(h, ncol(L$U)); t0 <- tnow()
        mo <- moments_delta(Tt, if (hh > 0) L$U[, seq_len(hh), drop = FALSE] else NULL)
        s <- pcg_nystrom(Tt$Tc, Tt$sc, Tt$p * mo$delta, L, tol = 1e-10)
        pr <- drop(Tt$ybar - sum(Tt$zbar * s$x) + Zte %*% s$x)
        sens_rows[[length(sens_rows) + 1]] <- data.frame(k = k, h = h, seed = sd, delta = mo$delta, delta_reml = fm$delta,
          truncated = mo$truncated, rel_det = mo$rel_det, iter = s$iter, time_nystrom = tn, time_rest = tnow() - t0,
          r = cor(pr, yte), r_onestage = cor(pm, yte), t(agree(pr, pm)))
      }
    }
  }
  rm(S1); gc(verbose = FALSE)
  ## YQ (all data), YQ-CVCAL (all data + cross-fitted calibration)
  t0 <- tnow(); S <- summarise_groups(Ztr, ytr, lab); f <- fit_wgblup(S$Zbar, S$ybar, S$m); rec("YQ", predict_wgblup(f, Zte), tnow() - t0, f$delta)
  t0 <- tnow(); cv <- fit_yq_cvcal(Ztr, ytr, Q, 5, seed); rec("YQ_CVCAL", cv$a + cv$b * predict_wgblup(cv$f, Zte), tnow() - t0, cv$f$delta)
  ## equal budget: 90% fit, 10% calibration (same split as evaluate_rep)
  n <- length(ytr); set.seed(seed); perm <- sample.int(n); nv <- floor(0.1 * n); iv <- perm[seq_len(nv)]; it <- perm[-seq_len(nv)]
  cal <- function(pv, pt) { bb <- cov(pv, ytr[iv]) / var(pv); aa <- mean(ytr[iv]) - bb * mean(pv); aa + bb * pt }
  t0 <- tnow(); S <- summarise_groups(Ztr[it, ], ytr[it], make_groups("YQ", ytr[it], Q)); f <- fit_wgblup(S$Zbar, S$ybar, S$m)
  rec("YQ_CAL", cal(predict_wgblup(f, Ztr[iv, ]), predict_wgblup(f, Zte)), tnow() - t0, f$delta)
  t0 <- tnow(); ny <- fit_nystrom(Ztr, ytr, Q, seed + Q); rec("NYS", predict_nystrom(ny, Zte)$pred, tnow() - t0, ny$delta)
  t0 <- tnow(); ny <- fit_nystrom(Ztr[it, ], ytr[it], Q, seed + Q)
  rec("NYS_CAL", cal(predict_nystrom(ny, Ztr[iv, ])$pred, predict_nystrom(ny, Zte)$pred), tnow() - t0, ny$delta); rm(ny)
  t0 <- tnow(); f9 <- fit_marker_gblup(Ztr[it, ], ytr[it])
  rec("EXACT_CAL", cal(predict_marker_gblup(f9, Ztr[iv, ]), predict_marker_gblup(f9, Zte)), tnow() - t0, f9$delta); rm(f9)
  set.seed(seed + Q); idx <- sample.int(n, Q); t0 <- tnow(); f <- fit_wgblup(Ztr[idx, ], ytr[idx], rep(1, Q))
  rec("SUB", predict_wgblup(f, Zte), tnow() - t0, f$delta)
  gc(verbose = FALSE)
  list(res = do.call(rbind, rows), sens = if (length(sens_rows)) do.call(rbind, sens_rows) else NULL)
}
save_part <- function(res, key, out, outs) {
  r <- res$res; for (nm in names(key)) r[[nm]] <- key[[nm]]
  old <- if (file.exists(out)) read.csv(out) else NULL; write.csv(rbind(old, r), out, row.names = FALSE)
  if (!is.null(res$sens)) { s <- res$sens; for (nm in names(key)) s[[nm]] <- key[[nm]]
    old <- if (file.exists(outs)) read.csv(outs) else NULL; write.csv(rbind(old, s), outs, row.names = FALSE) }
}
out <- file.path(resd, sprintf("rev_%s.csv", part)); outs <- file.path(resd, sprintf("rev_sens_%s.csv", part))
done <- if (file.exists(out)) unique(paste(read.csv(out)$scenario, read.csv(out)$rep)) else character(0)
if (part == "orig_sim") for (si in seq_along(SIM$SCEN)) for (rep in 0:4) {
  sc <- SIM$SCEN[si]; if (paste(sc, rep) %in% done) next; message(sc, " ", rep); t0 <- tnow()
  d <- make_sim_rep(si, rep)
  res <- eval_rev(d$Z[d$tr, ], d$y[d$tr], d$Z[d$te, ], d$y[d$te], d$g[d$te], d$seed, sens = TRUE)
  save_part(res, list(scenario = sc, rep = rep), out, outs); message(sprintf("  %.0f s", tnow() - t0)) }
if (part == "new_sim") for (si in seq_along(SIM2$SCEN)) for (rep in 0:4) {
  sc <- SIM2$SCEN[si]; if (paste(sc, rep) %in% done) next; message(sc, " ", rep); t0 <- tnow()
  d <- make_sim2_rep(si, rep)
  res <- eval_rev(d$Z[d$tr, ], d$y[d$tr], d$Z[d$te, ], d$y[d$te], d$g[d$te], d$seed, do_fe = TRUE)
  save_part(res, list(scenario = sc, rep = rep), out, outs); message(sprintf("  %.0f s", tnow() - t0)) }
if (part %in% c("orig_real", "lodo")) {
  D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds"))
  if (part == "orig_real") for (rep in 0:4) {
    if (paste("real_wheat", rep) %in% done) next; message("partition ", rep); t0 <- tnow()
    pt <- readRDS(file.path(resd, sprintf("real_partition_%d.rds", rep)))
    res <- eval_rev(D$Z[pt$tr, ], D$y[pt$tr], D$Z[pt$te, ], D$y[pt$te], NULL, pt$seed, sens = TRUE)
    save_part(res, list(scenario = "real_wheat", rep = rep), out, outs); message(sprintf("  %.0f s", tnow() - t0)) }
  if (part == "lodo") for (w in sort(unique(D$env))) {
    if (paste("lodo", w) %in% done) next; message("leave out ", w); t0 <- tnow()
    te <- which(D$env == w); tr <- which(D$env != w)
    res <- eval_rev(D$Z[tr, ], D$y[tr], D$Z[te, ], D$y[te], NULL, 2025 + match(w, sort(unique(D$env))), do_fe = TRUE)
    save_part(res, list(scenario = "lodo", rep = w), out, outs); message(sprintf("  %.0f s", tnow() - t0)) }
}
