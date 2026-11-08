###############################################################################
# 25_yq_fe_q.R -- YQ-FE-Q: YQ-FE written as a GBLUP on the Q group means
# (mixed-model equations of dimension Q + 1). Numerical examples of
# Supplementary Section S10:
#   part toy  : toy data of 13_yq_fe_toy.R (n = 12, p = 4, Q = 3), delta = 1 and REML
#   part sim  : scenarios S1-S3 (replicate 1), Q = 1,000; time versus Q in S1;
#               low-rank (approximate) versions YQ-FE-Q(k)
#   part kern : Gaussian kernel (non-additive) with the kernel form, S1 subset
#   part real : wheat data, partition 1, Q = 1,000; marker and kernel forms;
#               REML profiles of the one-stage and Q-dimensional criteria
# Usage: Rscript 25_yq_fe_q.R [toy|sim|kern|real]   (OPENBLAS_NUM_THREADS=1)
###############################################################################
source("00_functions.R")
part <- commandArgs(TRUE)[1]; if (is.na(part)) part <- "toy"
RES <- file.path("..", "outputs", "results"); TAB <- file.path("..", "outputs", "tables")
tm <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; list(v = v, t = proc.time()[["elapsed"]] - t0) }

if (part == "toy") {
  sink(file.path(TAB, "yq_fe_q_toy.txt"), split = TRUE)
  Z <- rbind(c(-0.01, 1.05, 0.74, 0.72), c(1.62, -1.21, -0.63, -1.32), c(-0.11, 1.00, -0.02, 0.50),
             c(-1.91, 0.15, -0.91, 1.78), c(0.89, 0.95, -0.06, 0.61), c(0.66, -0.34, -0.50, -0.11),
             c(-0.61, -0.59, -0.28, -0.73), c(0.77, -1.60, 0.82, -0.63), c(-0.55, -1.35, -0.14, -0.25),
             c(0.19, -0.53, 0.09, 1.82), c(0.41, -0.57, 0.95, -0.13), c(0.59, 0.61, -0.39, -1.93))
  y <- c(-0.62, 2.37, -1.01, -1.82, 0.99, 0.16, 0.68, 3.27, 1.12, 1.46, 1.97, 1.81)
  z0 <- rbind(c(0.50, -0.80, 0.30, 0.10), c(-1.00, 0.40, -0.20, 0.90))   # two new individuals
  n <- 12; p <- 4; delta <- 1; lab <- make_groups("YQ", y, 3)
  S1 <- yq_fe_stage1(Z, y, lab); P <- yq_fe_q_prep(S1); lam <- p * delta
  X1 <- cbind(1, Z); th1 <- solve(crossprod(X1) + diag(c(0, rep(lam, p))), crossprod(X1, y))
  cat("delta = 1 (lambda = p delta = 4)\n\none-stage (mu, b):", round(th1, 5), "\n")
  pt <- yq_fe_q_parts(S1, P, delta)
  mw <- solve(S1$W + lam * diag(p), S1$sW)
  cat("\nm_w = (W + p delta I)^-1 s_W:", round(mw, 5), "\n")
  cat("Zbar m_w:", round(drop(S1$Zbar %*% mw), 5), "\n")
  cat("ytilde = ybar - Zbar m_w:", round(pt$yt, 5), "\n")
  cat("\nK_G(delta) = delta Zbar (W + p delta I)^-1 Zbar':\n"); print(round(pt$K, 5))
  cat("\nG_bb = Zbar Zbar'/p (kernel of YQ):\n"); print(round(tcrossprod(S1$Zbar) / p, 5))
  Vb <- pt$K + diag(delta / S1$m); cat("\nV_b = K_G + delta D^-1:\n"); print(round(Vb, 5))
  # mixed-model equations of dimension Q + 1 in (mu, u)
  D <- diag(S1$m); Ki <- solve(pt$K)
  LHS <- rbind(c(sum(S1$m), S1$m), cbind(S1$m, D + delta * Ki)); RHS <- c(sum(S1$m * pt$yt), S1$m * pt$yt)
  cat("\nmixed-model equations (Q + 1 = 4):\nLHS =\n"); print(round(LHS, 5)); cat("RHS =", round(RHS, 5), "\n")
  sol <- solve(LHS, RHS); cat("solution (mu, u):", round(sol, 5), "\n")
  f <- fit_yq_fe_q(S1, delta = delta, P = P)
  cat("alpha = V_b^-1 (ytilde - 1 mu):", round(f$alpha, 5), "  K_G alpha =", round(drop(pt$K %*% f$alpha), 5), "\n")
  cat("\nYQ-FE-Q (mu, b) recovered:", round(c(f$mu, f$b), 5), "\n")
  cat("max |difference| with one-stage:", signif(max(abs(c(f$mu, f$b) - th1)), 3), "\n")
  K0 <- delta * z0 %*% solve(S1$W + lam * diag(p), t(S1$Zbar))
  cat("\nnew individuals: z0'm_w =", round(drop(z0 %*% mw), 5), "; K_0G =\n"); print(round(K0, 5))
  pr <- f$mu + drop(z0 %*% mw) + drop(K0 %*% f$alpha)
  cat("predictions YQ-FE-Q:", round(pr, 5), "  one-stage:", round(drop(cbind(1, z0) %*% th1), 5), "\n")
  # kernel form, REML
  G <- tcrossprod(Z) / p
  k1 <- fit_kernel_gblup(G, y, tcrossprod(z0, Z) / p); kq <- fit_yq_fe_q_kernel(G, y, lab, tcrossprod(z0, Z) / p)
  fr <- fit_yq_fe_q(S1)
  cat(sprintf("\nREML: delta one-stage = %.6f; YQ-FE-Q (marker form) = %.6f; YQ-FE-Q (kernel form) = %.6f\n", k1$delta, fr$delta, kq$delta))
  cat(sprintf("criterion one-stage minus Q-dimensional = %.6f ; log|D| = %.6f\n", k1$crit - kq$crit, sum(log(S1$m))))
  cat("predictions at REML: one-stage", round(k1$pred, 5), "; YQ-FE-Q", round(kq$pred, 5), "\n")
  sink()
}

if (part == "sim") {
  out <- NULL; prof <- NULL
  for (si in 1:3) {
    d <- make_sim_rep(si, 0); Ztr <- d$Z[d$tr, ]; ytr <- d$y[d$tr]; Zte <- d$Z[d$te, ]; yte <- d$y[d$te]; gte <- d$g[d$te]
    one <- tm(fit_marker_gblup(Ztr, ytr)); po <- predict_marker_gblup(one$v, Zte)
    Qs <- if (si == 1) c(250, 1000, 2000, 3000) else 1000
    for (Q in Qs) {
      lab <- make_groups("YQ", ytr, Q); s1 <- tm(yq_fe_stage1(Ztr, ytr, lab)); S1 <- s1$v
      fe <- tm(fit_yq_fe(S1)); pfe <- predict_lin(fe$v, Zte)
      fq <- fit_yq_fe_q(S1); pq <- predict_lin(fq, Zte)
      fq1 <- fit_yq_fe_q(S1, delta = one$v$delta, P = NULL); pq1 <- predict_lin(fq1, Zte)
      t_eval <- (fq$t_total - fq$t_prep) / max(fq$neval + 1, 1)
      out <- rbind(out, data.frame(data = SIM$SCEN[si], Q = Q, n = length(ytr), p = ncol(Ztr),
        delta_one = one$v$delta, delta_fe = fe$v$delta, delta_q = fq$delta,
        r_one = cor(po, yte), r_q = cor(pq, yte), rg_one = cor(po, gte), rg_q = cor(pq, gte),
        maxdiff_samedelta = max(abs(pq1 - po)), maxdiff_reml = max(abs(pq - po)),
        t_one = one$t, t_stage1 = s1$t, t_fe = fe$t, t_q_prep = fq$t_prep, t_q_total = fq$t_total,
        q_neval = fq$neval, t_q_per_eval = t_eval))
      print(tail(out, 1))
      rm(S1); gc(verbose = FALSE)
    }
  }
  write.csv(out, file.path(RES, "yq_fe_q_sim.csv"), row.names = FALSE)
}

if (part == "lowrank") {   # approximate YQ-FE-Q(k): W replaced by its k leading eigen-components
  d <- make_sim_rep(1, 0); Ztr <- d$Z[d$tr, ]; ytr <- d$y[d$tr]; Zte <- d$Z[d$te, ]; yte <- d$y[d$te]; gte <- d$g[d$te]
  one <- fit_marker_gblup(Ztr, ytr); po <- predict_marker_gblup(one, Zte)
  lab <- make_groups("YQ", ytr, 1000); S1 <- yq_fe_stage1(Ztr, ytr, lab); P <- yq_fe_q_prep(S1); prof <- NULL
  for (k in c(0, 100, 250, 500, 1000, 2000, 2500, P$r)) {
    kk <- min(k, P$r); Pk <- P; ix <- seq_len(kk)
    Pk$V <- P$V[, ix, drop = FALSE]; Pk$om <- P$om[ix]; Pk$Ab <- P$Ab[, ix, drop = FALSE]; Pk$vs <- P$vs[ix]; Pk$r <- kk
    for (rule in c("one-stage REML delta", "REML of the truncated criterion")) {
      a <- tm(fit_yq_fe_q(S1, P = Pk, delta = if (rule == "one-stage REML delta") one$delta else NULL)); pk <- predict_lin(a$v, Zte)
      prof <- rbind(prof, data.frame(k = kk, share_W = sum(P$om[ix]) / sum(P$om), delta_rule = rule, delta = a$v$delta,
        pct_r = 100 * cor(pk, yte) / cor(po, yte), rg = cor(pk, gte), slope = unname(coef(lm(yte ~ pk))[2]),
        maxdiff = max(abs(pk - po)), time = a$t))
      print(tail(prof, 1))
    }
  }
  write.csv(prof, file.path(RES, "yq_fe_q_lowrank.csv"), row.names = FALSE)
}

if (part == "kern") {   # Gaussian kernel: the Q-dimensional form holds for any kernel
  d <- make_sim_rep(1, 0); set.seed(11)
  tr <- d$tr[1:2000]; te <- d$te[1:500]; Ztr <- d$Z[tr, ]; Zte <- d$Z[te, ]
  sq <- function(A, B) outer(rowSums(A^2), rowSums(B^2), "+") - 2 * tcrossprod(A, B)
  D2 <- sq(Ztr, Ztr); h <- median(D2[upper.tri(D2)]); G <- exp(-D2 / h); G0 <- exp(-sq(Zte, Ztr) / h); rm(D2)
  ytr <- d$y[tr]; yte <- d$y[te]; out <- NULL
  one <- tm(fit_kernel_gblup(G, ytr, G0))
  for (Q in c(100, 200, 500, 1000)) {
    lab <- make_groups("YQ", ytr, Q)
    kq <- tm(fit_yq_fe_q_kernel(G, ytr, lab, G0))
    # YQ with the Gaussian kernel of the group means (within-group information ignored)
    A <- as.matrix(Matrix::sparseMatrix(i = as.integer(factor(lab)), j = seq_along(lab), x = 1 / tabulate(as.integer(factor(lab)))[as.integer(factor(lab))]))
    m <- tabulate(as.integer(factor(lab))); Kb <- A %*% G %*% t(A); K0b <- G0 %*% t(A); yb <- drop(A %*% ytr)
    cyq <- function(dl) { V <- Kb + diag(dl / m); L <- chol(V); Vi <- chol2inv(L); mu <- sum(Vi %*% yb) / sum(Vi)
      r <- yb - mu; s2 <- drop(r %*% Vi %*% r) / (Q - 1); (Q - 1) * log(s2) + 2 * sum(log(diag(L))) + log(sum(Vi)) }
    dl <- exp(optimize(function(l) cyq(exp(l)), log(c(1e-4, 1e4)))$minimum)
    V <- Kb + diag(dl / m); Vi <- solve(V); mu <- sum(Vi %*% yb) / sum(Vi); pyq <- mu + drop(K0b %*% Vi %*% (yb - mu))
    out <- rbind(out, data.frame(kernel = "Gaussian", n = 2000, Q = Q, delta_one = one$v$delta, delta_q = kq$v$delta,
      r_one = cor(one$v$pred, yte), r_q = cor(kq$v$pred, yte), r_yq = cor(pyq, yte),
      maxdiff = max(abs(kq$v$pred - one$v$pred)), crit_diff = one$v$crit - kq$v$crit, logdetD = sum(log(m)),
      t_one = one$t, t_q = kq$t))
    print(tail(out, 1))
  }
  write.csv(out, file.path(RES, "yq_fe_q_kernel.csv"), row.names = FALSE)
}

if (part == "real") {
  D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds")); pr <- readRDS(file.path(RES, "real_partition_0.rds"))
  Ztr <- D$Z[pr$tr, ]; ytr <- D$y[pr$tr]; Zte <- D$Z[pr$te, ]; yte <- D$y[pr$te]; Q <- 1000; p <- ncol(Ztr)
  one <- tm(fit_marker_gblup(Ztr, ytr)); po <- predict_marker_gblup(one$v, Zte)
  lab <- make_groups("YQ", ytr, Q); s1 <- tm(yq_fe_stage1(Ztr, ytr, lab)); S1 <- s1$v
  fe <- tm(fit_yq_fe(S1)); pfe <- predict_lin(fe$v, Zte)
  fq <- fit_yq_fe_q(S1); pq <- predict_lin(fq, Zte)
  fq1 <- fit_yq_fe_q(S1, delta = one$v$delta); pq1 <- predict_lin(fq1, Zte)
  rm(S1); gc(verbose = FALSE)
  G <- tcrossprod(Ztr) / p; G0 <- tcrossprod(Zte, Ztr) / p
  kq <- tm(fit_yq_fe_q_kernel(G, ytr, lab, G0)); gc(verbose = FALSE)
  k1 <- tm(fit_kernel_gblup(G, ytr, G0))
  out <- data.frame(data = "wheat, partition 1", Q = Q, n = length(ytr), p = p,
    delta_one = one$v$delta, delta_fe = fe$v$delta, delta_q = fq$delta, delta_q_kernel = kq$v$delta, delta_kernel_one = k1$v$delta,
    r_one = cor(po, yte), r_q = cor(pq, yte), maxdiff_samedelta = max(abs(pq1 - po)), maxdiff_reml = max(abs(pq - po)),
    maxdiff_kernel = max(abs(kq$v$pred - k1$v$pred)), crit_diff = k1$v$crit - kq$v$crit, logdetD = sum(log(table(lab))),
    t_one = one$t, t_stage1 = s1$t, t_fe = fe$t, t_q_prep = fq$t_prep, t_q_total = fq$t_total, q_neval = fq$neval,
    t_q_kernel = kq$t, t_q_kernel_prep = kq$v$t_prep, t_kernel_one = k1$t)
  print(out); write.csv(out, file.path(RES, "yq_fe_q_real.csv"), row.names = FALSE)
  # REML profiles (kernel forms) for Supplementary Figure S6B, and predictions for S6A
  grid <- exp(seq(log(0.5), log(8), length.out = 25))
  S <- list(m = as.vector(table(lab)))
  pf <- data.frame(delta = grid, crit_one = sapply(grid, k1$v$crit_fun))
  # Q-dimensional criterion through the marker form (same function used in fitting)
  S1 <- yq_fe_stage1(Ztr, ytr, lab); P <- yq_fe_q_prep(S1)
  pf$crit_q <- sapply(grid, function(dl) yq_fe_q_crit(S1, P, dl)$crit)
  write.csv(pf, file.path(RES, "yq_fe_q_real_profile.csv"), row.names = FALSE)
  write.csv(data.frame(one = po, yqfeq = pq, y = yte), file.path(RES, "yq_fe_q_real_pred.csv"), row.names = FALSE)
}

if (part == "fig") {   # Supplementary Figure S6 (YQ-FE-Q)
  suppressMessages({ library(ggplot2); library(patchwork) })
  th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
  pr <- read.csv(file.path(RES, "yq_fe_q_real_pred.csv")); pf <- read.csv(file.path(RES, "yq_fe_q_real_profile.csv"))
  rl <- read.csv(file.path(RES, "yq_fe_q_real.csv")); sm <- read.csv(file.path(RES, "yq_fe_q_sim.csv")); lr <- read.csv(file.path(RES, "yq_fe_q_lowrank.csv"))
  A <- ggplot(pr, aes(one, yqfeq)) + geom_abline(linetype = 2, colour = "grey50") + geom_point(size = 0.6, alpha = 0.5, colour = "#1f5fa8") +
    labs(x = "One-stage GBLUP prediction", y = "YQ-FE-Q prediction", title = "A  Wheat data, Q = 1,000") + th
  pf$crit_q_shift <- pf$crit_q + rl$logdetD
  pl <- rbind(data.frame(delta = pf$delta, crit = pf$crit_one, form = "One-stage GBLUP"),
              data.frame(delta = pf$delta, crit = pf$crit_q_shift, form = "YQ-FE-Q (S10.13) + log|D|"))
  B <- ggplot(pl, aes(delta, crit, colour = form, linetype = form)) + geom_line(linewidth = 0.8) +
    geom_vline(xintercept = rl$delta_one, colour = "grey40", linetype = 3) + scale_x_log10() +
    scale_colour_manual(values = c("#1f5fa8", "#e07b00")) + scale_linetype_manual(values = c(1, 2)) +
    labs(x = expression(delta ~ "(log scale)"), y = "-2 REML log-likelihood", colour = NULL, linetype = NULL, title = "B  REML criteria") +
    th + theme(legend.position = c(0.6, 0.85), legend.background = element_blank())
  s1 <- sm[sm$data == "S1_polygenic", ]
  tc <- rbind(data.frame(Q = s1$Q, t = s1$t_fe, method = "YQ-FE Stage 2 (p x p)"),
              data.frame(Q = s1$Q, t = s1$t_q_total, method = "YQ-FE-Q Stage 2 (Q + 1)"))
  C <- ggplot(tc, aes(Q, t, colour = method)) + geom_line() + geom_point() + scale_colour_manual(values = c("#1f5fa8", "#e07b00")) +
    labs(x = "Number of groups Q", y = "Stage 2 time (s, one core)", colour = NULL, title = "C  Cost versus Q (S1)") +
    th + theme(legend.position = c(0.4, 0.85), legend.background = element_blank())
  l2 <- lr[lr$delta_rule == "one-stage REML delta", ]
  D <- ggplot(l2, aes(100 * share_W, pct_r)) + geom_line(colour = "#1f5fa8") + geom_point(colour = "#1f5fa8") +
    geom_hline(yintercept = 100, linetype = 2, colour = "grey50") +
    labs(x = "Share of trace(W) retained (%)", y = "% of one-stage correlation", title = "D  Low-rank YQ-FE-Q(k)") + th
  g <- (A | B) / (C | D)
  ggsave(file.path("..", "outputs", "figures", "FigureS10_yq_fe_q.png"), g, width = 9, height = 7.2, dpi = 350)
}
