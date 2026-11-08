###############################################################################
# 06_equivalence_checks.R
# (A) YQ with Q = n reproduces the one-stage (full) GBLUP.
# (B) The YQ fit is the Stage 2 of a two-stage analysis: with replicated lines
#     (groups = lines), one-stage GBLUP on plot records (lme4) and two-stage
#     GBLUP (Stage 1 line BLUEs, Stage 2 weighted GBLUP, the same function used
#     for YQ) give the same GEBVs (Endelman, 2023; Piepho et al., 2012).
# (C) The weighted Stage 2 fitted with lme4 (REML) and with BGLR (Bayesian ridge
#     regression, weights = sqrt(m)) on the real data (YQ, Q = 1,000).
###############################################################################
source("00_functions.R")
res <- list()
## (A) -----------------------------------------------------------------------
set.seed(11)
d <- simulate_data("S1_polygenic", 2500, 1500); Z <- std_markers(d$X)
tr <- 1:2000; te <- 2001:2500
f1 <- fit_wgblup(Z[tr, ], d$y[tr], rep(1, 2000)); p1 <- predict_wgblup(f1, Z[te, ])
S <- summarise_groups(Z[tr, ], d$y[tr], make_groups("YQ", d$y[tr], 2000))
f2 <- fit_wgblup(S$Zbar, S$ybar, S$m); p2 <- predict_wgblup(f2, Z[te, ])
## (D) exact marker-space GBLUP (p x p) versus the n x n kernel GBLUP, same data as (A)
t0 <- tnow(); fm <- fit_marker_gblup(Z[tr, ], d$y[tr], block = 500); pm <- predict_marker_gblup(fm, Z[te, ]); tm <- tnow() - t0
res$D <- data.frame(check = "D: exact marker-space GBLUP (p x p) vs one-stage GBLUP (n x n)",
  quantity = c("delta n x n kernel", "delta marker space", "max |difference| in predictions (marker space)"),
  value = c(f1$delta, fm$delta, max(abs(p1 - pm))))
res$A <- data.frame(check = "A: YQ with Q = n vs one-stage GBLUP",
  quantity = c("delta one-stage", "delta YQ (Q = n)", "max |difference| in predictions", "correlation of predictions"),
  value = c(f1$delta, f2$delta, max(abs(p1 - p2)), cor(p1, p2)))
## (B) -----------------------------------------------------------------------
set.seed(12)
L <- 800; p <- 1500
d <- simulate_data("S1_polygenic", L + 200, p, h2 = 1)
Zl <- std_markers(d$X); g <- d$g / sd(d$g) * sqrt(0.3)      # line genetic values
trL <- 1:L; teL <- L + 1:200
reps <- sample(1:4, L, replace = TRUE)                       # unbalanced replication
line <- rep(trL, reps)
yrec <- 10 + g[line] + rnorm(length(line), sd = sqrt(0.7))   # plot records
## one-stage: y = mu + W g + e, g ~ N(0, s2g K), K = Z Z'/p = R'R (lme4, records)
K <- tcrossprod(Zl[trL, ]) / p
R <- chol(K + diag(1e-6, L))                                 # K = R'R
dat <- data.frame(y = yrec, line = factor(line, levels = trL))
lf <- lFormula(y ~ 1 + (1 | line), data = dat, REML = TRUE, control = lme4_ctrl)
lf$reTrms$Zt <- as(R %*% lf$reTrms$Zt, "CsparseMatrix")       # Zt = R W'
devf <- do.call(mkLmerDevfun, lf); opt <- optimizeLmer(devf)
m1 <- mkMerMod(environment(devf), opt, lf$reTrms, fr = lf$fr)
s2e1 <- sigma(m1)^2; s2g1 <- s2e1 * getME(m1, "theta")^2
gebv1 <- drop(crossprod(R, getME(m1, "b")))                   # g = R'u
## two-stage: Stage 1 BLUEs of line effects (fixed-effects model), Stage 2 weighted GBLUP
st1 <- lm(y ~ 0 + line, data = dat)
blue <- unname(coef(st1)); w <- as.vector(table(dat$line))    # Var(BLUE) = s2e / r_i
f2s <- fit_wgblup(Zl[trL, ], blue, w)
cf <- coef_wgblup(f2s); gebv2 <- drop(Zl[trL, ] %*% cf$b)
cf_plug <- coef_wgblup(f2s, s2e1 / s2g1); gebv2p <- drop(Zl[trL, ] %*% cf_plug$b)
new1 <- drop(Zl[teL, ] %*% coef_wgblup(f2s, s2e1 / s2g1)$b)
## one-stage prediction of new lines from the one-stage marker effects
alpha1 <- solve(K + diag(1e-6, L), gebv1); b1 <- drop(crossprod(Zl[trL, ], alpha1)) / p
new1s <- drop(Zl[teL, ] %*% b1)
res$B <- data.frame(check = "B: one-stage vs two-stage (replicated lines)",
  quantity = c("s2g one-stage", "s2e one-stage", "s2g two-stage", "s2e two-stage",
               "cor(GEBV one-stage, two-stage REML)", "max |diff| GEBV, two-stage REML",
               "max |diff| GEBV, two-stage with one-stage variance ratio",
               "max |diff| prediction of new lines (same variance ratio)",
               "cor(prediction, true g) new lines, one-stage", "cor(prediction, true g) new lines, two-stage"),
  value = c(s2g1, s2e1, f2s$s2g, f2s$s2e, cor(gebv1, gebv2), max(abs(gebv1 - gebv2)),
            max(abs(gebv1 - gebv2p)), max(abs(new1 - new1s)),
            cor(new1s, g[teL]), cor(predict_wgblup(f2s, Zl[teL, ]), g[teL])))
eq_scatter <- data.frame(one_stage = gebv1, two_stage = gebv2, reps = w)
write.csv(eq_scatter, file.path("..", "outputs", "results", "two_stage_scatter.csv"), row.names = FALSE)
## (C) -----------------------------------------------------------------------
if (requireNamespace("BGLR", quietly = TRUE)) {
  D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds"))
  pt <- readRDS(file.path("..", "outputs", "results", "real_partition_0.rds"))
  S <- summarise_groups(D$Z[pt$tr, ], D$y[pt$tr], make_groups("YQ", D$y[pt$tr], 1000))
  fl <- fit_wgblup(S$Zbar, S$ybar, S$m); pl <- predict_wgblup(fl, D$Z[pt$te, ])
  set.seed(1); t0 <- tnow()
  fb <- BGLR::BGLR(y = S$ybar, ETA = list(list(X = S$Zbar, model = "BRR")), weights = sqrt(S$m),
                   nIter = 6000, burnIn = 1000, thin = 5, verbose = FALSE,
                   saveAt = file.path(tempdir(), "yq_"))
  tb <- tnow() - t0
  pb <- drop(fb$mu + D$Z[pt$te, ] %*% fb$ETA[[1]]$b)
  res$C <- data.frame(check = "C: lme4 (REML) vs BGLR (BRR) for YQ, real data, Q = 1000",
    quantity = c("r test lme4", "r test BGLR", "cor(lme4, BGLR predictions)",
                 "SD ratio BGLR/lme4 predictions", "BGLR time (s)"),
    value = c(cor(pl, D$y[pt$te]), cor(pb, D$y[pt$te]), cor(pl, pb), sd(pb) / sd(pl), tb))
}
out <- do.call(rbind, res[intersect(c("A", "D", "B", "C"), names(res))])
write.csv(out, file.path("..", "outputs", "results", "equivalence_checks.csv"), row.names = FALSE)
print(out)
