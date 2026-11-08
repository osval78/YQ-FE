###############################################################################
# 00_functions.R
# Phenotype- and genomics-based data compression for efficient GBLUP
# Core library (R only, open-source packages: lme4, Matrix, irlba)
#
# Every variance-component model in this project is fitted by REML with lme4
# (Bates et al., 2015). Kernels are handled through a spectral rotation
# (Kang et al., 2008; Lippert et al., 2011): after rotating the data with the
# eigenvectors of the (weighted) genomic kernel, the model becomes a mixed
# model with a diagonal random-effect design matrix, which is passed to lme4
# by replacing the Zt slot of lFormula() (the same device used by lme4qtl,
# Ziyatdinov et al., 2018). The REML estimates are therefore exactly those of
# the weighted GBLUP model, and are reproducible with any recent R version.
###############################################################################
suppressPackageStartupMessages({
  library(lme4)
  library(Matrix)
})

## ---------------------------------------------------------------------------
## 1. Markers
## ---------------------------------------------------------------------------
std_markers <- function(X) {
  # Center and scale every marker column (genotypes only, no phenotype used)
  m <- colMeans(X)
  s <- sqrt(colSums(sweep(X, 2, m)^2) / (nrow(X) - 1))
  s[s < 1e-8] <- 1
  sweep(sweep(X, 2, m, "-"), 2, s, "/")
}

top_pcs <- function(Z, k = 30, seed = 1) {
  # First k principal-component scores of Z (randomized/Lanczos SVD, irlba)
  set.seed(seed)
  s <- irlba::irlba(Z, nv = k, nu = k)
  s$u %*% diag(s$d, k)
}

kmeans_lloyd <- function(P, K, seed = 1, iter = 25) {
  # Lloyd's k-means (Lloyd, 1982) on the PC scores. Empty clusters are dropped,
  # so the number of groups returned can be slightly smaller than K.
  set.seed(seed)
  n <- nrow(P)
  C <- P[sample.int(n, K), , drop = FALSE]
  lab <- integer(n)
  for (it in seq_len(iter)) {
    D <- -2 * tcrossprod(P, C) + rep(rowSums(C^2), each = n)
    new <- max.col(-D, ties.method = "first")
    if (it > 1 && all(new == lab)) break
    lab <- new
    C <- rowsum(P, lab) / as.vector(table(lab))
    lab <- as.integer(factor(lab))           # consecutive labels 1..K'
  }
  lab
}

## ---------------------------------------------------------------------------
## 2. Grouping rules: return an integer label (group) for every individual
##    RND random | GEN genomic clusters | YQ phenotype-quantile bins |
##    SQC genomic clusters x phenotype-quantile bins
## ---------------------------------------------------------------------------
make_groups <- function(method, y, Q, PC = NULL, B = 4, seed = 1) {
  n <- length(y)
  set.seed(seed)
  if (method == "RND") return(sample(rep_len(seq_len(Q), n)))
  if (method == "YQ") {
    lab <- integer(n)
    lab[order(y)] <- ceiling(seq_len(n) * Q / n)  # Q equal-count bins of sorted y
    return(lab)
  }
  if (method == "GEN") return(kmeans_lloyd(PC, Q, seed))
  if (method == "SQC") {
    C <- max(Q %/% B, 2)
    cl <- kmeans_lloyd(PC, C, seed)
    lab <- integer(n); nxt <- 0L
    for (c in sort(unique(cl))) {
      ix <- which(cl == c); ix <- ix[order(y[ix])]
      nb <- min(B, length(ix))
      lab[ix] <- nxt + ceiling(seq_along(ix) * nb / length(ix))
      nxt <- nxt + nb
    }
    return(lab)
  }
  stop("unknown method")
}

summarise_groups <- function(Z, y, lab) {
  # Stage 1: group means (BLUEs of the group effects) of markers and phenotype
  lab <- as.integer(factor(lab))
  m <- as.vector(table(lab))
  list(Zbar = rowsum(Z, lab, reorder = TRUE) / m,
       ybar = as.vector(rowsum(y, lab, reorder = TRUE)) / m,
       m = m)
}

## ---------------------------------------------------------------------------
## 3. REML with lme4 on a rotated (diagonal) model
##    ys = X beta + diag(sqrt(d)) u + e,  u ~ N(0, s2g I), e ~ N(0, s2e I)
## ---------------------------------------------------------------------------
lme4_ctrl <- lmerControl(check.nobs.vs.nlev = "ignore",
                         check.nobs.vs.nRE  = "ignore",
                         check.nlev.gtreq.5 = "ignore",
                         check.rankX = "ignore",
                         calc.derivs = FALSE)

reml_rotated_lme4 <- function(ys, d, xs = NULL) {
  nobs <- length(ys); nre <- length(d)
  dat <- data.frame(ys = ys, id = factor(c(seq_len(nre), rep(1L, nobs - nre)),
                                         levels = seq_len(nre)))
  if (is.null(xs)) {
    lf <- lFormula(ys ~ 0 + (1 | id), data = dat, REML = TRUE, control = lme4_ctrl)
  } else {
    dat$xs <- xs
    lf <- lFormula(ys ~ 0 + xs + (1 | id), data = dat, REML = TRUE, control = lme4_ctrl)
  }
  # Replace the indicator design by the diagonal sqrt(eigenvalue) design
  lf$reTrms$Zt <- sparseMatrix(i = seq_len(nre), j = seq_len(nre),
                               x = pmax(sqrt(pmax(d, 0)), 1e-10), dims = c(nre, nobs))
  devfun <- do.call(mkLmerDevfun, lf)
  opt <- optimizeLmer(devfun)
  fit <- mkMerMod(environment(devfun), opt, lf$reTrms, fr = lf$fr)
  theta <- as.numeric(getME(fit, "theta"))
  s2e <- sigma(fit)^2
  list(fit = fit, s2e = s2e, s2g = s2e * theta^2,
       delta = if (theta > 1e-5) 1 / theta^2 else 1e10,
       mu = if (is.null(xs)) 0 else unname(fixef(fit)[1]))
}

## ---------------------------------------------------------------------------
## 4. Weighted GBLUP on (group) means = Stage 2 of the two-stage analysis
##    ybar = mu 1 + Zbar b + ebar,  b ~ N(0, s2b I),  ebar ~ N(0, s2e D^-1)
##    D = diag(m). With m = 1 and Zbar = Z this is the full (one-stage) GBLUP.
## ---------------------------------------------------------------------------
fit_wgblup <- function(Zbar, ybar, m) {
  p <- ncol(Zbar); s <- sqrt(m)
  K <- tcrossprod(Zbar) / p
  if (any(m != 1)) K <- K * tcrossprod(s)          # K* = D^1/2 K D^1/2
  e <- eigen(K, symmetric = TRUE); rm(K)
  d <- pmax(e$values, 0); U <- e$vectors
  ys <- drop(crossprod(U, ybar * s)); xs <- drop(crossprod(U, s))
  vc <- reml_rotated_lme4(ys, d, xs)
  structure(list(U = U, d = d, ys = ys, xs = xs, s = s, Zbar = Zbar, p = p,
                 delta = vc$delta, s2g = vc$s2g, s2e = vc$s2e, mu_reml = vc$mu,
                 Q = length(ybar)), class = "wgblup")
}

coef_wgblup <- function(f, delta = f$delta) {
  w <- 1 / (f$d + delta)
  mu <- sum(w * f$xs * f$ys) / sum(w * f$xs^2)          # GLS intercept
  alpha <- f$U %*% (w * (f$ys - mu * f$xs))
  b <- drop(crossprod(f$Zbar, f$s * alpha)) / f$p       # marker effects
  list(mu = mu, b = b)
}

predict_wgblup <- function(f, Znew, delta = f$delta) {
  cf <- coef_wgblup(f, delta)
  drop(cf$mu + Znew %*% cf$b)
}

## ---------------------------------------------------------------------------
## 5. Nystrom low-rank GBLUP (Williams & Seeger, 2001), REML with lme4
##    K ~ K_nQ K_QQ^+ K_Qn = Phi Phi',  Phi = K_nQ K_QQ^-1/2  (n x Q)
## ---------------------------------------------------------------------------
fit_nystrom <- function(Ztr, ytr, Q, seed = 1) {
  n <- nrow(Ztr); p <- ncol(Ztr)
  set.seed(seed); idx <- sample.int(n, Q)
  Zm <- Ztr[idx, , drop = FALSE]
  em <- eigen(tcrossprod(Zm) / p, symmetric = TRUE)
  keep <- em$values > 1e-8 * max(em$values)
  Tm <- em$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(em$values[keep]), sum(keep))
  Phi <- (tcrossprod(Ztr, Zm) / p) %*% Tm               # n x r feature matrix
  cap_train <- mean(rowSums(Phi^2)) / mean(rowSums(Ztr^2) / p)
  mphi <- colMeans(Phi); Phi <- sweep(Phi, 2, mphi)
  ym <- mean(ytr); yc <- ytr - ym
  eg <- eigen(crossprod(Phi), symmetric = TRUE)
  ks <- seq_len(min(length(eg$values), n - 2))        # at most n - 2 non-null directions
  S2 <- pmax(eg$values[ks], 1e-12); V <- eg$vectors[, ks, drop = FALSE]
  Vq <- drop(crossprod(V, crossprod(Phi, yc))); rm(Phi)
  ysr <- Vq / sqrt(S2); r <- length(S2)
  rss <- max(sum(yc^2) - sum(ysr^2), 1e-12)
  # n-1 error contrasts: r rotated observations + (n-1-r) pure-residual
  # contrasts, whose likelihood depends only on their sum of squares (rss)
  ys <- c(ysr, sqrt(rss), rep(0, n - 2 - r))
  vc <- reml_rotated_lme4(ys, S2, xs = NULL)
  beta <- drop(V %*% (Vq / (S2 + vc$delta)))
  structure(list(Zm = Zm, Tm = Tm, mphi = mphi, ym = ym, beta = beta, p = p,
                 delta = vc$delta, cap_train = cap_train), class = "nys")
}

predict_nystrom <- function(f, Znew) {
  Phin <- (tcrossprod(Znew, f$Zm) / f$p) %*% f$Tm
  list(pred = drop(f$ym + sweep(Phin, 2, f$mphi) %*% f$beta),
       cap = mean(rowSums(Phin^2)) / mean(rowSums(Znew^2) / f$p))
}

## ---------------------------------------------------------------------------
## 6. One replicate: every method at every Q
## ---------------------------------------------------------------------------
tnow <- function() proc.time()[["elapsed"]]
cor2 <- function(a, b) suppressWarnings(cor(a, b))

metrics <- function(pred, yte, gte) {
  c(r = cor2(pred, yte),
    r_g = if (is.null(gte)) NA else cor2(pred, gte),
    slope = cov(pred, yte) / var(pred),
    mse = mean((pred - yte)^2))
}

evaluate_rep <- function(Ztr, ytr, Zte, yte, gte = NULL, Qs, seed,
                         val_frac = 0.1, do_full = TRUE, logf = message) {
  rows <- list(); n <- length(ytr)
  rec <- function(method, Q, shrink, pred, t, groups, delta = NA, cap_tr = NA, cap_te = NA)
    rows[[length(rows) + 1]] <<- data.frame(method = method, Q = Q, shrink = shrink,
      t(metrics(pred, yte, gte)), time = t, groups = groups, delta = delta,
      cap_train = cap_tr, cap_test = cap_te)
  if (do_full) {
    t0 <- tnow(); f <- fit_wgblup(Ztr, ytr, rep(1, n)); pr <- predict_wgblup(f, Zte)
    rec("FULL", n, "REML", pr, tnow() - t0, n, f$delta); rm(f); gc(verbose = FALSE)
    logf("   FULL done")
  }
  set.seed(seed)
  perm <- sample.int(n); nv <- floor(val_frac * n)
  iv <- perm[seq_len(nv)]; it <- perm[-seq_len(nv)]
  t0 <- tnow(); P <- top_pcs(Ztr, 30, seed); tP <- tnow() - t0
  Pt <- P[it, ]; Zt <- Ztr[it, ]; yt <- ytr[it]; Zv <- Ztr[iv, ]; yv <- ytr[iv]
  for (Q in Qs) {
    # SUB: random subsample of Q individuals (the rest discarded)
    set.seed(seed + Q)
    t0 <- tnow(); idx <- sample.int(n, Q)
    f <- fit_wgblup(Ztr[idx, ], ytr[idx], rep(1, Q))
    rec("SUB", Q, "REML", predict_wgblup(f, Zte), tnow() - t0, Q, f$delta)
    # NYS: Nystrom with Q random landmarks, all individuals
    t0 <- tnow(); ny <- fit_nystrom(Ztr, ytr, Q, seed + Q); pn <- predict_nystrom(ny, Zte)
    rec("NYS", Q, "REML", pn$pred, tnow() - t0, Q, ny$delta, ny$cap_train, pn$cap)
    rm(ny); gc(verbose = FALSE)
    for (meth in c("RND", "GEN", "YQ", "SQC")) {
      tadd <- if (meth %in% c("GEN", "SQC")) tP else 0
      t0 <- tnow()
      lab <- make_groups(meth, ytr, Q, P, seed = seed)
      S <- summarise_groups(Ztr, ytr, lab)
      f <- fit_wgblup(S$Zbar, S$ybar, S$m)
      rec(meth, Q, "REML", predict_wgblup(f, Zte), tnow() - t0 + tadd, length(S$m), f$delta)
      if (meth != "RND") {
        # 90% compressed, 10% held out: VAL (delta tuned) and CAL (recalibrated)
        t0 <- tnow()
        lab <- make_groups(meth, yt, Q, Pt, seed = seed)
        S <- summarise_groups(Zt, yt, lab)
        f <- fit_wgblup(S$Zbar, S$ybar, S$m)
        t_fit <- tnow() - t0
        grid <- exp(seq(log(f$delta) - 5, log(f$delta) + 5, length.out = 21))
        rv <- sapply(grid, function(dl) cor2(predict_wgblup(f, Zv, dl), yv))
        best <- grid[which.max(rv)]
        rec(meth, Q, "VAL", predict_wgblup(f, Zte, best), tnow() - t0 + tadd, length(S$m), best)
        t0 <- tnow()
        pv <- predict_wgblup(f, Zv); pt <- predict_wgblup(f, Zte)
        bb <- cov(pv, yv) / var(pv); aa <- mean(yv) - bb * mean(pv)
        rec(meth, Q, "CAL", aa + bb * pt, t_fit + tnow() - t0 + tadd, length(S$m), f$delta)
      }
      rm(f, S); gc(verbose = FALSE)
    }
    logf(sprintf("   Q = %d done", Q))
  }
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------------
## 7. Simulation of genotypes and phenotypes
## ---------------------------------------------------------------------------
simulate_data <- function(scenario, n, p, rho = 0.95, h2 = 0.5) {
  maf <- runif(p, 0.05, 0.50)
  nsub <- if (scenario == "S3_structure") 5 else 1
  sub <- (sample.int(n) - 1) %% nsub + 1
  if (nsub > 1) {   # Balding-Nichols subpopulation frequencies, F = 0.10
    Fst <- 0.10
    fk <- matrix(rbeta(nsub * p, rep(maf * (1 - Fst) / Fst, each = nsub),
                       rep((1 - maf) * (1 - Fst) / Fst, each = nsub)), nsub, p)
    thr <- qnorm(1 - pmin(pmax(fk, 0.01), 0.99))
    thr <- matrix(thr, nsub, p)
  } else thr <- matrix(qnorm(1 - maf), 1, p)
  X <- matrix(0, n, p)
  for (h in 1:2) {               # two haplotypes from an AR(1) latent process
    z <- matrix(0, n, p); z[, 1] <- rnorm(n)
    E <- matrix(rnorm(n * p, sd = sqrt(1 - rho^2)), n, p)
    for (j in 2:p) z[, j] <- rho * z[, j - 1] + E[, j]
    X <- X + (z > thr[sub, , drop = FALSE])
    rm(z, E)
  }
  nq <- if (scenario == "S2_sparse") 30 else 500
  q <- sample.int(p, nq)
  g <- drop((X[, q] - 1) %*% rnorm(nq))
  g <- (g - mean(g)) / sd(g) * sqrt(h2)
  y <- g + rnorm(n, sd = sqrt(1 - h2))
  list(X = X, y = y, g = g)
}

SIM <- list(N = 8000, N_TRAIN = 6400, P = 3000, QS = c(250, 1000, 2000, 3000),
            SCEN = c("S1_polygenic", "S2_sparse", "S3_structure"), NREP = 5, SEED = 2025)

make_sim_rep <- function(si, rep) {
  # Reproducible data and partition of replicate `rep` of scenario `si`
  set.seed(SIM$SEED * 1000 + si * 100 + rep)
  d <- simulate_data(SIM$SCEN[si], SIM$N, SIM$P)
  Z <- std_markers(d$X)
  perm <- sample.int(SIM$N)
  list(Z = Z, y = d$y, g = d$g, tr = perm[seq_len(SIM$N_TRAIN)],
       te = perm[-seq_len(SIM$N_TRAIN)], seed = sample.int(1e8, 1))
}

## ===========================================================================
## 8. Methods added in the revision
## ---------------------------------------------------------------------------
## 8a. Min_GRM training-set selection (Fernandez-Gonzalez et al., 2024)
##     Keeps Q individual records (no averaging). For every candidate training
##     individual i, the score is its minimum genomic relationship with the
##     testing set (TS):  s_i = min_{j in TS} G_ij,  G_TRS,TS = Z_TRS Z_TS'/p.
##     The Q candidates with the largest scores are kept, so the selected
##     training set (TRS) maximizes the minimum relationship between each of its
##     individuals and any individual of the TS. Uses the genotypes of the
##     training and testing individuals only (no phenotypes), so it is a
##     targeted (TS-aware) selection. Ties are broken by the row order.
##     Cost: O(n n_TS p) for the cross-relationships, computed in blocks of the
##     training individuals so that only a block x n_TS matrix is held in memory.
## ---------------------------------------------------------------------------
mingrm_score <- function(Ztr, Zte, block = 2000) {
  p <- ncol(Ztr); s <- numeric(nrow(Ztr))
  for (st in seq(1, nrow(Ztr), by = block)) {
    ix <- st:min(nrow(Ztr), st + block - 1)
    Gb <- tcrossprod(Ztr[ix, , drop = FALSE], Zte) / p      # block x n_TS
    s[ix] <- apply(Gb, 1, min)
  }
  s
}

select_mingrm <- function(Ztr, Zte, Q, score = mingrm_score(Ztr, Zte)) {
  order(-score, seq_along(score))[seq_len(Q)]
}

## ---------------------------------------------------------------------------
## 8b. Exact GBLUP in marker space (p x p), useful when p < n
##     One pass over the data in blocks accumulates Z'Z, Z'y, 1'Z, y'y and n.
##     With intercept, REML uses the centered cross-products
##     Zc'Zc = Z'Z - n zbar zbar',  Zc'yc = Z'y - n zbar ybar.
##     The non-null eigenvalues of Zc Zc'/p are those of Zc'Zc/p = V S V'/p;
##     the rotated observations are S^-1/2 V' Zc'yc, and the n - 1 - r
##     null-space contrasts enter the REML likelihood only through their sum
##     of squares, yc'yc - ||S^-1/2 V' Zc'yc||^2 (same device as fit_nystrom).
##     Time O(n p^2), memory O(p^2), independent of n. The REML estimate and
##     the predictions are identical to those of the n x n kernel GBLUP.
## ---------------------------------------------------------------------------
crossprod_blocks <- function(Z, y, block = 2000) {
  p <- ncol(Z); ZtZ <- matrix(0, p, p); Zty <- numeric(p); s1 <- numeric(p)
  syy <- 0; sy <- 0; n <- 0
  for (st in seq(1, nrow(Z), by = block)) {
    ix <- st:min(nrow(Z), st + block - 1); Zb <- Z[ix, , drop = FALSE]; yb <- y[ix]
    ZtZ <- ZtZ + crossprod(Zb); Zty <- Zty + drop(crossprod(Zb, yb))
    s1 <- s1 + colSums(Zb); syy <- syy + sum(yb^2); sy <- sy + sum(yb); n <- n + length(ix)
  }
  list(ZtZ = ZtZ, Zty = Zty, s1 = s1, syy = syy, sy = sy, n = n, p = p)
}

fit_marker_gblup <- function(Z, y, block = 2000) {
  cp <- crossprod_blocks(Z, y, block); n <- cp$n; p <- cp$p
  zbar <- cp$s1 / n; ybar <- cp$sy / n
  ZcZc <- cp$ZtZ - n * tcrossprod(zbar); Zcy <- cp$Zty - n * zbar * ybar
  yy <- cp$syy - n * ybar^2; rm(cp)
  e <- eigen(ZcZc, symmetric = TRUE); rm(ZcZc)
  keep <- e$values > 1e-9 * max(e$values)
  keep[seq_along(keep) > n - 2] <- FALSE                  # at most n - 2 non-null directions
  S <- e$values[keep]; V <- e$vectors[, keep, drop = FALSE]; rm(e)
  Vz <- drop(crossprod(V, Zcy)); ysr <- Vz / sqrt(S); r <- length(S)
  rss <- max(yy - sum(ysr^2), 1e-12)
  ys <- c(ysr, sqrt(rss), rep(0, n - 2 - r))
  vc <- reml_rotated_lme4(ys, S / p, xs = NULL)
  b <- drop(V %*% (Vz / (S + p * vc$delta)))             # (Zc'Zc + p delta I)^-1 Zc'yc
  structure(list(b = b, mu = ybar - sum(zbar * b), delta = vc$delta,
                 s2g = vc$s2g, s2e = vc$s2e, rank = r), class = "mgblup")
}
predict_marker_gblup <- function(f, Znew) drop(f$mu + Znew %*% f$b)

## ---------------------------------------------------------------------------
## 8c. One-pass marginal-regression predictor (limit of ridge regression when
##     delta -> infinity): b proportional to Z'(y - ybar). Correlations do not
##     depend on the scale of b; intercept and scale are set by regressing y on
##     the training scores, so that MSE and slope are also meaningful.
## ---------------------------------------------------------------------------
fit_marginal <- function(Z, y) {
  b <- drop(crossprod(Z, y - mean(y))) / length(y)
  s <- drop(Z %*% b); cf <- coef(lm(y ~ s))
  list(b = b * cf[[2]], mu = cf[[1]])
}
predict_marginal <- function(f, Znew) drop(f$mu + Znew %*% f$b)

## ---------------------------------------------------------------------------
## 8d. One replicate of the added methods (same data, partitions and Q values
##     as evaluate_rep; deterministic, no random numbers involved).
##     `methods` selects which of EXACT, MARG and MINGRM are computed.
## ---------------------------------------------------------------------------
evaluate_added <- function(Ztr, ytr, Zte, yte, gte = NULL, Qs, logf = message,
                           methods = c("EXACT", "MARG", "MINGRM")) {
  rows <- list()
  rec <- function(method, Q, pred, t, groups, delta = NA)
    rows[[length(rows) + 1]] <<- data.frame(method = method, Q = Q, shrink = "REML",
      t(metrics(pred, yte, gte)), time = t, groups = groups, delta = delta,
      cap_train = NA, cap_test = NA)
  n <- length(ytr)
  if ("EXACT" %in% methods) {
    t0 <- tnow(); fm <- fit_marker_gblup(Ztr, ytr); pm <- predict_marker_gblup(fm, Zte)
    rec("EXACT", n, pm, tnow() - t0, n, fm$delta); rm(fm)
  }
  if ("MARG" %in% methods) {
    t0 <- tnow(); fg <- fit_marginal(Ztr, ytr)
    rec("MARG", n, predict_marginal(fg, Zte), tnow() - t0, n)
  }
  if ("MINGRM" %in% methods) {
    t0 <- tnow(); sc <- mingrm_score(Ztr, Zte); t_score <- tnow() - t0  # computed once for all Q
    for (Q in Qs) {
      t0 <- tnow(); idx <- select_mingrm(Ztr, Zte, Q, sc)
      f <- fit_wgblup(Ztr[idx, , drop = FALSE], ytr[idx], rep(1, Q))
      rec("MINGRM", Q, predict_wgblup(f, Zte), tnow() - t0 + t_score, Q, f$delta)
      rm(f); gc(verbose = FALSE)
    }
  }
  logf("   added methods done")
  if (length(rows)) do.call(rbind, rows) else NULL
}

## ===========================================================================
## 9. Fully efficient (exact) phenotype-quantile compression  -- YQ-FE (v8)
## ---------------------------------------------------------------------------
## For any grouping (in particular the phenotype-quantile bins of YQ), the
## block-diagonal orthogonal matrix H = diag(H_q), H_q = [1/sqrt(m_q) 1, C_q]
## (C_q: m_q x (m_q - 1) orthonormal within-group contrasts) splits the data of
## each group into its scaled mean sqrt(m_q) ybar_q and m_q - 1 within-group
## contrasts C_q'y_q. Because H is orthogonal,
##   H'y = [D^1/2 ybar ; C'y],   H'Z = [D^1/2 Zbar ; C'Z],   H'1 = [D^1/2 1 ; 0]
## and the residuals stay iid N(0, s2e I): the likelihood and the REML criterion
## of the transformed data are those of the one-stage model, whatever the groups.
## Stage 1 (per group): ybar_q, zbar_q, m_q and the within-group sums
##   W = sum_q (Z_q - 1 zbar_q')'(Z_q - 1 zbar_q')   (p x p, within-group SSCP)
##   sW = sum_q (Z_q - 1 zbar_q')'(y_q - ybar_q 1)   (p)
##   tW = sum_q ||y_q - ybar_q 1||^2,                nW = n - Q
## Stage 2: the weighted between-group model of YQ PLUS the within-group
## information; dropping W, sW, tW gives back the usual (approximate) YQ.
## ---------------------------------------------------------------------------
yq_fe_stage1 <- function(Z, y, lab, block = 2000) {
  lab <- as.integer(factor(lab)); m <- as.vector(table(lab)); p <- ncol(Z)
  Zbar <- rowsum(Z, lab, reorder = TRUE) / m
  ybar <- as.vector(rowsum(y, lab, reorder = TRUE)) / m
  W <- matrix(0, p, p); sW <- numeric(p); tW <- 0
  for (st in seq(1, nrow(Z), by = block)) {        # one pass, in blocks of individuals
    ix <- st:min(nrow(Z), st + block - 1)
    Zw <- Z[ix, , drop = FALSE] - Zbar[lab[ix], , drop = FALSE]
    yw <- y[ix] - ybar[lab[ix]]
    W <- W + crossprod(Zw); sW <- sW + drop(crossprod(Zw, yw)); tW <- tW + sum(yw^2)
  }
  list(Zbar = Zbar, ybar = ybar, m = m, W = W, sW = sW, tW = tW, n = length(y), Q = length(m), p = p)
}

fit_from_crossprods <- function(ZtZ, Zty, s1, sy, syy, n, p) {
  # REML (lme4) and BLUP of the marker-space GBLUP from its sufficient statistics
  zbar <- s1 / n; ybar <- sy / n
  ZcZc <- ZtZ - n * tcrossprod(zbar); Zcy <- Zty - n * zbar * ybar; yy <- syy - n * ybar^2
  e <- eigen(ZcZc, symmetric = TRUE)
  keep <- e$values > 1e-9 * max(e$values); keep[seq_along(keep) > n - 2] <- FALSE
  S <- e$values[keep]; V <- e$vectors[, keep, drop = FALSE]
  Vz <- drop(crossprod(V, Zcy)); ysr <- Vz / sqrt(S); r <- length(S)
  rss <- max(yy - sum(ysr^2), 1e-12)
  vc <- reml_rotated_lme4(c(ysr, sqrt(rss), rep(0, n - 2 - r)), S / p, xs = NULL)
  b <- drop(V %*% (Vz / (S + p * vc$delta)))
  structure(list(b = b, mu = ybar - sum(zbar * b), delta = vc$delta, s2g = vc$s2g, s2e = vc$s2e, rank = r),
            class = "mgblup")
}

## (a) marker form: Stage 2 = mixed-model equations built from the Stage 1 summaries
fit_yq_fe <- function(S1) {
  m <- S1$m
  ZtZ <- crossprod(S1$Zbar * sqrt(m)) + S1$W                 # Zbar'D Zbar + W  = Z'Z
  Zty <- drop(crossprod(S1$Zbar, m * S1$ybar)) + S1$sW        # Zbar'D ybar + sW = Z'y
  s1 <- drop(crossprod(S1$Zbar, m)); sy <- sum(m * S1$ybar)
  syy <- sum(m * S1$ybar^2) + S1$tW
  fit_from_crossprods(ZtZ, Zty, s1, sy, syy, S1$n, S1$p)
}

## generic GBLUP on transformed data rows X (design), y (data), x (intercept column),
## iid residuals: REML through the spectral rotation of XX'/p
fit_rows_gblup <- function(X, y, x, rss_extra = NULL, n_extra = 0) {
  # rss_extra / n_extra: n_extra additional pure-residual observations (zero design,
  # zero intercept) represented by their sum of squares, as in fit_nystrom
  p <- ncol(X); e <- eigen(tcrossprod(X) / p, symmetric = TRUE)
  d <- pmax(e$values, 0); U <- e$vectors
  ys <- drop(crossprod(U, y)); xs <- drop(crossprod(U, x))
  if (n_extra > 0) {
    ys2 <- c(ys, sqrt(rss_extra), rep(0, n_extra - 1)); d2 <- c(d, rep(0, n_extra)); xs2 <- c(xs, rep(0, n_extra))
  } else { ys2 <- ys; d2 <- d; xs2 <- xs }
  vc <- reml_rotated_lme4(ys2, d2, xs2)
  w <- 1 / (d + vc$delta); mu <- sum(w * xs * ys) / sum(w * xs^2)
  alpha <- U %*% (w * (ys - mu * xs)); b <- drop(crossprod(X, alpha)) / p
  list(b = b, mu = mu, delta = vc$delta, s2g = vc$s2g, s2e = vc$s2e)
}

## orthonormal within-group contrasts (Helmert), returned as a sparse n x (n - Q) matrix
within_contrasts <- function(lab) {
  lab <- as.integer(factor(lab)); n <- length(lab); cols <- list(); k <- 0
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (q in sort(unique(lab))) {
    ix <- which(lab == q); mq <- length(ix); if (mq < 2) next
    Hm <- stats::contr.helmert(mq); Hm <- sweep(Hm, 2, sqrt(colSums(Hm^2)), "/")
    for (j in seq_len(mq - 1)) { nz <- which(Hm[, j] != 0)
      ii <- c(ii, ix[nz]); jj <- c(jj, rep(k + j, length(nz))); xx <- c(xx, Hm[nz, j]) }
    k <- k + mq - 1
  }
  Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, k))
}

## (b) kernel (two-stage, data-level) form: Stage 2 fitted to the Q scaled group
## means and the n - Q within-group contrasts jointly, with their full genomic
## covariance matrix (between groups, within groups and between means and contrasts)
fit_yq_fe_kernel <- function(Z, y, lab) {
  lab <- as.integer(factor(lab)); m <- as.vector(table(lab)); Q <- length(m)
  A <- Matrix::sparseMatrix(i = lab, j = seq_along(lab), x = 1 / m[lab], dims = c(Q, length(lab)))
  C <- within_contrasts(lab)
  X <- rbind(as.matrix(sqrt(m) * (A %*% Z)), as.matrix(Matrix::crossprod(C, Z)))
  yy <- c(sqrt(m) * as.vector(A %*% y), as.vector(Matrix::crossprod(C, y)))
  xx <- c(sqrt(m), rep(0, ncol(C)))
  fit_rows_gblup(X, yy, xx)
}

## (c) compact form: Q scaled group means + the r_W non-null within-group
## components (rotated with the eigenvectors of W) + one pseudo-observation for
## the within-group residual sum of squares (and zeros to keep n observations).
## With k < r_W leading components and no residual pseudo-observation it gives
## the reduced (approximate) versions YQ-FE(k); k = 0 is the original YQ.
fit_yq_fe_compact <- function(S1, k = NULL, exact = TRUE) {
  m <- S1$m; p <- S1$p; ew <- eigen(S1$W, symmetric = TRUE)
  rW <- sum(ew$values > 1e-9 * max(ew$values)); rW <- min(rW, S1$n - S1$Q - 1)
  if (is.null(k)) k <- rW; k <- min(k, rW)
  Xb <- sqrt(m) * S1$Zbar; yb <- sqrt(m) * S1$ybar; xb <- sqrt(m)
  if (k > 0) {
    V <- ew$vectors[, seq_len(k), drop = FALSE]; s <- sqrt(ew$values[seq_len(k)])
    Xw <- t(V) * s                                        # k x p  (Sigma_k V_k')
    yw <- drop(crossprod(V, S1$sW)) / s                   # U_k' y^w
  } else { Xw <- matrix(0, 0, p); yw <- numeric(0) }
  X <- rbind(Xb, Xw); yy <- c(yb, yw); xx <- c(xb, rep(0, k))
  if (exact && k == rW) {                                 # remaining within contrasts: pure residual
    nres <- S1$n - S1$Q - rW; rss <- max(S1$tW - sum(yw^2), 1e-12)
    f <- fit_rows_gblup(X, yy, xx, rss_extra = rss, n_extra = nres)
  } else { nres <- 0; f <- fit_rows_gblup(X, yy, xx) }
  f$k <- k; f$rW <- rW; f$nobs <- length(yy) + nres; f
}
predict_lin <- function(f, Znew) drop(f$mu + Znew %*% f$b)

## ===========================================================================
## 10. YQ-FEN: Nystrom-accelerated fully efficient phenotype-quantile method (v10)
## ---------------------------------------------------------------------------
## Input: the YQ-FE Stage 1 summaries (yq_fe_stage1). From them the total centred
## cross-products  T = Zc'Zc = Zbar'D Zbar + W - n zbar zbar',  s = Zc'yc  and
## y'y are formed exactly as in YQ-FE. YQ-FE then needs the eigen-decomposition
## of the p x p matrix T (O(p^3)) to run REML. YQ-FEN avoids it:
##  (i)  variance components by the method of moments (Haseman-Elston regression)
##       computed from the same summaries, after projecting out the h leading
##       Nystrom directions (leading genomic PCs) of step (ii); without projection,
##       with K = ZcZc'/p on the n - 1 contrasts,
##         y'y   = s2g tr(K)   + s2e (n - 1),    tr(K)   = tr(T)/p
##         y'Ky  = s2g tr(K^2) + s2e tr(K),      tr(K^2) = ||T||_F^2 / p^2,  y'Ky = s's/p
##       (O(p^2), no decomposition);
##  (ii) a rank-k Nystrom approximation of T from k landmark markers (columns J):
##         T~ = T[,J] T[J,J]^+ T[J,] = U Lambda U'          (O(p k^2 + k^3))
##  (iii) the YQ-FE mixed-model equations (T + p delta I) b = s are solved exactly
##       by conjugate gradients preconditioned with the Nystrom approximation
##       (Frangella et al. 2023):
##         M^-1 = (lambda_k + a) U (Lambda + a I)^-1 U' + (I - U U'),  a = p delta
##       each iteration costs one product T v (O(p^2)).
## For a given delta the solution is the YQ-FE / one-stage BLUP (to the CG
## tolerance); the only approximation is the moment estimate of delta.
## Predictions do not depend on Q (the summaries T, s are those of all records);
## k changes the number of CG iterations and, through the h projected directions,
## slightly the moment estimate of delta.
## ---------------------------------------------------------------------------
yq_fe_totals <- function(S1) {
  m <- S1$m; n <- S1$n
  ZtZ <- crossprod(S1$Zbar * sqrt(m)) + S1$W                  # Zbar'D Zbar + W = Z'Z
  Zty <- drop(crossprod(S1$Zbar, m * S1$ybar)) + S1$sW         # Zbar'D ybar + sW = Z'y
  zbar <- drop(crossprod(S1$Zbar, m)) / n; ybar <- sum(m * S1$ybar) / n
  syy <- sum(m * S1$ybar^2) + S1$tW
  list(Tc = ZtZ - n * tcrossprod(zbar), sc = Zty - n * zbar * ybar, yy = syy - n * ybar^2,
       zbar = zbar, ybar = ybar, n = n, p = S1$p)
}
moments_delta <- function(Tt, U = NULL) {
  ## moment (Haseman-Elston) equations on the n - 1 centred contrasts; if U (p x h,
  ## leading Nystrom directions of T) is given, the h leading genomic principal
  ## components B = Zc U are first projected out of the data (as fixed covariates),
  ## which removes the high-variance leading directions (population structure):
  ##   P = I - B (B'B)^-1 B',  B'B = U'TU,  Zc'P Zc = T - TU (U'TU)^-1 U'T = R
  p <- Tt$p; N <- Tt$n - 1
  if (is.null(U) || ncol(U) == 0) { R <- Tt$Tc; r <- Tt$sc; yy <- Tt$yy; h <- 0 } else {
    C1 <- Tt$Tc %*% U; Gi <- solve(crossprod(U, C1)); us <- drop(crossprod(U, Tt$sc))
    R <- Tt$Tc - C1 %*% Gi %*% t(C1); r <- Tt$sc - drop(C1 %*% (Gi %*% us))
    yy <- Tt$yy - sum(us * drop(Gi %*% us)); h <- ncol(U) }
  trK <- sum(diag(R)) / p; trK2 <- sum(R^2) / p^2; yKy <- sum(r^2) / p
  A <- matrix(c(trK, trK2, N - h, trK), 2)
  ## identifiability: det(A) = trK^2 - (N - h) trK2 < 0 unless all eigenvalues of
  ## the projected kernel are equal (Cauchy-Schwarz); relative det reported.
  ## Negative solutions are truncated (s2g <= 0 -> delta = 1e8, no genomic signal;
  ## s2e <= 0 -> delta = 1e-8) and flagged.
  vc <- solve(A, c(yy, yKy)); s2g <- max(vc[1], 1e-8); s2e <- max(vc[2], 1e-8)
  list(s2g = s2g, s2e = s2e, delta = min(max(s2e / s2g, 1e-8), 1e8), h = h,
       truncated = any(vc <= 0), rel_det = (trK^2 - (N - h) * trK2) / ((N - h) * trK2))
}
nystrom_markers <- function(Tc, k, seed = 1) {
  p <- ncol(Tc); k <- min(k, p); set.seed(seed); J <- sort(sample.int(p, k))
  ej <- eigen(Tc[J, J], symmetric = TRUE); kp <- ej$values > 1e-9 * ej$values[1]
  Phi <- Tc[, J] %*% (ej$vectors[, kp, drop = FALSE] %*% (t(ej$vectors[, kp, drop = FALSE]) / sqrt(ej$values[kp])))
  ef <- eigen(crossprod(Phi), symmetric = TRUE); kk <- ef$values > 1e-9 * ef$values[1]
  lam <- ef$values[kk]
  U <- Phi %*% sweep(ef$vectors[, kk, drop = FALSE], 2, sqrt(lam), "/")
  list(U = U, lam = lam, J = J, k = k)
}
pcg_nystrom <- function(Tc, rhs, a, L = NULL, x0 = NULL, tol = 1e-10, maxit = 2000) {
  if (is.null(L)) Minv <- function(v) v else {
    U <- L$U; lam <- L$lam; mu <- lam[length(lam)] + a; f <- mu / (lam + a) - 1
    Minv <- function(v) drop(U %*% (f * drop(crossprod(U, v)))) + v }
  x <- if (is.null(x0)) numeric(length(rhs)) else x0
  r <- rhs - (drop(Tc %*% x) + a * x); z <- Minv(r); d <- z; rz <- sum(r * z)
  nb <- sqrt(sum(rhs^2)); it <- 0
  while (sqrt(sum(r^2)) > tol * nb && it < maxit) {
    Ad <- drop(Tc %*% d) + a * d; al <- rz / sum(d * Ad)
    x <- x + al * d; r <- r - al * Ad; z <- Minv(r); rz2 <- sum(r * z)
    d <- z + (rz2 / rz) * d; rz <- rz2; it <- it + 1
  }
  list(x = x, iter = it, relres = sqrt(sum(r^2)) / nb)
}
## delta: NULL = moment estimate; a number = given delta (e.g. the REML estimate,
## used to check that the CG solution reproduces YQ-FE)
fit_yq_fen <- function(S1, k = 500, h = 200, seed = 1, delta = NULL, tol = 1e-10) {
  ## k: landmark markers of the Nystrom approximation (k = 0: no Nystrom, plain
  ## moments and plain CG); h: leading Nystrom directions projected out in the
  ## moment equations (h <= k)
  t0 <- tnow(); Tt <- yq_fe_totals(S1)
  L <- if (k > 0) nystrom_markers(Tt$Tc, k, seed) else NULL; t1 <- tnow()
  hh <- if (is.null(L)) 0 else min(h, ncol(L$U))
  mo <- moments_delta(Tt, if (hh > 0) L$U[, seq_len(hh), drop = FALSE] else NULL)
  d <- if (is.null(delta)) mo$delta else delta; t2 <- tnow()
  s <- pcg_nystrom(Tt$Tc, Tt$sc, Tt$p * d, L, tol = tol); t3 <- tnow()
  structure(list(b = s$x, mu = Tt$ybar - sum(Tt$zbar * s$x), delta = d, s2g = mo$s2g, s2e = mo$s2e,
                 iter = s$iter, relres = s$relres, k = k, h = hh, truncated = mo$truncated, rel_det = mo$rel_det,
                 time_nystrom = t1 - t0, time_moments = t2 - t1, time_cg = t3 - t2), class = "yqfen")
}

## ===========================================================================
## 11. Revision (v12): additional evaluation tools
## ---------------------------------------------------------------------------
## 11a. Simulated breeding population with full-sib families: founders from the
##      AR(1) haplotype model, Nf full-sib families of size m from random pairs
##      of founders; each offspring gamete is a recombinant of the two parental
##      haplotypes (Poisson(1.5) crossovers along the marker sequence).
simulate_families <- function(n = 8000, p = 3000, nfam = 400, nfounders = 200, rho = 0.95, h2 = 0.5, nq = 500) {
  maf <- runif(p, 0.05, 0.50); thr <- qnorm(1 - maf)
  hap <- function(nn) { z <- matrix(0, nn, p); z[, 1] <- rnorm(nn)
    E <- matrix(rnorm(nn * p, sd = sqrt(1 - rho^2)), nn, p)
    for (j in 2:p) z[, j] <- rho * z[, j - 1] + E[, j]
    (z > rep(thr, each = nn)) * 1L }
  H1 <- hap(nfounders); H2 <- hap(nfounders)
  gamete <- function(i) { nc <- rpois(1, 1.5); pts <- sort(sample.int(p - 1, nc))
    src <- rep(sample(0:1, 1), p); if (nc > 0) for (b in pts) src[(b + 1):p] <- 1 - src[(b + 1):p]
    ifelse(src == 0, H1[i, ], H2[i, ]) }
  m <- n / nfam; X <- matrix(0L, n, p); fam <- integer(n); r <- 0
  for (f in seq_len(nfam)) { par <- sample.int(nfounders, 2)
    for (o in seq_len(m)) { r <- r + 1; X[r, ] <- gamete(par[1]) + gamete(par[2]); fam[r] <- f } }
  q <- sample.int(p, nq); g <- drop((X[, q] - 1) %*% rnorm(nq))
  g <- (g - mean(g)) / sd(g) * sqrt(h2); y <- g + rnorm(n, sd = sqrt(1 - h2))
  list(X = X, y = y, g = g, fam = fam)
}
## additional scenarios: S4 low heritability (h2 = 0.2), S5 high heritability (0.8),
## S6 full-sib families (random split), S7 full-sib families with whole families
## in the testing set (prediction of new families)
SIM2 <- list(SCEN = c("S4_h2_0.2", "S5_h2_0.8", "S6_families", "S7_new_families"), NREP = 5)
make_sim2_rep <- function(si, rep) {
  set.seed(SIM$SEED * 1000 + 500 + si * 100 + rep)
  sc <- SIM2$SCEN[si]
  if (si <= 2) { d <- simulate_data("S1_polygenic", SIM$N, SIM$P, h2 = c(0.2, 0.8)[si]); fam <- NULL }
  else { d <- simulate_families(SIM$N, SIM$P); fam <- d$fam }
  Z <- std_markers(d$X)
  if (si == 4) { tf <- sample(unique(fam), round(0.2 * length(unique(fam)))); te <- which(fam %in% tf); tr <- which(!fam %in% tf) }
  else { perm <- sample.int(SIM$N); tr <- perm[seq_len(SIM$N_TRAIN)]; te <- perm[-seq_len(SIM$N_TRAIN)] }
  list(Z = Z, y = d$y, g = d$g, tr = tr, te = te, seed = sample.int(1e8, 1), scenario = sc)
}

## 11b. direct iterative solution in marker space (no grouping, no Nystrom):
##      Z'Z, Z'y, y'y accumulated directly, plain moments, plain CG
fit_direct_cg <- function(Z, y, block = 2000, tol = 1e-10) {
  cp <- crossprod_blocks(Z, y, block); n <- cp$n; p <- cp$p
  zbar <- cp$s1 / n; ybar <- cp$sy / n
  Tt <- list(Tc = cp$ZtZ - n * tcrossprod(zbar), sc = cp$Zty - n * zbar * ybar, yy = cp$syy - n * ybar^2,
             zbar = zbar, ybar = ybar, n = n, p = p)
  mo <- moments_delta(Tt); s <- pcg_nystrom(Tt$Tc, Tt$sc, p * mo$delta, NULL, tol = tol)
  list(b = s$x, mu = ybar - sum(zbar * s$x), delta = mo$delta, iter = s$iter)
}

## 11c. agreement of a prediction vector with the one-stage predictions
agree <- function(pr, pm, top = 0.10) {
  nt <- ceiling(top * length(pm)); a <- order(-pm)[seq_len(nt)]; b <- order(-pr)[seq_len(nt)]
  c(max_abs_diff = max(abs(pr - pm)), rms_diff = sqrt(mean((pr - pm)^2)), rel_rms = sqrt(mean((pr - pm)^2)) / sd(pm),
    top10_overlap = length(intersect(a, b)) / nt, cor_onestage = cor(pr, pm))
}

## 11d. YQ with cross-fitted calibration (all training individuals used for the fit;
##      the calibration intercept and slope come from K-fold out-of-fold predictions)
fit_yq_cvcal <- function(Ztr, ytr, Q, K = 5, seed = 1) {
  set.seed(seed); fold <- sample(rep_len(seq_len(K), length(ytr))); oof <- numeric(length(ytr))
  for (k in seq_len(K)) { it <- fold != k
    S <- summarise_groups(Ztr[it, ], ytr[it], make_groups("YQ", ytr[it], Q)); f <- fit_wgblup(S$Zbar, S$ybar, S$m)
    oof[!it] <- predict_wgblup(f, Ztr[!it, , drop = FALSE]) }
  bb <- cov(oof, ytr) / var(oof); aa <- mean(ytr) - bb * mean(oof)
  S <- summarise_groups(Ztr, ytr, make_groups("YQ", ytr, Q)); f <- fit_wgblup(S$Zbar, S$ybar, S$m)
  list(f = f, a = aa, b = bb)
}

## ===========================================================================
## 12. YQ-FE-Q: YQ-FE written as a GBLUP on the Q group means (v13)
## ---------------------------------------------------------------------------
## With the within-group information absorbed into the distribution of b,
##   ytilde = ybar - Zbar m_w = 1 mu + u + ebar,   u ~ N(0, s2g K_G(delta)),
##   ebar ~ N(0, s2e D^-1) = N(0, s2g delta D^-1),
##   K_G(delta) = delta Zbar (W + p delta I)^-1 Zbar'                (marker form)
##              = G_bb - G_bw (G_ww + delta I)^-1 G_wb               (kernel form)
##   m_w = (W + p delta I)^-1 s_W ,   Zbar m_w = G_bw (G_ww + delta I)^-1 y_w.
## Mixed-model equations of dimension Q + 1; predictions and REML identical to
## the one-stage GBLUP. -2 REML = (n-1) log s2g + log|V_b| + log(1'V_b^-1 1)
##   + log|G_ww + delta I|, V_b = K_G + delta D^-1, with
##   s2g = (r'V_b^-1 r + y_w'(G_ww + delta I)^-1 y_w) / (n - 1).
## ---------------------------------------------------------------------------
yq_fe_q_prep <- function(S1) {                     # marker form: one eigen-decomposition of W
  e <- eigen(S1$W, symmetric = TRUE)
  keep <- e$values > 1e-9 * max(e$values); keep[seq_along(keep) > S1$n - S1$Q] <- FALSE
  V <- e$vectors[, keep, drop = FALSE]; om <- e$values[keep]
  list(V = V, om = om, Ab = S1$Zbar %*% V, vs = drop(crossprod(V, S1$sW)), r = sum(keep))
}
## the Q-dimensional GBLUP of YQ-FE-Q for a given delta (marker form)
yq_fe_q_parts <- function(S1, P, delta) {
  p <- S1$p; lam <- p * delta
  # (W + lam I)^-1 = V diag(1/(om+lam)) V' + (I - V V')/lam
  ZZ <- tcrossprod(S1$Zbar); ZV <- P$Ab
  K <- delta * (ZZ / lam + ZV %*% (t(ZV) * (1 / (P$om + lam) - 1 / lam)))   # K_G(delta)
  mw_proj <- P$vs / (P$om + lam)                    # V' m_w (s_W lies in the span of V)
  yt <- S1$ybar - drop(ZV %*% mw_proj)              # ytilde
  qw <- (S1$tW - sum(P$vs * mw_proj)) / delta       # y_w'(G_ww + delta I)^-1 y_w
  ldw <- sum(log(P$om / p + delta)) + (S1$n - S1$Q - P$r) * log(delta)   # log|G_ww + delta I|
  list(K = K, yt = yt, qw = qw, ldw = ldw, mw_proj = mw_proj)
}
yq_fe_q_crit <- function(S1, P, delta, parts = NULL) {
  if (is.null(parts)) parts <- yq_fe_q_parts(S1, P, delta)
  Vb <- parts$K + diag(delta / S1$m, S1$Q); L <- chol(Vb)
  Vi1 <- backsolve(L, forwardsolve(t(L), rep(1, S1$Q))); Viy <- backsolve(L, forwardsolve(t(L), parts$yt))
  mu <- sum(Viy) / sum(Vi1); alpha <- Viy - mu * Vi1   # alpha = V_b^-1 (ytilde - 1 mu)
  qb <- sum((parts$yt - mu) * alpha); n <- S1$n
  s2g <- (qb + parts$qw) / (n - 1)
  list(crit = (n - 1) * log(s2g) + 2 * sum(log(diag(L))) + log(sum(Vi1)) + parts$ldw,
       mu = mu, alpha = alpha, s2g = s2g, parts = parts)
}
fit_yq_fe_q <- function(S1, delta = NULL, P = NULL, interval = c(1e-4, 1e4), tol = 1e-8) {
  t0 <- proc.time()[["elapsed"]]
  if (is.null(P)) P <- yq_fe_q_prep(S1)
  t1 <- proc.time()[["elapsed"]]; neval <- 0
  if (is.null(delta)) {
    f <- function(ld) { neval <<- neval + 1; yq_fe_q_crit(S1, P, exp(ld))$crit }
    delta <- exp(optimize(f, log(interval), tol = tol)$minimum)
  }
  cr <- yq_fe_q_crit(S1, P, delta); lam <- S1$p * delta
  # marker effects: b = m_w + delta (W + lam I)^-1 Zbar' alpha
  za <- drop(crossprod(S1$Zbar, cr$alpha))
  vza <- drop(crossprod(P$V, za))
  b <- drop(P$V %*% cr$parts$mw_proj) + delta * (za / lam + drop(P$V %*% (vza * (1 / (P$om + lam) - 1 / lam))))
  structure(list(b = b, mu = cr$mu, delta = delta, s2g = cr$s2g, s2e = cr$s2g * delta,
                 K = cr$parts$K, ytilde = cr$parts$yt, alpha = cr$alpha, crit = cr$crit, neval = neval,
                 t_prep = t1 - t0, t_total = proc.time()[["elapsed"]] - t0), class = "yqfeq")
}
## kernel (individual-space) form: works with any training kernel G (n x n) and
## cross kernel G0 (n0 x n); needs (n - Q) x (n - Q) eigen-decomposition of G_ww
fit_yq_fe_q_kernel <- function(G, y, lab, G0 = NULL, delta = NULL, interval = c(1e-4, 1e4), tol = 1e-8) {
  t0 <- proc.time()[["elapsed"]]
  lab <- as.integer(factor(lab)); m <- as.vector(table(lab)); Q <- length(m); n <- length(y)
  A <- as.matrix(Matrix::sparseMatrix(i = lab, j = seq_len(n), x = 1 / m[lab], dims = c(Q, n)))
  C <- within_contrasts(lab)
  GA <- G %*% t(A); Gbb <- A %*% GA; Gbw <- as.matrix(Matrix::crossprod(C, GA)); Gbw <- t(Gbw)   # Q x (n-Q)
  Gww <- as.matrix(Matrix::crossprod(C, Matrix::crossprod(G, C)))                               # C'GC
  ybar <- drop(A %*% y); yw <- drop(as.matrix(Matrix::crossprod(C, y)))
  e <- eigen(Gww, symmetric = TRUE); th <- pmax(e$values, 0); U <- e$vectors
  B <- Gbw %*% U; uy <- drop(crossprod(U, yw))
  t1 <- proc.time()[["elapsed"]]
  parts <- function(delta) {
    w <- 1 / (th + delta)
    list(K = Gbb - B %*% (t(B) * w), yt = ybar - drop(B %*% (w * uy)), qw = sum(w * uy^2),
         ldw = sum(log(th + delta)))
  }
  S <- list(m = m, Q = Q, n = n)
  crit <- function(delta) yq_fe_q_crit(S, NULL, delta, parts(delta))
  neval <- 0
  if (is.null(delta)) delta <- exp(optimize(function(ld) { neval <<- neval + 1; crit(exp(ld))$crit },
                                            log(interval), tol = tol)$minimum)
  cr <- crit(delta)
  pred <- NULL
  if (!is.null(G0)) {
    G0b <- G0 %*% t(A); G0w <- as.matrix(G0 %*% C); B0 <- G0w %*% U; w <- 1 / (th + delta)
    K0 <- G0b - B0 %*% (t(B) * w)
    pred <- drop(cr$mu + B0 %*% (w * uy) + K0 %*% cr$alpha)
  }
  list(mu = cr$mu, delta = delta, s2g = cr$s2g, crit = cr$crit, K = cr$parts$K, ytilde = cr$parts$yt,
       alpha = cr$alpha, pred = pred, neval = neval, t_prep = t1 - t0, t_total = proc.time()[["elapsed"]] - t0)
}
## one-stage kernel GBLUP (REML by the same profiled criterion, n x n eigen-decomposition)
fit_kernel_gblup <- function(G, y, G0 = NULL, delta = NULL, interval = c(1e-4, 1e4), tol = 1e-8) {
  t0 <- proc.time()[["elapsed"]]; n <- length(y)
  e <- eigen(G, symmetric = TRUE); d <- pmax(e$values, 0); U <- e$vectors
  uy <- drop(crossprod(U, y)); u1 <- colSums(U)
  crit <- function(delta) {
    w <- 1 / (d + delta); mu <- sum(w * u1 * uy) / sum(w * u1^2)
    q <- sum(w * (uy - mu * u1)^2); s2g <- q / (n - 1)
    list(crit = (n - 1) * log(s2g) + sum(log(d + delta)) + log(sum(w * u1^2)), mu = mu, s2g = s2g)
  }
  neval <- 0
  if (is.null(delta)) delta <- exp(optimize(function(ld) { neval <<- neval + 1; crit(exp(ld))$crit },
                                            log(interval), tol = tol)$minimum)
  cr <- crit(delta); w <- 1 / (d + delta)
  alpha <- drop(U %*% (w * (uy - cr$mu * u1)))
  pred <- if (is.null(G0)) NULL else drop(cr$mu + G0 %*% alpha)
  list(mu = cr$mu, delta = delta, s2g = cr$s2g, crit = cr$crit, pred = pred, crit_fun = function(d) crit(d)$crit,
       neval = neval, t_total = proc.time()[["elapsed"]] - t0)
}
