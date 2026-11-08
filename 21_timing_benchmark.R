###############################################################################
# 21_timing_benchmark.R -- computing cost of all core pipelines measured in ONE
# session on ONE machine (one core, OPENBLAS_NUM_THREADS=1, same BLAS, same
# convergence criteria), with 3 repetitions, peak R memory and a breakdown into
# grouping (phenotype sorting), accumulation (Stage 1 or cross-products),
# estimation + solution, and prediction of the testing individuals.
# Settings: (A) p = 3,000, n = 2,000 ... 20,000 (500 new individuals);
#           (B) n = 8,000, p = 1,000 ... 8,000; (C) first simulated replicate of
#           S1 and first real-data partition. Q = 1,000; YQ-FEN k = 500, h = 200.
# The n x n one-stage GBLUP is measured only up to n = 6,000 (A) and in (C).
# Output: outputs/results/timing_benchmark.csv
###############################################################################
source("00_functions.R")
out <- file.path("..", "outputs", "results", "timing_benchmark.csv")
Q <- 1000
mem0 <- function() { g <- gc(reset = TRUE); sum(g[, 2]) }
memp <- function(base) { g <- gc(); sum(g[, 6]) - base }
pipe <- list(
  FULL = function(Z, y, Zn) { t <- c(); t0 <- tnow(); f <- fit_wgblup(Z, y, rep(1, nrow(Z))); t["estimate"] <- tnow() - t0
    t0 <- tnow(); p <- predict_wgblup(f, Zn); t["predict"] <- tnow() - t0; t },
  EXACT = function(Z, y, Zn) { t <- c(); t0 <- tnow(); cp <- crossprod_blocks(Z, y); t["accumulate"] <- tnow() - t0
    t0 <- tnow(); f <- fit_from_crossprods(cp$ZtZ, cp$Zty, cp$s1, cp$sy, cp$syy, cp$n, cp$p); t["estimate"] <- tnow() - t0
    t0 <- tnow(); p <- predict_lin(f, Zn); t["predict"] <- tnow() - t0; t },
  DIRECT_CG = function(Z, y, Zn) { t <- c(); t0 <- tnow(); cp <- crossprod_blocks(Z, y); t["accumulate"] <- tnow() - t0
    t0 <- tnow(); n <- cp$n; zb <- cp$s1 / n; yb <- cp$sy / n
    Tt <- list(Tc = cp$ZtZ - n * tcrossprod(zb), sc = cp$Zty - n * zb * yb, yy = cp$syy - n * yb^2, zbar = zb, ybar = yb, n = n, p = cp$p)
    mo <- moments_delta(Tt); s <- pcg_nystrom(Tt$Tc, Tt$sc, Tt$p * mo$delta, NULL); t["estimate"] <- tnow() - t0
    t0 <- tnow(); p <- drop(yb - sum(zb * s$x) + Zn %*% s$x); t["predict"] <- tnow() - t0; t },
  YQFE = function(Z, y, Zn) { t <- c(); t0 <- tnow(); lab <- make_groups("YQ", y, Q); t["group"] <- tnow() - t0
    t0 <- tnow(); S1 <- yq_fe_stage1(Z, y, lab); t["accumulate"] <- tnow() - t0
    t0 <- tnow(); f <- fit_yq_fe(S1); t["estimate"] <- tnow() - t0
    t0 <- tnow(); p <- predict_lin(f, Zn); t["predict"] <- tnow() - t0; t },
  YQFEN = function(Z, y, Zn) { t <- c(); t0 <- tnow(); lab <- make_groups("YQ", y, Q); t["group"] <- tnow() - t0
    t0 <- tnow(); S1 <- yq_fe_stage1(Z, y, lab); t["accumulate"] <- tnow() - t0
    t0 <- tnow(); f <- fit_yq_fen(S1, k = 500, h = 200); t["estimate"] <- tnow() - t0
    t0 <- tnow(); p <- predict_lin(f, Zn); t["predict"] <- tnow() - t0; t },
  YQ = function(Z, y, Zn) { t <- c(); t0 <- tnow(); lab <- make_groups("YQ", y, Q); t["group"] <- tnow() - t0
    t0 <- tnow(); S <- summarise_groups(Z, y, lab); t["accumulate"] <- tnow() - t0
    t0 <- tnow(); f <- fit_wgblup(S$Zbar, S$ybar, S$m); t["estimate"] <- tnow() - t0
    t0 <- tnow(); p <- predict_wgblup(f, Zn); t["predict"] <- tnow() - t0; t },
  NYS = function(Z, y, Zn) { t <- c(); t0 <- tnow(); f <- fit_nystrom(Z, y, Q, 1); t["estimate"] <- tnow() - t0
    t0 <- tnow(); p <- predict_nystrom(f, Zn); t["predict"] <- tnow() - t0; t })
run <- function(setting, n, p, Z, y, Zn, methods, reps = 3) {
  rows <- list()
  for (m in methods) for (r in seq_len(reps)) {
    base <- mem0(); t <- pipe[[m]](Z, y, Zn); mp <- memp(base)
    rows[[length(rows) + 1]] <- data.frame(setting = setting, n = n, p = p, method = m, rep = r,
      group = if ("group" %in% names(t)) t[["group"]] else 0, accumulate = if ("accumulate" %in% names(t)) t[["accumulate"]] else 0,
      estimate = t[["estimate"]], predict = t[["predict"]], total = sum(t), peak_mem_MB = mp)
    gc(verbose = FALSE)
  }
  x <- do.call(rbind, rows); old <- if (file.exists(out)) read.csv(out) else NULL
  write.csv(rbind(old, x), out, row.names = FALSE); print(aggregate(total ~ method, x, median))
}
core <- c("EXACT", "YQFE", "YQFEN", "DIRECT_CG", "YQ", "NYS")
## (C) data of the study
d <- make_sim_rep(1, 0)
run("C_sim", length(d$tr), ncol(d$Z), d$Z[d$tr, ], d$y[d$tr], d$Z[d$te, ], c("FULL", core)); rm(d); gc()
D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds")); pt <- readRDS(file.path("..", "outputs", "results", "real_partition_0.rds"))
run("C_real", length(pt$tr), ncol(D$Z), D$Z[pt$tr, ], D$y[pt$tr], D$Z[pt$te, ], c("FULL", core)); rm(D); gc()
## (A) n varies, p = 3,000
for (n in c(2000, 4000, 6000, 8000, 12000, 20000)) {
  set.seed(n); Z <- matrix(rnorm(n * 3000), n, 3000); y <- rowSums(Z[, 1:50]) * 0.3 + rnorm(n)
  run("A_n", n, 3000, Z, y, Z[1:500, ], c(if (n <= 4000) "FULL", core))
  if (n == 6000) run("A_n", n, 3000, Z, y, Z[1:500, ], "FULL", reps = 1)
  rm(Z); gc() }
## (B) p varies, n = 8,000
for (P in c(1000, 2000, 4000, 6000, 8000)) {
  set.seed(P); Z <- matrix(rnorm(8000 * P), 8000, P); y <- rowSums(Z[, 1:50]) * 0.3 + rnorm(8000)
  run("B_p", 8000, P, Z, y, Z[1:500, ], core); rm(Z); gc() }
