###############################################################################
# 15_yq_fe_scatter.R -- predictions of the testing individuals from the one-stage
# GBLUP, YQ-FE, YQ-FEN (k = 500) and YQ (Q = 1,000; first simulated replicate of S1
# and first real-data partition), used in Figure 12.
# Output: outputs/results/yq_fe_scatter.csv
###############################################################################
source("00_functions.R")
D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds"))
pt <- readRDS(file.path("..", "outputs", "results", "real_partition_0.rds"))
ds <- make_sim_rep(1, 0)
L <- list(`Simulation (S1, replicate 1)` = list(ds$Z[ds$tr, ], ds$y[ds$tr], ds$Z[ds$te, ], ds$y[ds$te]),
          `Real data (partition 1)` = list(D$Z[pt$tr, ], D$y[pt$tr], D$Z[pt$te, ], D$y[pt$te]))
rm(D, ds); gc(verbose = FALSE)
out <- NULL
for (nm in names(L)) {
  X <- L[[nm]]; Q <- 1000; lab <- make_groups("YQ", X[[2]], Q)
  fm <- fit_marker_gblup(X[[1]], X[[2]]); p1 <- predict_marker_gblup(fm, X[[3]])
  S1 <- yq_fe_stage1(X[[1]], X[[2]], lab)
  f <- fit_yq_fe(S1); p2 <- predict_lin(f, X[[3]])
  fn <- fit_yq_fen(S1, k = 500); p4 <- predict_lin(fn, X[[3]])
  S <- summarise_groups(X[[1]], X[[2]], lab); fy <- fit_wgblup(S$Zbar, S$ybar, S$m); p3 <- predict_wgblup(fy, X[[3]])
  out <- rbind(out, data.frame(data = nm, one_stage = p1, yq_fe = p2, yq_fen = p4, yq = p3, y = X[[4]]))
  message(nm, " done"); gc(verbose = FALSE)
}
write.csv(out, file.path("..", "outputs", "results", "yq_fe_scatter.csv"), row.names = FALSE)
