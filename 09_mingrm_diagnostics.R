###############################################################################
# 09_mingrm_diagnostics.R -- describes the training sets chosen by Min_GRM
# (Section 2.9) and compares them with random subsamples of the same size:
#   * check of the rule: the kept individuals are exactly the Q candidates with
#     the largest s_i = min_{j in TS} G_ij (no candidate outside the set has a
#     larger score than any candidate inside it);
#   * mean diagonal element of G = ZZ'/p of the kept individuals;
#   * mean over testing individuals of their largest and of their average
#     genomic relationship with the kept individuals;
#   * mean genomic relationship among the kept individuals (diversity);
#   * real data: number of kept lines from each data set (W1-W6).
# Data: the five real-data partitions and the first replicate of each simulated
# scenario; Q = 250, 1000, 2000, 3000. The random subsample uses the same seeds
# as evaluate_rep() (set.seed(seed + Q)).
# Output: outputs/results/mingrm_diagnostics.csv, mingrm_by_dataset.csv
###############################################################################
source("00_functions.R")
QS <- c(250, 1000, 2000, 3000)
desc <- function(Ztr, Zte, idx, sc) {
  p <- ncol(Ztr); Zk <- Ztr[idx, , drop = FALSE]
  Gtk <- tcrossprod(Zte, Zk) / p                       # n_TS x Q
  Gkk <- tcrossprod(Zk) / p
  c(mean_diag_kept = mean(diag(Gkk)),
    mean_offdiag_kept = (sum(Gkk) - sum(diag(Gkk))) / (length(idx) * (length(idx) - 1)),
    mean_max_rel_test = mean(apply(Gtk, 1, max)),
    mean_avg_rel_test = mean(Gtk),
    mean_score_kept = mean(sc[idx]))
}
rows <- list(); drows <- list()
run_one <- function(data, rep, Ztr, Zte, seed, env = NULL) {
  sc <- mingrm_score(Ztr, Zte); n <- nrow(Ztr)
  for (Q in QS) {
    im <- select_mingrm(Ztr, Zte, Q, sc)
    ok <- min(sc[im]) >= max(sc[-im])                   # top-Q rule
    set.seed(seed + Q); ir <- sample.int(n, Q)
    rows[[length(rows) + 1]] <<- rbind(
      data.frame(data = data, rep = rep, Q = Q, method = "MINGRM", top_Q_rule = ok, t(desc(Ztr, Zte, im, sc))),
      data.frame(data = data, rep = rep, Q = Q, method = "SUB", top_Q_rule = NA, t(desc(Ztr, Zte, ir, sc))))
    if (!is.null(env)) drows[[length(drows) + 1]] <<- data.frame(rep = rep, Q = Q,
      t(c(table(factor(env[im], levels = sort(unique(env)))))), all = Q)
  }
  message(sprintf("%s rep %d done", data, rep))
}
## real data
D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds"))
for (rep in 0:4) {
  pt <- readRDS(file.path("..", "outputs", "results", sprintf("real_partition_%d.rds", rep)))
  # seed of the random subsample: the one stored with the partition (as in 03_real_data.R)
  run_one("real_wheat", rep, D$Z[pt$tr, ], D$Z[pt$te, ], pt$seed, D$env[pt$tr])
}
## simulation: first replicate of each scenario
for (si in seq_along(SIM$SCEN)) {
  d <- make_sim_rep(si, 0)
  run_one(SIM$SCEN[si], 0, d$Z[d$tr, ], d$Z[d$te, ], d$seed)
}
out <- do.call(rbind, rows); print(out)
write.csv(out, file.path("..", "outputs", "results", "mingrm_diagnostics.csv"), row.names = FALSE)
dd <- do.call(rbind, drows); envtr <- table(D$env) / length(D$env); print(dd); print(round(envtr, 3))
write.csv(dd, file.path("..", "outputs", "results", "mingrm_by_dataset.csv"), row.names = FALSE)
