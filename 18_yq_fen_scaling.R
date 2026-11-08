###############################################################################
# 18_yq_fen_scaling.R -- computing time of YQ-FEN in the two timing experiments
# of 04_scaling.R (same simulated matrices, same seeds): (A) p = 3,000 markers,
# n = 2,000 ... 20,000; (B) n = 8,000, p = 1,000 ... 8,000; Q = 1,000 groups,
# k = 500 landmark markers, 500 new individuals predicted. YQ-FE is re-timed in
# the same session (same machine, one core), so the two columns are comparable;
# the columns t_yqfe and t_yqfen (and the stage split of YQ-FEN) are written
# into scale_results.csv and scale_p_results.csv.
###############################################################################
source("00_functions.R")
Q <- 1000
time_both <- function(Z, y, Znew) {
  t0 <- tnow(); S1 <- yq_fe_stage1(Z, y, make_groups("YQ", y, Q)); ts1 <- tnow() - t0
  t0 <- tnow(); f <- fit_yq_fe(S1); pr <- predict_lin(f, Znew); tfe <- ts1 + tnow() - t0
  t0 <- tnow(); g <- fit_yq_fen(S1, k = 500); pr2 <- predict_lin(g, Znew); tfen2 <- tnow() - t0
  c(t_yqfe = tfe, t_yqfen = ts1 + tfen2, t_stage1 = ts1, t_yqfen_stage2 = tfen2,
    maxdiff = max(abs(pr - pr2)), cor = cor(pr, pr2))
}
out <- file.path("..", "outputs", "results", "scale_results.csv"); A <- read.csv(out); P <- 3000
for (i in seq_len(nrow(A))) {
  n <- A$n[i]; set.seed(n); Z <- matrix(rnorm(n * P), n, P); y <- rowSums(Z[, 1:50]) * 0.3 + rnorm(n); Znew <- Z[1:500, ]
  r <- time_both(Z, y, Znew); for (nm in names(r)) A[i, nm] <- r[[nm]]
  print(A[i, c("n", names(r))]); write.csv(A, out, row.names = FALSE); rm(Z); gc(verbose = FALSE)
}
out <- file.path("..", "outputs", "results", "scale_p_results.csv"); B <- read.csv(out); n <- 8000
for (i in seq_len(nrow(B))) {
  P <- B$p[i]; set.seed(P); Z <- matrix(rnorm(n * P), n, P); y <- rowSums(Z[, 1:50]) * 0.3 + rnorm(n); Znew <- Z[1:500, ]
  r <- time_both(Z, y, Znew); for (nm in names(r)) B[i, nm] <- r[[nm]]
  print(B[i, c("p", names(r))]); write.csv(B, out, row.names = FALSE); rm(Z); gc(verbose = FALSE)
}
