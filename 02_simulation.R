###############################################################################
# 02_simulation.R  -- simulation study (resumable)
# 3 scenarios x 5 replicates; N = 8,000 individuals (6,400 training,
# 1,600 testing), p = 3,000 markers, h2 = 0.5, Q = 250, 1000, 2000, 3000.
###############################################################################
source("00_functions.R")
out <- file.path("..", "outputs", "results", "sim_results.csv")
rows <- if (file.exists(out)) read.csv(out) else NULL
done <- if (is.null(rows)) character(0) else unique(paste(rows$scenario, rows$rep))
t00 <- tnow()
for (si in seq_along(SIM$SCEN)) for (rep in seq_len(SIM$NREP) - 1) {
  sc <- SIM$SCEN[si]
  if (paste(sc, rep) %in% done) next
  message(sprintf("%s rep %d", sc, rep + 1))
  d <- make_sim_rep(si, rep)
  res <- evaluate_rep(d$Z[d$tr, ], d$y[d$tr], d$Z[d$te, ], d$y[d$te], d$g[d$te],
                      SIM$QS, seed = d$seed)
  res$scenario <- sc; res$rep <- rep; res$n_train <- length(d$tr); res$p <- SIM$P
  rows <- rbind(rows, res); write.csv(rows, out, row.names = FALSE)
  rm(d); gc(verbose = FALSE)
  message(sprintf("  finished [%.0f s total]", tnow() - t00))
}
