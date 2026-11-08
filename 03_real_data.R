###############################################################################
# 03_real_data.R -- pooled wheat data W1-W6 (n = 8,042), 5 random 80/20
# partitions, Q = 250, 1000, 2000, 3000 (resumable)
###############################################################################
source("00_functions.R")
D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds"))
Z <- D$Z; y <- D$y; n <- length(y)
QS <- c(250, 1000, 2000, 3000); NREP <- 5
out <- file.path("..", "outputs", "results", "real_results.csv")
rows <- if (file.exists(out)) read.csv(out) else NULL
done <- if (is.null(rows)) integer(0) else unique(rows$rep)
for (rep in seq_len(NREP) - 1) {
  if (rep %in% done) next
  set.seed(2025 * 1000 + 990 + rep)
  perm <- sample.int(n); nte <- round(0.2 * n)
  te <- perm[seq_len(nte)]; tr <- perm[-seq_len(nte)]; seed <- sample.int(1e8, 1)
  message(sprintf("real partition %d", rep + 1))
  res <- evaluate_rep(Z[tr, ], y[tr], Z[te, ], y[te], NULL, QS, seed = seed)
  res$scenario <- "real_wheat"; res$rep <- rep; res$n_train <- length(tr); res$p <- ncol(Z)
  rows <- rbind(rows, res); write.csv(rows, out, row.names = FALSE)
  saveRDS(list(tr = tr, te = te, seed = seed), file.path("..", "outputs", "results", sprintf("real_partition_%d.rds", rep)))
  gc(verbose = FALSE)
}
