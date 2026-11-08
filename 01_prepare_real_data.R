###############################################################################
# 01_prepare_real_data.R
# Pool the six wheat data sets (W1-W6): same 5,740 markers, grain yield (GY).
# Markers with minor-allele frequency > 0.05 in the pooled data are kept and
# standardized; GY is standardized within each data set (W1-W6) to remove
# differences in mean and scale among them (environment/cohort effect).
###############################################################################
library(data.table)
raw <- file.path("..", "data", "raw")
Xs <- list(); ys <- list(); env <- c(); ids <- c()
for (k in 1:6) {
  X <- fread(file.path(raw, sprintf("X_W%d.csv", k)))
  Y <- fread(file.path(raw, sprintf("Y_W%d.csv", k)))
  stopifnot(all(X[[1]] == Y$Line))                 # same order of lines
  Xs[[k]] <- as.matrix(X[, -1])
  ys[[k]] <- (Y$GY - mean(Y$GY)) / sd(Y$GY)
  env <- c(env, rep(sprintf("W%d", k), nrow(Y))); ids <- c(ids, Y$Line)
  summ_k <- data.frame(set = sprintf("W%d", k), lines = nrow(Y), mean_GY = mean(Y$GY), sd_GY = sd(Y$GY))
  summ <- if (k == 1) summ_k else rbind(summ, summ_k)
  cat(sprintf("W%d: %d lines, mean GY = %.3f, SD = %.3f\n", k, nrow(Y), mean(Y$GY), sd(Y$GY)))
}
X <- do.call(rbind, Xs); y <- unlist(ys); rm(Xs)
stopifnot(!any(duplicated(ids)))
f <- colMeans(X) / 2; keep <- pmin(f, 1 - f) > 0.05
X <- X[, keep]
source("00_functions.R")
Z <- std_markers(X); rm(X)
saveRDS(list(Z = Z, y = y, env = env, id = ids), file.path("..", "data", "processed", "real_pooled.rds"))
cat("Pooled data:", nrow(Z), "lines x", ncol(Z), "markers (of", length(keep), ")\n")
fwrite(summ, file.path("..", "outputs", "tables", "real_data_summary.csv"))
