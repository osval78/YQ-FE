###############################################################################
# 24_real_data_description.R -- description of the six wheat data sets used in
# the pooled analysis: lines, missing marker calls in the raw files, coding,
# (non-integer dosages = values imputed before distribution),
# genomic relationships within and between data sets (mean off-diagonal G and,
# for every line, its largest relationship with lines of its own set and of the
# other sets). Output: outputs/tables/real_data_description.csv
###############################################################################
suppressPackageStartupMessages(library(data.table)); source("00_functions.R")
raw <- file.path("..", "data", "raw"); out <- NULL
for (k in 1:6) { X <- fread(file.path(raw, sprintf("X_W%d.csv", k))); M <- as.matrix(X[, -1])
  out <- rbind(out, data.frame(set = sprintf("W%d", k), lines = nrow(M), markers = ncol(M), missing_calls = sum(is.na(M)),
    prop_non_integer_dosages = mean(abs(M - round(M)) > 1e-8), n_distinct_values = length(unique(as.vector(M))))) }
D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds")); p <- ncol(D$Z); rel <- NULL
for (w in sort(unique(D$env))) { i <- which(D$env == w); G <- tcrossprod(D$Z[i, ], D$Z) / p
  own <- D$env == w; Gi <- G[, own]; diag(Gi) <- NA; Go <- G[, !own]
  rel <- rbind(rel, data.frame(set = w, mean_G_within = mean(Gi, na.rm = TRUE), mean_G_between = mean(Go),
    mean_max_G_within = mean(apply(Gi, 1, max, na.rm = TRUE)), mean_max_G_other_sets = mean(apply(Go, 1, max)))) }
out <- merge(out, rel, by = "set"); print(out)
fwrite(out, file.path("..", "outputs", "tables", "real_data_description.csv"))
