###############################################################################
# 22_distributed_aggregation.R -- practical value of the two-stage organization:
# each of the six wheat data sets computes its own Stage 1 summaries (its own
# phenotype-quantile groups); the summaries are combined in Stage 2 (YQ-FE and
# YQ-FEN) without sharing individual records. Predictions of the first real-data
# partition's testing lines are compared with those of the pooled one-stage GBLUP.
# Also checks that YQ-FE and YQ-FEN do not depend on Q.
# Output: outputs/results/distributed_aggregation.csv
###############################################################################
source("00_functions.R")
D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds")); pt <- readRDS(file.path("..", "outputs", "results", "real_partition_0.rds"))
tr <- pt$tr; te <- pt$te; fm <- fit_marker_gblup(D$Z[tr, ], D$y[tr]); pm <- predict_marker_gblup(fm, D$Z[te, ])
combine <- function(parts) {   # Stage 1 summaries of several sources -> one summary object
  list(Zbar = do.call(rbind, lapply(parts, `[[`, "Zbar")), ybar = unlist(lapply(parts, `[[`, "ybar")), m = unlist(lapply(parts, `[[`, "m")),
       W = Reduce(`+`, lapply(parts, `[[`, "W")), sW = Reduce(`+`, lapply(parts, `[[`, "sW")), tW = sum(sapply(parts, `[[`, "tW")),
       n = sum(sapply(parts, `[[`, "n")), Q = sum(sapply(parts, `[[`, "Q")), p = parts[[1]]$p) }
res <- NULL
for (Qs in c(100, 170)) {   # groups per data set (about 600 and 1,000 in total)
  parts <- lapply(sort(unique(D$env)), function(w) { i <- tr[D$env[tr] == w]; yq_fe_stage1(D$Z[i, ], D$y[i], make_groups("YQ", D$y[i], Qs)) })
  S1 <- combine(parts)
  fe <- fit_yq_fe(S1); pfe <- predict_lin(fe, D$Z[te, ]); fn <- fit_yq_fen(S1, 500, 200); pfn <- predict_lin(fn, D$Z[te, ])
  res <- rbind(res, data.frame(groups_per_set = Qs, total_groups = S1$Q, method = c("YQ-FE", "YQ-FEN"),
    max_abs_diff_vs_pooled_onestage = c(max(abs(pfe - pm)), max(abs(pfn - pm))), delta = c(fe$delta, fn$delta), delta_onestage = fm$delta))
}
for (Q in c(250, 3000)) {
  S1 <- yq_fe_stage1(D$Z[tr, ], D$y[tr], make_groups("YQ", D$y[tr], Q)); fe <- fit_yq_fe(S1); fn <- fit_yq_fen(S1, 500, 200)
  res <- rbind(res, data.frame(groups_per_set = NA, total_groups = Q, method = c("YQ-FE (pooled)", "YQ-FEN (pooled)"),
    max_abs_diff_vs_pooled_onestage = c(max(abs(predict_lin(fe, D$Z[te, ]) - pm)), max(abs(predict_lin(fn, D$Z[te, ]) - pm))),
    delta = c(fe$delta, fn$delta), delta_onestage = fm$delta))
}
print(res); write.csv(res, file.path("..", "outputs", "results", "distributed_aggregation.csv"), row.names = FALSE)
