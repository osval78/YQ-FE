###############################################################################
# 23_revision_summary.R -- summaries of the revision analyses (20, 21, 22):
#  * Table R1: by data group (S1-S7, real random partitions, leave-one-set-out)
#    and method at Q = 1,000: % of the one-stage test correlation, calibration
#    slope, MSE ratio, correlation with true genetic values, RMS difference from
#    the one-stage predictions (in SD units) and overlap of the top 10%.
#  * Table R2: pre-specified paired comparisons with intervals that respect the
#    sampling unit (t over replicates within each scenario; corrected resampled t
#    of Nadeau and Bengio (2003) for the overlapping random partitions; t over the
#    six held-out data sets), Holm adjustment within each data group, and the
#    two one-sided tests (TOST) of equivalence of YQ-FEN and the one-stage GBLUP
#    with margin +-0.005 in correlation.
#  * Table R3: YQ-FEN sensitivity to k, h and the landmark seed.
#  * Table R4: timing benchmark (median and range over 3 repetitions, peak memory).
# Outputs: outputs/tables/TableR1_by_group.csv ... TableR4_timing.csv, rev_key_numbers.txt
###############################################################################
suppressPackageStartupMessages(library(data.table))
res <- file.path("..", "outputs", "results"); tabd <- file.path("..", "outputs", "tables")
X <- rbindlist(lapply(c("orig_sim", "orig_real", "new_sim", "lodo"), function(p) {
  f <- file.path(res, sprintf("rev_%s.csv", p)); if (file.exists(f)) fread(f, colClasses = list(character = "rep")) }), fill = TRUE)
X[, group := ifelse(scenario == "real_wheat", "Real (random partitions)", ifelse(scenario == "lodo", "Real (leave one set out)", scenario))]
ex <- X[method == "EXACT", .(group, rep, r_ex = r, mse_ex = mse)]
X <- merge(X, ex, by = c("group", "rep"))
X[, `:=`(pct = 100 * r / r_ex, mse_ratio = mse / mse_ex)]
R1 <- X[, .(n = .N, pct = mean(pct), pct_min = min(pct), slope = mean(slope), mse_ratio = mean(mse_ratio), r_g = mean(r_g),
            rel_rms = mean(rel_rms), top10 = mean(top10_overlap), time = mean(time), r = mean(r)), by = .(group, method)]
fwrite(R1, file.path(tabd, "TableR1_by_group.csv"))
## paired comparisons
ntr <- 6434; nte <- 1608
ci_fun <- function(d, grp) { J <- length(d); m <- mean(d); s2 <- var(d)
  f <- if (grp == "Real (random partitions)") (1 / J + nte / ntr) else 1 / J
  se <- sqrt(f * s2); tq <- qt(0.975, J - 1); tt <- m / se
  c(diff = m, lo = m - tq * se, hi = m + tq * se, t = tt, df = J - 1, P = 2 * pt(-abs(tt), J - 1), se = se) }
cmp <- list(c("EXACT", "YQ"), c("EXACT", "NYS"), c("EXACT", "YQ_CAL"), c("NYS_CAL", "YQ_CAL"), c("YQ_CVCAL", "YQ_CAL"), c("EXACT", "YQFEN"))
R2 <- NULL
for (g in unique(X$group)) {
  W <- dcast(X[group == g], rep ~ method, value.var = "r")
  for (c2 in cmp) { if (!all(c2 %in% names(W))) next
    v <- ci_fun(W[[c2[1]]] - W[[c2[2]]], g)
    R2 <- rbind(R2, data.table(group = g, comparison = paste(c2[1], "minus", c2[2]), t(v))) }
}
R2[, P_holm := p.adjust(P, "holm"), by = group]
## TOST for YQ-FEN vs one-stage: H0 |diff| >= 0.005
R2[comparison == "EXACT minus YQFEN", `:=`(P_tost = pmax(pt((diff + 0.005) / se, df, lower.tail = FALSE), pt((diff - 0.005) / se, df)),
                                           equivalent = (lo > -0.005 & hi < 0.005))]
fwrite(R2, file.path(tabd, "TableR2_paired.csv"))
## sensitivity
S <- rbindlist(lapply(c("orig_sim", "orig_real"), function(p) { f <- file.path(res, sprintf("rev_sens_%s.csv", p)); if (file.exists(f)) fread(f, colClasses = list(character = "rep")) }))
S[, data := ifelse(scenario == "real_wheat", "Real", "Simulation")]
R3 <- S[, .(n = .N, pct = mean(100 * r / r_onestage), pct_min = min(100 * r / r_onestage), mae_log_delta = mean(abs(log(delta / delta_reml))),
            max_ratio_delta = max(exp(abs(log(delta / delta_reml)))), max_abs_diff = max(max_abs_diff), rel_rms = mean(rel_rms),
            top10 = mean(top10_overlap), iter = mean(iter), time = mean(time_nystrom + time_rest), truncated = sum(truncated),
            rel_det_max = max(rel_det)), by = .(data, k, h, seed_set = ifelse(seed == 1, "seed 1", "seeds 2-5"))]
fwrite(R3, file.path(tabd, "TableR3_yqfen_sensitivity.csv"))
## timing
if (file.exists(file.path(res, "timing_benchmark.csv"))) {
  Tm <- fread(file.path(res, "timing_benchmark.csv"))
  R4 <- Tm[, .(total_median = median(total), total_min = min(total), total_max = max(total), group = median(group), accumulate = median(accumulate),
               estimate = median(estimate), predict = median(predict), peak_mem_MB = max(peak_mem_MB), reps = .N), by = .(setting, n, p, method)]
  fwrite(R4, file.path(tabd, "TableR4_timing.csv"))
} else R4 <- NULL
da <- if (file.exists(file.path(res, "distributed_aggregation.csv"))) fread(file.path(res, "distributed_aggregation.csv")) else NULL
options(width = 250)
capture.output({ print(R1, digits = 4); print(R2, digits = 3); print(R3, digits = 4); print(R4, digits = 4); print(da) },
               file = file.path(tabd, "rev_key_numbers.txt"))
cat("revision summary written\n")
