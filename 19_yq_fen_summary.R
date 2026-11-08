###############################################################################
# 19_yq_fen_summary.R -- key numbers of YQ-FEN used in the text: accuracy,
# calibration and MSE relative to the one-stage GBLUP, moment versus REML
# estimates of delta, paired t-tests against YQ-FE, YQ and Nystrom, iterations
# and time by number of landmark markers k, and the timing experiments.
# Output: outputs/tables/yq_fen_key_numbers.txt, TableC10_yq_fen_paired.csv
###############################################################################
suppressPackageStartupMessages(library(data.table))
res <- file.path("..", "outputs", "results"); tabd <- file.path("..", "outputs", "tables")
rd <- function(d) {
  sim <- fread(file.path(res, sprintf("%s_results.csv", c(sim = "sim", real = "real")[d])))
  fen <- fread(file.path(res, sprintf("yq_fen_%s_results.csv", d)))
  fe <- fread(file.path(res, sprintf("yq_fe_%s_results.csv", d)))
  if (d == "real") { sim[, scenario := "real_wheat"]; sim <- sim[rep < 5] } else sim <- sim[rep < 5]
  full <- sim[method == "FULL", .(scenario, rep, r_full = r, mse_full = mse)]
  list(sim = sim, fen = merge(fen, full, by = c("scenario", "rep")), fe = fe, full = full)
}
out <- list()
pt <- function(x, y) { d <- x - y; tt <- t.test(d); c(diff = mean(d), lo = tt$conf.int[1], hi = tt$conf.int[2], t = unname(tt$statistic), df = unname(tt$parameter), P = tt$p.value) }
P <- NULL
for (d in c("sim", "real")) {
  X <- rd(d); f <- X$fen
  f[, `:=`(pct = 100 * r / r_full, mse_ratio = mse / mse_full, lr = log(delta / delta_onestage))]
  out[[paste0(d, "_by_variant")]] <- f[, .(n = .N, pct = mean(pct), min_pct = min(pct), slope = mean(slope), mse_ratio = mean(mse_ratio),
       r_g = mean(r_g), delta = mean(delta), delta_reml = mean(delta_onestage), mae_log_delta = mean(abs(lr)), max_ratio_delta = max(exp(abs(lr))),
       max_abs_diff_pred = max(max_abs_diff_pred_vs_onestage), med_abs_diff_pred = median(max_abs_diff_pred_vs_onestage),
       iter = mean(iter), time = mean(time), time_stage1 = mean(time_stage1), time_stage2 = mean(time_stage2),
       time_nystrom = mean(time_nystrom), time_moments = mean(time_moments), time_cg = mean(time_cg),
       time_yqfe = mean(time_yqfe), max_diff_vs_yqfe = max(max_abs_diff_vs_yqfe)), by = .(variant, Q)]
  k5 <- f[variant == "k500" & Q == 1000, .(scenario, rep, r_fen = r)]
  sm <- merge(X$sim[shrink == "REML" & method %in% c("YQ", "NYS", "SQC", "SUB"), .(scenario, rep, method, Q, r)], k5, by = c("scenario", "rep"))
  fe1 <- merge(X$fe[Q == 1000, .(scenario, rep, r_fe = r)], k5, by = c("scenario", "rep"))
  P <- rbind(P, data.table(data = d, comparison = "YQ-FE minus YQ-FEN", Q = NA, t(pt(fe1$r_fe, fe1$r_fen))))
  for (m in c("YQ", "NYS", "SQC", "SUB")) for (q in c(250, 1000, 2000, 3000)) {
    z <- sm[method == m & Q == q]; P <- rbind(P, data.table(data = d, comparison = paste("YQ-FEN minus", m), Q = q, t(pt(z$r_fen, z$r)))) }
  out[[paste0(d, "_by_scenario")]] <- f[variant == "k500" & Q == 1000, .(pct = mean(pct), r_g = mean(r_g), max_ratio_delta = max(exp(abs(lr)))), by = scenario]
}
fwrite(P, file.path(tabd, "TableC10_yq_fen_paired.csv")); out$paired <- P
sc <- fread(file.path(res, "scale_results.csv")); scp <- fread(file.path(res, "scale_p_results.csv"))
out$scale_n <- sc[, intersect(c("n", "t_exact", "t_yqfe", "t_yqfen", "t_stage1", "t_yqfen_stage2", "t_yq", "t_nys", "maxdiff", "cor"), names(sc)), with = FALSE]
out$scale_p <- scp[, intersect(c("p", "t_exact", "t_yqfe", "t_yqfen", "t_stage1", "t_yqfen_stage2", "t_yq", "t_nys", "maxdiff", "cor"), names(scp)), with = FALSE]
options(width = 250); capture.output(print(out, digits = 4), file = file.path(tabd, "yq_fen_key_numbers.txt"))
cat("YQ-FEN summary written\n")
