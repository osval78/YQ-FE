###############################################################################
# 07_figures_tables.R -- all figures (ggplot2, 350 dpi) and result tables
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(data.table) })
res <- file.path("..", "outputs", "results"); figd <- file.path("..", "outputs", "figures")
tabd <- file.path("..", "outputs", "tables")
sim <- fread(file.path(res, "sim_results.csv"))[rep < 5]; real <- fread(file.path(res, "real_results.csv"))
## methods added in the revision (08_added_methods.R): MINGRM (Min_GRM selection rule),
## EXACT (marker-space GBLUP, identical to FULL) and MARG (marginal regression)
asim <- fread(file.path(res, "added_sim_results.csv")); areal <- fread(file.path(res, "added_real_results.csv"))
cols0 <- names(sim)
sim <- rbind(sim, asim[method == "MINGRM", ..cols0]); real <- rbind(real, areal[method == "MINGRM", ..cols0])
## v8: fully efficient phenotype-quantile compression (14_yq_fe_all.R) as method "YQFE"
fe_as_method <- function(f) { x <- fread(f); x[, `:=`(cap_train = NA_real_, cap_test = NA_real_)]; x[, ..cols0] }
if (file.exists(file.path(res, "yq_fe_sim_results.csv"))) {
  sim <- rbind(sim, fe_as_method(file.path(res, "yq_fe_sim_results.csv")))
  real <- rbind(real, fe_as_method(file.path(res, "yq_fe_real_results.csv"))) }
## v10: Nystrom-accelerated YQ-FE (17_yq_fen_all.R) as method "YQFEN" (k = 500 landmark
## markers; variance components by moments; shown with the default variant of every method)
fen_as_method <- function(f) { x <- fread(f)[variant == "k500"]
  x[, `:=`(method = "YQFEN", shrink = "REML", groups = Q, cap_train = NA_real_, cap_test = NA_real_)]; x[, ..cols0] }
if (file.exists(file.path(res, "yq_fen_sim_results.csv"))) {
  sim <- rbind(sim, fen_as_method(file.path(res, "yq_fen_sim_results.csv")))
  real <- rbind(real, fen_as_method(file.path(res, "yq_fen_real_results.csv"))) }
LAB <- c(YQFE = "YQ-FE (proposed)", YQFEN = "YQ-FEN (proposed)", YQ = "YQ (proposed)",
         NYS = "Nystr\u00F6m (rank Q)", SQC = "Cluster \u00D7 quantile means (SQC)", GEN = "Genomic-cluster means",
         RND = "Random group means", SUB = "Random subsample", MINGRM = "Min_GRM selection (MINGRM)")
COL <- c(YQFE = "#00315E", YQFEN = "#C51B7D", YQ = "#0072B2", NYS = "#7B3294", SQC = "#D55E00", GEN = "#009E73",
         RND = "#E69F00", SUB = "#8C8C8C", MINGRM = "#CC79A7")
## marginal-regression reference (does not depend on Q)
margS <- asim[method == "MARG"]; margR <- areal[method == "MARG"]
SCL <- c(S1_polygenic = "S1 polygenic", S2_sparse = "S2 sparse", S3_structure = "S3 population structure")
QS <- c(250, 1000, 2000, 3000)
th <- theme_classic(base_size = 10) + theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10), legend.position = "top",
        legend.title = element_blank(), plot.tag = element_text(face = "bold"))
ci <- function(x) { x <- x[is.finite(x)]; n <- length(x); if (n < 2) return(c(mean(x), NA, NA))
  h <- qt(0.975, n - 1) * sd(x) / sqrt(n); c(mean(x), mean(x) - h, mean(x) + h) }
sumci <- function(dt, val, by) dt[, { v <- ci(get(val)); .(mean = v[1], lo = v[2], hi = v[3], nrep = .N) }, by = by]
save_fig <- function(p, name, w, h) ggsave(file.path(figd, name), p, width = w, height = h, dpi = 350, bg = "white", device = ragg::agg_png)

## percentage of the full model, replicate by replicate
addpct <- function(dt) { fr <- dt[method == "FULL", .(scenario, rep, r_full = r, mse_full = mse, t_full = time)]
  dt <- merge(dt[method != "FULL"], fr, by = c("scenario", "rep"))
  dt[, `:=`(pct = 100 * r / r_full, mse_ratio = mse / mse_full, speedup = t_full / time)]; dt }
simp <- addpct(sim); realp <- addpct(real)

## ---------------- Figure 1: workflow and grouping rules (toy) ---------------
bx <- data.frame(x = 1:5, lab = c("Training data\nn individuals\n(n \u00D7 p markers,\nn phenotypes)",
  "1. Grouping rule\nn \u2192 Q groups\n(Q \u226A n)", "2. Stage 1\ngroup means of\nmarkers and\nphenotype (BLUEs,\nVar = \u03C3\u00B2/m)",
  "3. Stage 2\nweighted GBLUP\non Q group means\n(lme4 REML\nor BGLR)", "4. Marker effects\npredict new\nindividuals"),
  fill = c("#2B6CB0", "#E69F00", "#009E73", "#7B3294", "#D55E00"))
pA <- ggplot(bx) + geom_tile(aes(x, 0, fill = fill), width = 0.88, height = 1, colour = NA) +
  geom_text(aes(x, 0, label = lab), colour = "white", size = 2.75, lineheight = 0.95, fontface = "bold") +
  annotate("segment", x = 1:4 + 0.42, xend = 2:5 - 0.42, y = 0, yend = 0,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  annotate("text", x = 3, y = -0.8, lineheight = 1.1, size = 3, fontface = "italic", colour = "grey25",
    label = "Testing individuals are never compressed or used for grouping.\nFull GBLUP: n \u00D7 n kernel, O(n\u00B3) time, O(n\u00B2) memory;  compressed model: Q \u00D7 Q kernel, O(Q\u00B3 + Q\u00B2p) after one pass over the data.") +
  scale_fill_identity() + coord_cartesian(xlim = c(0.5, 5.5), ylim = c(-1.05, 0.55)) + theme_void() +
  labs(title = "Compress \u2192 fit \u2192 predict: a two-stage analysis in which only Q group means enter the model") +
  theme(plot.title = element_text(face = "bold", size = 11))
set.seed(4); nt <- 60; pc1 <- rnorm(nt); yy <- 0.6 * pc1 + rnorm(nt, sd = 0.8)
Qt <- 6; toy <- data.frame(pc1, yy)
toy[, c("RND","GEN","YQ","SQC")] <- 0
toy$RND <- sample(rep_len(1:Qt, nt)); toy$GEN <- kmeans(pc1, Qt, nstart = 20)$cluster
toy$YQ <- 0; toy$YQ[order(yy)] <- ceiling(seq_len(nt) * Qt / nt)
cl3 <- kmeans(pc1, 3, nstart = 20)$cluster; toy$SQC <- 0
for (c in 1:3) { ix <- which(cl3 == c); ix <- ix[order(yy[ix])]; toy$SQC[ix] <- (c - 1) * 2 + ceiling(seq_along(ix) * 2 / length(ix)) }
toy[] <- lapply(toy, as.numeric); tl <- melt(as.data.table(toy), id.vars = c("pc1", "yy"))
tl[, variable := factor(variable, levels = c("RND", "GEN", "YQ", "SQC"),
  labels = c("Random groups", "Genomic clusters", "Quantile bins (YQ)", "Cluster \u00D7 quantile (SQC)"))]
pB <- ggplot(tl, aes(pc1, yy, colour = factor(value))) + geom_point(size = 1.8) + facet_wrap(~variable, nrow = 1) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9")) +
  labs(x = "genomic axis (PC1)", y = "phenotype", title = "Four grouping rules (toy example; colour = group)") +
  th + theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(),
             plot.title = element_text(face = "bold", size = 11))
save_fig(pA / pB + plot_layout(heights = c(1, 1.15)), "Figure1_workflow.png", 10, 6.2)

## ---------------- Figure 2: simulation, r versus Q ---------------------------
s2 <- sumci(sim[method != "FULL" & shrink == "REML"], "r", c("scenario", "method", "Q"))
fb <- sumci(sim[method == "FULL"], "r", "scenario")
s2[, method := factor(method, levels = names(LAB))]
p2 <- ggplot(s2, aes(Q, mean, colour = method)) +
  geom_rect(data = fb, aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey85", alpha = 0.8) +
  geom_hline(data = fb, aes(yintercept = mean), linetype = 2) +
  geom_hline(data = sumci(margS, "r", "scenario"), aes(yintercept = mean), linetype = 3, colour = "grey30") +
  geom_line(position = position_dodge(60)) + geom_point(size = 1.6, position = position_dodge(60)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 70, position = position_dodge(60)) +
  facet_wrap(~scenario, labeller = as_labeller(SCL), scales = "free_y") +
  scale_colour_manual(values = COL, labels = LAB) + scale_x_continuous(breaks = QS) +
  labs(x = "Number of groups / landmarks Q", y = "Test correlation") + th +
  guides(colour = guide_legend(nrow = 3))
save_fig(p2, "Figure2_sim_r_vs_Q.png", 10, 4.6)

## ---------------- Figure 3: % of full, simulation ----------------------------
s3 <- sumci(simp[shrink == "REML" & Q %in% c(250, 1000, 3000)], "pct", c("method", "Q"))
s3[, method := factor(method, levels = rev(names(LAB)))]
p3 <- ggplot(s3, aes(mean, method, colour = method)) + geom_vline(xintercept = 100, linetype = 2) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25) + geom_point(size = 2.2) +
  facet_wrap(~Q, labeller = labeller(Q = function(q) paste0("Q = ", format(as.numeric(q), big.mark = ",")))) +
  scale_colour_manual(values = COL, guide = "none") + scale_y_discrete(labels = LAB) +
  labs(x = "% of the full-GBLUP test correlation", y = NULL) + th
save_fig(p3, "Figure3_sim_pct_full.png", 10, 3.9)

## ---------------- Figure 4: calibration ---------------------------------------
cal <- rbind(cbind(simp[Q == 1000], data = "Simulation"), cbind(realp[Q == 1000], data = "Real data"))
cal[, data := factor(data, levels = c("Simulation", "Real data"))]
s4a <- sumci(cal[shrink == "REML"], "slope", c("data", "method"))
s4a[, method := factor(method, levels = rev(names(LAB)))]
p4a <- ggplot(s4a, aes(mean, method, shape = data)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = method), height = 0.2, position = position_dodge(0.5)) +
  geom_point(aes(colour = method), size = 2.3, position = position_dodge(0.5), fill = "white") +
  scale_shape_manual(values = c(Simulation = 16, `Real data` = 21)) +
  scale_colour_manual(values = COL, guide = "none") + scale_y_discrete(labels = LAB) +
  labs(x = "Calibration slope (observed on predicted), Q = 1,000", y = NULL, tag = "A") + th
s4b <- sumci(cal[method %in% c("GEN", "YQ", "SQC") & shrink %in% c("REML", "CAL")], "mse_ratio", c("data", "method", "shrink"))
s4b[, shrink := factor(shrink, levels = c("REML", "CAL"), labels = c("before recalibration", "after recalibration"))]
p4b <- ggplot(s4b, aes(method, mean, fill = shrink)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.8), width = 0.2) +
  geom_hline(yintercept = 1, linetype = 2) + facet_wrap(~data) +
  scale_fill_manual(values = c("#9ECAE1", "#08519C")) +
  labs(x = NULL, y = "MSE relative to full GBLUP", tag = "B") + th
save_fig(p4a + p4b + plot_layout(widths = c(1, 1.1)), "Figure4_calibration.png", 10.5, 4.4)

## ---------------- Figure 5: computing cost ------------------------------------
sc <- fread(file.path(res, "scale_results.csv"))
scl <- melt(sc, id.vars = "n", measure.vars = intersect(c("t_full", "t_exact", "t_yqfe", "t_yqfen", "t_yq", "t_sqc", "t_nys", "t_mingrm"), names(sc)), na.rm = TRUE)
scl[, variable := factor(variable, levels = c("t_full", "t_exact", "t_yqfe", "t_yqfen", "t_yq", "t_nys", "t_sqc", "t_mingrm"),
  labels = c("Full GBLUP (n \u00D7 n)", "Exact GBLUP, marker space (p \u00D7 p)", "YQ-FE, Q = 1,000", "YQ-FEN, Q = 1,000, k = 500",
             "YQ, Q = 1,000", "Nystr\u00F6m, Q = 1,000", "SQC, Q = 1,000", "Min_GRM, Q = 1,000"))]
fitpow <- lm(log(t_full) ~ log(n), data = sc[!is.na(t_full) & n >= 2000])
ext <- data.frame(n = c(8000, 12000, 20000)); ext$t <- exp(predict(fitpow, ext))
p5a <- ggplot(scl, aes(n, value, colour = variable)) + geom_line() + geom_point(size = 1.8) +
  geom_line(data = ext, aes(n, t), inherit.aes = FALSE, linetype = 3) +
  scale_x_log10(breaks = c(1000, 2000, 4000, 8000, 20000), labels = scales::comma) +
  scale_y_log10(labels = scales::label_comma(accuracy = 0.1), breaks = c(0.1, 1, 10, 100, 1000, 10000)) +
  scale_colour_manual(values = c("black", "grey55", COL[["YQFE"]], COL[["YQFEN"]], COL[["YQ"]], COL[["NYS"]], COL[["SQC"]], COL[["MINGRM"]]), drop = FALSE) +
  labs(x = "Training individuals n (p = 3,000)", y = "Time (s)", tag = "A") + th +
  theme(legend.position = "bottom", legend.text = element_text(size = 7.5)) + guides(colour = guide_legend(ncol = 2))
s5b <- simp[shrink == "REML" & Q %in% c(1000, 3000), .(r = mean(r), time = mean(time)), by = .(method, Q)]
fl <- sim[method == "FULL", .(r = mean(r), time = mean(time))]
p5b <- ggplot(s5b, aes(time, r, colour = method, shape = factor(Q))) + geom_point(size = 2.6) +
  annotate("point", x = fl$time, y = fl$r, size = 3, shape = 8) +
  annotate("text", x = fl$time, y = fl$r, label = "Full GBLUP", hjust = 1.2, vjust = 0.4, size = 3) +
  scale_x_log10() + scale_colour_manual(values = COL, labels = LAB) +
  scale_shape_manual(values = c(`1000` = 16, `3000` = 17), labels = c("Q = 1,000", "Q = 3,000")) +
  labs(x = "Time per fit (s, log scale)", y = "Test correlation (simulation)", tag = "B") + th + theme(legend.position = "right", legend.text = element_text(size = 7.5)) + guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
scp <- fread(file.path(res, "scale_p_results.csv"))
spl <- melt(scp, id.vars = "p", measure.vars = intersect(c("t_exact", "t_yqfe", "t_yqfen", "t_yq", "t_nys"), names(scp)))
spl[, variable := factor(variable, levels = c("t_exact", "t_yqfe", "t_yqfen", "t_yq", "t_nys"), labels = c("Exact GBLUP in marker space (p \u00D7 p)", "YQ-FE (Q = 1,000)",
  "YQ-FEN (Q = 1,000, k = 500)", "YQ (Q = 1,000)", "Nystr\u00F6m (Q = 1,000)"))]
p5c <- ggplot(spl, aes(p, value, colour = variable)) + geom_line() + geom_point(size = 1.8) +
  scale_x_continuous(breaks = scp$p, labels = scales::comma) + scale_y_log10(labels = scales::label_comma(accuracy = 0.1)) +
  scale_colour_manual(values = c("grey35", COL[["YQFE"]], COL[["YQFEN"]], COL[["YQ"]], COL[["NYS"]])) +
  labs(x = "Markers p (n = 8,000)", y = "Time (s)", tag = "C") + th + theme(legend.position = "bottom", legend.text = element_text(size = 7.5)) + guides(colour = guide_legend(ncol = 1))
save_fig((p5a | p5b) / (p5c | plot_spacer()) + plot_layout(heights = c(1.2, 1)), "Figure5_cost.png", 11.5, 9)

## ---------------- Figure 6: real data, r versus Q -----------------------------
s6 <- sumci(real[method != "FULL" & shrink == "REML"], "r", c("method", "Q"))
fr6 <- sumci(real[method == "FULL"], "r", "scenario")
s6[, method := factor(method, levels = names(LAB))]
p6 <- ggplot(s6, aes(Q, mean, colour = method)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = fr6$lo, ymax = fr6$hi, fill = "grey85", alpha = 0.8) +
  geom_hline(yintercept = fr6$mean, linetype = 2) +
  geom_hline(yintercept = mean(margR$r), linetype = 3, colour = "grey30") +
  geom_line(position = position_dodge(60)) + geom_point(size = 1.8, position = position_dodge(60)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 70, position = position_dodge(60)) +
  scale_colour_manual(values = COL, labels = LAB) + scale_x_continuous(breaks = QS) +
  labs(x = "Number of groups / landmarks Q", y = "Test correlation") + th + guides(colour = guide_legend(nrow = 3))
save_fig(p6, "Figure6_real_r_vs_Q.png", 6.8, 5.2)

## ---------------- Figure 7: real data, % of full and speed-up -----------------
s7 <- sumci(realp[shrink == "REML" & Q %in% c(250, 1000, 3000)], "pct", c("method", "Q"))
s7[, method := factor(method, levels = rev(names(LAB)))]
p7a <- ggplot(s7, aes(mean, method, colour = method, shape = factor(Q))) + geom_vline(xintercept = 100, linetype = 2) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, position = position_dodge(0.6)) +
  geom_point(size = 2.2, position = position_dodge(0.6), fill = "white") +
  scale_shape_manual(values = c(`250` = 21, `1000` = 16, `3000` = 17), labels = c("Q = 250", "Q = 1,000", "Q = 3,000")) +
  scale_colour_manual(values = COL, guide = "none") + scale_y_discrete(labels = LAB) +
  labs(x = "% of the full-GBLUP correlation", y = NULL, tag = "A") + th
s7b <- realp[shrink == "REML" & Q %in% c(1000, 3000), .(speedup = mean(speedup)), by = .(method, Q)]
s7b[, method := factor(method, levels = rev(names(LAB)))]
p7b <- ggplot(s7b, aes(speedup, method, fill = method, alpha = factor(Q))) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_alpha_manual(values = c(`1000` = 1, `3000` = 0.45), labels = c("Q = 1,000", "Q = 3,000")) +
  scale_fill_manual(values = COL, guide = "none") + scale_x_log10() + scale_y_discrete(labels = NULL) +
  labs(x = "Speed-up relative to full GBLUP (log scale)", y = NULL, tag = "B") + th
save_fig(p7a + p7b + plot_layout(widths = c(1.2, 1)), "Figure7_real_pct_speed.png", 10.5, 4.8)

## ---------------- Figure 8: sensitivity ---------------------------------------
if (file.exists(file.path(res, "sens_results.csv"))) {
  se <- fread(file.path(res, "sens_results.csv"))
  s8a <- sumci(se[factor == "bins"], "pct", c("scenario", "B"))
  p8a <- ggplot(s8a, aes(B, mean, shape = scenario)) + geom_line(aes(group = scenario), colour = COL[["SQC"]]) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.08, colour = COL[["SQC"]], position = position_dodge(0.1)) +
    geom_point(size = 2.3, colour = COL[["SQC"]], position = position_dodge(0.1)) +
    scale_x_log10(breaks = c(1, 2, 4, 8, 16, 32)) + scale_shape_manual(values = c(16, 15), labels = SCL) +
    labs(x = "Phenotype bins per genomic cluster B (SQC)", y = "% of full-GBLUP correlation", tag = "A") + th
  s8b <- sumci(se[factor == "PCs"], "pct", c("scenario", "method", "k"))
  p8b <- ggplot(s8b, aes(k, mean, colour = method, shape = scenario)) + geom_line(aes(group = interaction(method, scenario))) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.05, position = position_dodge(0.06)) +
    geom_point(size = 2.3, position = position_dodge(0.06)) + scale_x_log10(breaks = c(5, 10, 30, 60)) +
    scale_colour_manual(values = COL, labels = LAB) + scale_shape_manual(values = c(16, 15), labels = SCL) +
    labs(x = "Principal components k", y = "% of full-GBLUP correlation", tag = "B") + th + guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  save_fig(p8a + p8b, "Figure8_sensitivity.png", 10.5, 4.2)
  st6 <- se[, .(pct = round(mean(pct), 2)), by = .(factor, method, B, k, scenario)]
  fwrite(dcast(st6, factor + method + B + k ~ scenario, value.var = "pct"), file.path(tabd, "Table6_sensitivity.csv"))
}

## ---------------- Figure 9: two-stage equivalence -----------------------------
if (file.exists(file.path(res, "two_stage_scatter.csv"))) {
  eq <- fread(file.path(res, "two_stage_scatter.csv"))
  p9 <- ggplot(eq, aes(one_stage, two_stage, colour = factor(reps))) + geom_abline(linetype = 2) +
    geom_point(size = 1.4, alpha = 0.8) + scale_colour_manual(values = c("#9ECAE1", "#4292C6", "#08519C", "#08306B"), name = "Records per line") +
    labs(x = "GEBV, one-stage GBLUP on plot records (lme4)", y = "GEBV, two-stage GBLUP on line means\n(Stage 2 = function used for YQ)") +
    th + theme(legend.title = element_text())
  save_fig(p9, "Figure9_two_stage_equivalence.png", 5.6, 4.6)
}

## ---------------- Figure 10: Nystr\u00F6m captured genomic variance ---------------
ny <- rbind(cbind(sim[method == "NYS"], data = "Simulation (p = 3,000)"), cbind(real[method == "NYS"], data = "Real data (p = 5,527)"))
ny[, data := factor(data, levels = c("Simulation (p = 3,000)", "Real data (p = 5,527)"))]
nyl <- melt(ny, id.vars = c("data", "Q", "rep", "scenario"), measure.vars = c("cap_train", "cap_test"))
s10 <- sumci(nyl, "value", c("data", "Q", "variable"))
s10[, variable := factor(variable, levels = c("cap_train", "cap_test"), labels = c("training individuals", "testing individuals"))]
p10 <- ggplot(s10, aes(Q, 100 * mean, colour = variable)) + geom_line() + geom_point(size = 2) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = 60) + facet_wrap(~data) +
  scale_x_continuous(breaks = QS) + scale_colour_manual(values = c("#7B3294", "#C2A5CF")) +
  labs(x = "Number of landmarks Q", y = "% of genomic self-relationship\nretained by the rank-Q kernel") + th
save_fig(p10, "Figure10_nystrom_capture.png", 8.5, 3.8)

## ---------------- Tables -------------------------------------------------------
t5 <- rbind(simp[shrink == "REML", .(data = "Sim", pct = mean(pct)), by = .(method, Q)],
            realp[shrink == "REML", .(data = "Real", pct = mean(pct)), by = .(method, Q)])
fwrite(dcast(t5, method ~ data + Q, value.var = "pct"), file.path(tabd, "Table5_pct_full.csv"))
t7 <- simp[shrink == "REML" & Q == 1000, .(pct = mean(pct)), by = .(method, scenario)]
t7 <- dcast(t7, method ~ scenario, value.var = "pct")
t7 <- merge(t7, sim[shrink == "REML" & Q == 1000, .(r_g = mean(r_g)), by = method], by = "method")
fwrite(t7, file.path(tabd, "Table7_by_scenario_Q1000.csv"))
C1 <- sim[, .(r = mean(r), r_g = mean(r_g), slope = mean(slope)), by = .(scenario, method, shrink, Q)]
fwrite(C1, file.path(tabd, "TableC1_sim_full.csv"))
C2 <- realp[, .(r = mean(r), pct = mean(pct), slope = mean(slope), mse_ratio = mean(mse_ratio),
                time = mean(time), groups = mean(groups)), by = .(method, shrink, Q)]
fwrite(C2, file.path(tabd, "TableC2_real_full.csv"))
fwrite(sc, file.path(tabd, "TableC3_scaling.csv")); fwrite(scp, file.path(tabd, "TableC4_scaling_p.csv"))
## added references: marginal regression and exact marker-space GBLUP
fs <- sim[method == "FULL", .(scenario, rep, r_full = r, t_full = time)]; frl <- real[method == "FULL", .(rep, r_full = r, t_full = time)]
xs <- merge(asim[method %in% c("MARG", "EXACT")], fs, by = c("scenario", "rep")); xr <- merge(areal[method %in% c("MARG", "EXACT")], frl, by = "rep")
fwrite(rbind(xs[, .(data = "Sim", r = mean(r), pct = mean(100 * r / r_full), slope = mean(slope), time = mean(time),
                    max_abs_diff_r_vs_full = max(abs(r - r_full))), by = method],
             xr[, .(data = "Real", r = mean(r), pct = mean(100 * r / r_full), slope = mean(slope), time = mean(time),
                    max_abs_diff_r_vs_full = max(abs(r - r_full))), by = method]), file.path(tabd, "TableC5_marginal_exact.csv"))
std <- fread(file.path(res, "standardization_check.csv"))
fwrite(rbind(std, as.list(c(rep = NA, colMeans(std[, -1])))), file.path(tabd, "TableC6_standardization.csv"))
key <- list(
  sim_full = sim[method == "FULL", .(r = mean(r), r_g = mean(r_g), time = mean(time)), by = scenario],
  real_full = real[method == "FULL", .(r = mean(r), sd = sd(r), time = mean(time))],
  sim_pct = simp[shrink == "REML", .(pct = mean(pct), slope = mean(slope), mse = mean(mse_ratio), time = mean(time), speed = mean(speedup)), by = .(method, Q)],
  real_pct = realp[shrink == "REML", .(pct = mean(pct), slope = mean(slope), mse = mean(mse_ratio), time = mean(time), speed = mean(speedup)), by = .(method, Q)],
  shrink_cmp = rbind(simp[, .(data = "sim", pct = mean(pct), slope = mean(slope), mse = mean(mse_ratio)), by = .(method, shrink, Q)],
                     realp[, .(data = "real", pct = mean(pct), slope = mean(slope), mse = mean(mse_ratio)), by = .(method, shrink, Q)]))
capture.output(print(key), file = file.path(tabd, "key_numbers.txt"))

## ---------------- Figure 12 and Table C8: fully efficient YQ (v8) ----------------
if (file.exists(file.path(res, "yq_fe_scatter.csv")) && file.exists(file.path(res, "yq_fe_sim_results.csv"))) {
  sc12 <- fread(file.path(res, "yq_fe_scatter.csv"))
  sc12[, data := factor(data, levels = unique(data))]
  l12 <- melt(sc12, id.vars = c("data", "one_stage", "y"), measure.vars = intersect(c("yq_fe", "yq_fen", "yq"), names(sc12)))
  l12[, variable := factor(variable, levels = c("yq_fe", "yq_fen", "yq"), labels = c("YQ-FE (exact)", "YQ-FEN (k = 500)", "YQ (between groups only)"))]
  p12a <- ggplot(l12, aes(one_stage, value, colour = variable)) + geom_abline(linetype = 2) +
    geom_point(size = 0.7, alpha = 0.6) + facet_wrap(~data, scales = "free") +
    scale_colour_manual(values = c(COL[["YQFE"]], COL[["YQFEN"]], "#56B4E9")) +
    labs(x = "Prediction, one-stage GBLUP", y = "Prediction, two-stage method (Q = 1,000)", tag = "A") + th +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
  fe <- rbind(fread(file.path(res, "yq_fe_sim_results.csv")), fread(file.path(res, "yq_fe_real_results.csv")))
  fe[, data := ifelse(scenario == "real_wheat", "Real data (5 partitions)", "Simulation (15 replicates)")]
  fe[, meth := "YQ-FE"]
  fn <- rbind(fread(file.path(res, "yq_fen_sim_results.csv")), fread(file.path(res, "yq_fen_real_results.csv")))[variant == "k500"]
  fn[, data := ifelse(scenario == "real_wheat", "Real data (5 partitions)", "Simulation (15 replicates)")]; fn[, meth := "YQ-FEN"]
  fb <- rbind(fe[, .(Q, data, meth, max_abs_diff_pred_vs_onestage)], fn[, .(Q, data, meth, max_abs_diff_pred_vs_onestage)])
  p12b <- ggplot(fb, aes(factor(Q), max_abs_diff_pred_vs_onestage, colour = meth, shape = data)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 1), size = 1.6) +
    scale_y_log10(limits = c(1e-14, 1e2), breaks = 10^seq(-14, 0, 2), labels = scales::label_math(10^.x, format = log10)) +
    { ck <- fread(file.path(res, "yq_fe_checks.csv"))[method == "YQ (original, approximate)"]
      ck[, data := ifelse(data == "Real_wheat", "Real data (5 partitions)", "Simulation (15 replicates)")]
      geom_point(data = ck, aes(factor(Q), max_abs_diff_vs_full), shape = 4, size = 2.6, stroke = 1.1, colour = "grey45",
                 inherit.aes = FALSE) } +
    annotate("text", x = 0.55, y = 30, hjust = 0, size = 3, colour = "grey35",
             label = "\u00D7 original YQ (first replicate / partition)") +
    scale_colour_manual(values = c(`YQ-FE` = COL[["YQFE"]], `YQ-FEN` = COL[["YQFEN"]])) + scale_shape_manual(values = c(16, 17)) +
    labs(x = "Number of phenotype-quantile groups Q", y = "Largest |difference| from the\none-stage predictions", tag = "B", colour = NULL, shape = NULL) + th
  save_fig(p12a / p12b + plot_layout(heights = c(1.1, 1)), "Figure12_yq_fe_exactness.png", 9, 8)
  c8 <- fe[, .(fits = .N, max_abs_diff_pred = max(max_abs_diff_pred_vs_onestage), max_abs_diff_delta = max(abs(delta - delta_onestage)),
               max_abs_diff_r = max(abs(r - r_onestage)), time = mean(time), time_stage1 = mean(time_stage1), time_stage2 = mean(time_stage2)), by = .(data, Q)]
  fwrite(c8, file.path(tabd, "TableC8_yq_fe.csv"))
}

## =============== v9 main-text figures (Heredity version) =====================
## Main Figure 1: YQ-FE workflow (A) and information split between/within groups (B)
bx2 <- data.frame(x = c(1, 2.4, 3.85, 3.85, 5.55, 5.55, 5.55, 7.15), y = c(0, 0, 0.62, -0.62, 1.15, 0, -1.15, 0),
  w = c(1.1, 1.1, 1.25, 1.25, 1.7, 1.7, 1.7, 1.1), h = c(1.1, 1.1, 1.0, 1.0, 0.95, 0.95, 0.95, 1.1),
  lab = c("Training data\nn individuals\nmarkers Z\nphenotypes y",
          "Q phenotype-\nquantile groups\n(sorted by y)",
          "Stage 1, between:\ngroup means of Z, y\nand sizes m_q",
          "Stage 1, within:\npooled SSCP\nW, s_W, t_W",
          "YQ-FE (exact)\nbetween + within\nREML, p \u00D7 p eigen",
          "YQ-FEN (Nystr\u00F6m)\nbetween + within\nmoments + PCG",
          "YQ (approximate)\nbetween only\nGBLUP on Q means",
          "Marker effects\npredict new\nindividuals"),
  fill = c("#4D4D4D", "#0072B2", "#56B4E9", "#00315E", COL[["YQFE"]], COL[["YQFEN"]], COL[["YQ"]], "#D55E00"))
seg <- data.frame(x = c(1.57, 4.49, 4.49, 4.49, 4.49, 4.49, 6.42, 6.42, 6.42),
                  xend = c(1.83, 4.68, 4.68, 4.68, 4.68, 4.68, 6.58, 6.58, 6.58),
                  y = c(0, 0.75, 0.62, 0.5, -0.5, -0.75, 1.15, 0, -1.15), yend = c(0, 1.15, 0.05, -1.05, 1.0, -0.12, 0.25, 0, -0.25))
pf1a <- ggplot(bx2) + geom_tile(aes(x, y, fill = fill, width = w, height = h), colour = NA) +
  geom_text(aes(x, y, label = lab), colour = "white", size = 2.5, lineheight = 0.9, fontface = "bold") +
  annotate("segment", x = c(2.96, 2.96), xend = c(3.2, 3.2), y = c(0.15, -0.15), yend = c(0.55, -0.55),
           arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend), arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = 4.1, y = -2.0, size = 2.6, fontface = "italic", colour = "grey25", lineheight = 1.05,
           label = "YQ-FE = one-stage GBLUP for any grouping.  YQ-FEN = YQ-FE solution for a moment estimate of \u03B4 (no O(p\u00B3) decomposition).\nYQ = YQ-FE with W = 0, s_W = 0, t_W = 0 (fastest, approximate).") +
  scale_fill_identity() + coord_cartesian(xlim = c(0.4, 7.75), ylim = c(-2.3, 1.7)) + theme_void() +
  labs(tag = "A") + theme(plot.tag = element_text(face = "bold"))
if (file.exists(file.path(res, "information_split.csv"))) {
  inf <- fread(file.path(res, "information_split.csv"))
  il <- melt(inf, id.vars = c("data", "scenario", "rep", "Q"), measure.vars = c("between_markers", "between_phenotype"))
  si <- il[, .(mean = 100 * mean(value), lo = 100 * min(value), hi = 100 * max(value)), by = .(data, Q, variable)]
  si[, variable := factor(variable, labels = c("markers, tr(Z\u2032Z): between groups", "phenotype, y\u2032y: between groups"))]
  si[, data := factor(data, levels = c("Simulation", "Real data"))]
  pf1b <- ggplot(si, aes(Q, mean, colour = variable, linetype = data)) + geom_line() + geom_point(size = 1.6) +
    scale_y_continuous(limits = c(0, 100)) + scale_x_continuous(breaks = c(250, 1000, 2000, 3000, 4000)) +
    scale_colour_manual(values = c("#56B4E9", "#D55E00")) +
    labs(x = "Number of phenotype-quantile groups Q", y = "% of the total (centred) sum of squares\nkept by the group means (YQ)", tag = "B",
         colour = NULL, linetype = NULL) + th + theme(legend.position = "right", legend.text = element_text(size = 7.5))
  save_fig(pf1a / pf1b + plot_layout(heights = c(1.15, 1)), "Main_Figure1_YQFE_workflow.png", 10, 7)
}
## Main Figure 2: percentage of the full-GBLUP correlation versus Q, simulation and real data
m2 <- rbind(cbind(simp[shrink == "REML"], d = "Simulation (15 replicates)"), cbind(realp[shrink == "REML"], d = "Real wheat data (5 partitions)"), fill = TRUE)
s2m <- sumci(m2, "pct", c("d", "method", "Q"))
s2m[, d := factor(d, levels = c("Simulation (15 replicates)", "Real wheat data (5 partitions)"))]
s2m[, method := factor(method, levels = names(LAB))]
pf2 <- ggplot(s2m, aes(Q, mean, colour = method)) + geom_hline(yintercept = 100, linetype = 2) +
  geom_line(position = position_dodge(70)) + geom_point(size = 1.7, position = position_dodge(70)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 70, position = position_dodge(70)) +
  facet_wrap(~d) + scale_colour_manual(values = COL, labels = LAB) + scale_x_continuous(breaks = QS) +
  labs(x = "Number of groups, landmarks or individuals Q", y = "% of the one-stage GBLUP test correlation") + th +
  guides(colour = guide_legend(nrow = 3))
save_fig(pf2, "Main_Figure2_accuracy_vs_Q.png", 10, 5.2)
## Table C9: YQ-FEN (v10) by data, Q and number of landmark markers k
if (file.exists(file.path(res, "yq_fen_sim_results.csv"))) {
  fnr <- rbind(fread(file.path(res, "yq_fen_sim_results.csv")), fread(file.path(res, "yq_fen_real_results.csv")))
  fnr[, data := ifelse(scenario == "real_wheat", "Real", "Sim")]
  c9 <- fnr[, .(fits = .N, pct = mean(100 * r / r_onestage), min_pct = min(100 * r / r_onestage), slope = mean(slope),
                delta = mean(delta), delta_reml = mean(delta_onestage), max_abs_diff_pred = max(max_abs_diff_pred_vs_onestage),
                iter = mean(iter), time = mean(time), time_stage1 = mean(time_stage1), time_stage2 = mean(time_stage2),
                time_yqfe = mean(time_yqfe), time_onestage_exact = mean(time_onestage_exact),
                max_abs_diff_vs_yqfe = max(max_abs_diff_vs_yqfe)), by = .(data, variant, Q, k)]
  fwrite(c9, file.path(tabd, "TableC9_yq_fen.csv"))
}
cat("figures and tables written\n")

## =============== v12 (revision) main-text figures ===========================
## Figure 4 (revised): computing cost measured in one session (21_timing_benchmark.R)
if (file.exists(file.path(res, "timing_benchmark.csv"))) {
  tb <- fread(file.path(res, "timing_benchmark.csv"))
  ML <- c(FULL = "One-stage GBLUP (n × n)", EXACT = "EXACT (p × p, REML)", YQFE = "YQ-FE", YQFEN = "YQ-FEN (k = 500)",
          DIRECT_CG = "DIRECT-CG (no Nyström)", YQ = "YQ (Q = 1,000)", NYS = "Nyström (Q = 1,000)")
  MC <- c(FULL = "black", EXACT = "grey55", YQFE = COL[["YQFE"]], YQFEN = COL[["YQFEN"]], DIRECT_CG = "#E6AB02", YQ = COL[["YQ"]], NYS = COL[["NYS"]])
  ag <- tb[, .(med = median(total), lo = min(total), hi = max(total), mem = max(peak_mem_MB)), by = .(setting, n, p, method)]
  ag[, method := factor(method, levels = names(ML))]
  An <- ag[setting == "A_n"]
  fp <- lm(log(med) ~ log(n), data = An[method == "FULL"]); ext <- data.frame(n = c(6000, 8000, 12000, 20000)); ext$t <- exp(predict(fp, ext))
  q4a <- ggplot(An, aes(n, med, colour = method)) + geom_line() + geom_point(size = 1.6) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.03) +
    geom_line(data = ext, aes(n, t), inherit.aes = FALSE, linetype = 3) +
    annotate("text", x = 14000, y = max(ext$t) / 3, label = "extrapolated", size = 2.6, colour = "grey30") +
    scale_x_log10(breaks = c(2000, 4000, 8000, 20000), labels = scales::comma) + scale_y_log10(labels = scales::label_comma(accuracy = 0.1)) +
    scale_colour_manual(values = MC, labels = ML, drop = FALSE) + labs(x = "Training individuals n (p = 3,000)", y = "Time (s), median and range of 3 runs", tag = "A") +
    th + theme(legend.position = "none")
  Bp <- ag[setting == "B_p"]
  q4b <- ggplot(Bp, aes(p, med, colour = method)) + geom_line() + geom_point(size = 1.6) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 150) +
    scale_y_log10(labels = scales::label_comma(accuracy = 0.1)) + scale_x_continuous(breaks = unique(Bp$p), labels = scales::comma) +
    scale_colour_manual(values = MC, labels = ML, drop = FALSE) + labs(x = "Markers p (n = 8,000)", y = "Time (s)", tag = "B") + th +
    theme(legend.position = "right", legend.text = element_text(size = 7.5)) + guides(colour = guide_legend(ncol = 1))
  st <- tb[setting %in% c("C_sim", "C_real"), .(group = median(group), accumulate = median(accumulate), estimate = median(estimate), predict = median(predict)), by = .(setting, method)]
  stl <- melt(st, id.vars = c("setting", "method")); stl[, setting := factor(setting, levels = c("C_sim", "C_real"), labels = c("Simulation (n = 6,400, p = 3,000)", "Real data (n = 6,434, p = 5,527)"))]
  stl[, method := factor(method, levels = rev(names(ML)))]
  stl[, variable := factor(variable, levels = c("group", "accumulate", "estimate", "predict"), labels = c("sorting/grouping", "accumulation (Stage 1)", "estimation + solution", "prediction"))]
  q4c <- ggplot(stl, aes(value, method, fill = variable)) + geom_col(width = 0.7) + facet_wrap(~setting, scales = "free_x") +
    scale_y_discrete(labels = ML) + scale_fill_manual(values = c("#BDBDBD", "#6BAED6", "#08519C", "#FDAE6B")) +
    labs(x = "Time (s), median of 3 runs", y = NULL, tag = "C", fill = NULL) + th + theme(legend.position = "bottom")
  save_fig((q4a | q4b) / q4c + plot_layout(heights = c(1, 1.1)), "Main_Figure4_cost_benchmark.png", 11.5, 9)
  fwrite(ag, file.path(tabd, "TableR4b_timing_aggregated.csv"))
}
## Figure 5 (revised): calibration with equal training and validation budgets (20_revision_eval.R)
revf <- file.path(res, c("rev_orig_sim.csv", "rev_orig_real.csv", "rev_new_sim.csv", "rev_lodo.csv"))
if (all(file.exists(revf))) {
  RV <- rbindlist(lapply(revf, fread, colClasses = list(character = "rep")), fill = TRUE)
  RV[, grp := fcase(scenario %in% c("S1_polygenic", "S2_sparse", "S3_structure"), "Simulation S1–S3", scenario %in% c("S4_h2_0.2", "S5_h2_0.8"), "Simulation S4–S5 (h² 0.2, 0.8)",
                    scenario %in% c("S6_families", "S7_new_families"), "Simulation S6–S7 (families)", scenario == "real_wheat", "Real, random partitions", default = "Real, leave one set out")]
  ex <- RV[method == "EXACT", .(scenario, rep, mse_ex = mse)]; RV <- merge(RV, ex, by = c("scenario", "rep"))
  RV[, mse_ratio := mse / mse_ex]
  mm <- c("YQ", "YQ_CAL", "YQ_CVCAL", "NYS", "NYS_CAL", "EXACT_CAL", "YQFEN")
  ML5 <- c(YQ = "YQ (REML)", YQ_CAL = "YQ-CAL (90% + 10%)", YQ_CVCAL = "YQ-CVCAL (cross-fitted)", NYS = "Nyström", NYS_CAL = "Nyström-CAL (90% + 10%)",
           EXACT_CAL = "One-stage-CAL (90% + 10%)", YQFEN = "YQ-FEN")
  s5 <- RV[method %in% mm, { v <- ci(mse_ratio); s <- ci(slope); .(mse = v[1], mlo = v[2], mhi = v[3], slope = s[1], slo = s[2], shi = s[3]) }, by = .(grp, method)]
  s5[, method := factor(method, levels = rev(mm))]
  q5a <- ggplot(s5, aes(slope, method, colour = grp)) + geom_vline(xintercept = 1, linetype = 2) +
    geom_pointrange(aes(xmin = slo, xmax = shi), position = position_dodge(0.6), size = 0.25) + scale_y_discrete(labels = ML5) +
    labs(x = "Calibration slope (Q = 1,000)", y = NULL, tag = "A", colour = NULL) + th + theme(legend.position = "none")
  q5b <- ggplot(s5, aes(mse, method, colour = grp)) + geom_vline(xintercept = 1, linetype = 2) +
    geom_pointrange(aes(xmin = mlo, xmax = mhi), position = position_dodge(0.6), size = 0.25) + scale_y_discrete(labels = NULL) + scale_x_log10() +
    labs(x = "MSE relative to the one-stage GBLUP (log scale)", y = NULL, tag = "B", colour = NULL) + th + theme(legend.position = "right", legend.text = element_text(size = 7.5))
  save_fig(q5a + q5b + plot_layout(widths = c(1, 1.15)), "Main_Figure5_calibration_budget.png", 11.5, 4.6)
}
cat("v12 figures written\n")
