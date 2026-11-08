###############################################################################
# 05_sensitivity.R -- design choices of the cluster-based rules at Q = 1,000:
# phenotype bins per genomic cluster (B) in SQC and number of principal
# components (k) in SQC and GEN. Same simulated data as 02_simulation.R
# (scenarios S1 and S3, 5 replicates). Revision: the principal components and
# the k-means start use the same seed as 02_simulation.R (d$seed), so that the
# configurations shared with the main analysis (SQC with B = 4 and k = 30; GEN
# with k = 30) reproduce exactly the values of Table 8.
###############################################################################
source("00_functions.R")
sim <- read.csv(file.path("..", "outputs", "results", "sim_results.csv"))
out <- file.path("..", "outputs", "results", "sens_results.csv")
old <- if (file.exists(out)) read.csv(out) else NULL
rows <- if (is.null(old)) list() else list(old); Q <- 1000
done <- if (is.null(old)) character(0) else unique(paste(old$scenario, old$rep))
for (si in c(1, 3)) for (rep in seq_len(SIM$NREP) - 1) {
  if (paste(SIM$SCEN[si], rep) %in% done) next
  sc <- SIM$SCEN[si]; d <- make_sim_rep(si, rep); sd0 <- d$seed
  Ztr <- d$Z[d$tr, ]; ytr <- d$y[d$tr]; Zte <- d$Z[d$te, ]; yte <- d$y[d$te]; rm(d)
  rf <- sim$r[sim$scenario == sc & sim$rep == rep & sim$method == "FULL"]
  for (k in c(5, 10, 30, 60)) {
    P <- top_pcs(Ztr, k, sd0)
    for (meth in c("GEN", "SQC")) {
      S <- summarise_groups(Ztr, ytr, make_groups(meth, ytr, Q, P, B = 4, seed = sd0))
      pr <- predict_wgblup(fit_wgblup(S$Zbar, S$ybar, S$m), Zte)
      rows[[length(rows) + 1]] <- data.frame(scenario = sc, rep = rep, factor = "PCs",
        method = meth, B = ifelse(meth == "SQC", 4, NA), k = k, pct = 100 * cor(pr, yte) / rf)
    }
  }
  P <- top_pcs(Ztr, 30, sd0)
  for (B in c(1, 2, 4, 8, 16, 32)) {
    S <- summarise_groups(Ztr, ytr, make_groups("SQC", ytr, Q, P, B = B, seed = sd0))
    pr <- predict_wgblup(fit_wgblup(S$Zbar, S$ybar, S$m), Zte)
    rows[[length(rows) + 1]] <- data.frame(scenario = sc, rep = rep, factor = "bins",
      method = "SQC", B = B, k = 30, pct = 100 * cor(pr, yte) / rf)
  }
  write.csv(do.call(rbind, rows), out, row.names = FALSE)
  message(sc, " rep ", rep + 1); gc(verbose = FALSE)
}
