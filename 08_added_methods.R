###############################################################################
# 08_added_methods.R -- methods added in the revision, evaluated on exactly the
# same simulated data sets and real-data partitions as 02_simulation.R and
# 03_real_data.R (the data and partitions are regenerated from the same seeds):
#   MINGRM Min_GRM training-set selection: the Q training individuals whose
#          minimum genomic relationship with the testing set is largest
#          (Fernandez-Gonzalez et al., 2024); genotypes only, no phenotypes
#   EXACT exact GBLUP in marker space (p x p cross-products; identical to FULL)
#   MARG  one-pass marginal-regression predictor, b ~ Z'(y - ybar)
# Also (real data): effect of standardizing GY with training lines only.
# Resumable at the level of replicate x method (only missing methods are run). Output: outputs/results/added_sim_results.csv, added_real_results.csv,
# standardization_check.csv
###############################################################################
source("00_functions.R")
resd <- file.path("..", "outputs", "results")
## ---- simulation ----
out <- file.path(resd, "added_sim_results.csv")
rows <- if (file.exists(out)) read.csv(out) else NULL
METH <- c("EXACT", "MARG", "MINGRM")
todo <- function(rows, sc, rep) {
  if (is.null(rows)) return(METH)
  setdiff(METH, unique(rows$method[rows$scenario == sc & rows$rep == rep]))
}
for (si in seq_along(SIM$SCEN)) for (rep in seq_len(SIM$NREP) - 1) {
  sc <- SIM$SCEN[si]; mm <- todo(rows, sc, rep); if (!length(mm)) next
  message(sprintf("%s rep %d: %s", sc, rep + 1, paste(mm, collapse = ", "))); t0 <- tnow()
  d <- make_sim_rep(si, rep)
  res <- evaluate_added(d$Z[d$tr, ], d$y[d$tr], d$Z[d$te, ], d$y[d$te], d$g[d$te], SIM$QS, methods = mm)
  res$scenario <- sc; res$rep <- rep; res$n_train <- length(d$tr); res$p <- SIM$P
  rows <- rbind(rows, res); write.csv(rows, out, row.names = FALSE)
  rm(d); gc(verbose = FALSE); message(sprintf("  %.0f s", tnow() - t0))
}
## ---- real data ----
if (Sys.getenv("SKIP_REAL") == "1") quit(save = "no")
D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds"))
Z <- D$Z; y <- D$y; QS <- c(250, 1000, 2000, 3000)
out <- file.path(resd, "added_real_results.csv"); outs <- file.path(resd, "standardization_check.csv")
rows <- if (file.exists(out)) read.csv(out) else NULL
srows <- if (file.exists(outs)) read.csv(outs) else NULL
for (rep in 0:4) {
  mm <- todo(rows, "real_wheat", rep); if (!length(mm)) next
  pt <- readRDS(file.path(resd, sprintf("real_partition_%d.rds", rep))); tr <- pt$tr; te <- pt$te
  message(sprintf("real partition %d: %s", rep + 1, paste(mm, collapse = ", "))); t0 <- tnow()
  res <- evaluate_added(Z[tr, ], y[tr], Z[te, ], y[te], NULL, QS, methods = mm)
  res$scenario <- "real_wheat"; res$rep <- rep; res$n_train <- length(tr); res$p <- ncol(Z)
  rows <- rbind(rows, res); write.csv(rows, out, row.names = FALSE)
  message(sprintf("  %.0f s", tnow() - t0)); gc(verbose = FALSE)
}
sdone <- if (is.null(srows)) integer(0) else unique(srows$rep)
for (rep in 0:4) {
  if (rep %in% sdone) next
  pt <- readRDS(file.path(resd, sprintf("real_partition_%d.rds", rep))); tr <- pt$tr; te <- pt$te
  message(sprintf("standardization check, partition %d", rep + 1)); t0 <- tnow()
  # GY standardized within data set with training lines only (no test information)
  y2 <- y
  for (e in unique(D$env)) { ie <- which(D$env == e); it <- intersect(ie, tr)
    y2[ie] <- (y[ie] - mean(y[it])) / sd(y[it]) }
  f1 <- fit_marker_gblup(Z[tr, ], y[tr]); f2 <- fit_marker_gblup(Z[tr, ], y2[tr])
  S1 <- summarise_groups(Z[tr, ], y[tr], make_groups("YQ", y[tr], 1000))
  S2 <- summarise_groups(Z[tr, ], y2[tr], make_groups("YQ", y2[tr], 1000))
  q1 <- predict_wgblup(fit_wgblup(S1$Zbar, S1$ybar, S1$m), Z[te, ])
  q2 <- predict_wgblup(fit_wgblup(S2$Zbar, S2$ybar, S2$m), Z[te, ])
  srows <- rbind(srows, data.frame(rep = rep,
    r_full_allstd = cor(predict_marker_gblup(f1, Z[te, ]), y[te]),
    r_full_trainstd = cor(predict_marker_gblup(f2, Z[te, ]), y2[te]),
    r_yq1000_allstd = cor(q1, y[te]), r_yq1000_trainstd = cor(q2, y2[te])))
  write.csv(srows, outs, row.names = FALSE)
  message(sprintf("  %.0f s", tnow() - t0)); gc(verbose = FALSE)
}
