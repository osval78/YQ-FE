###############################################################################
# 16_information_split.R -- how the information of the training data is split
# between the phenotype-quantile groups (used by YQ) and within them (added by
# YQ-FE): proportion of the total centred marker sum of squares, tr(Zc'Zc), and
# of the phenotypic sum of squares that lies between groups, for Q = 250 ...
# 3,000, in every simulated replicate and real-data partition (Figure 1B).
# Output: outputs/results/information_split.csv
###############################################################################
source("00_functions.R")
QS <- c(250, 500, 1000, 2000, 3000, 4000)
split_info <- function(Z, y) {
  zc <- colMeans(Z); tot_z <- sum(sweep(Z, 2, zc)^2); tot_y <- sum((y - mean(y))^2)
  do.call(rbind, lapply(QS, function(Q) {
    lab <- make_groups("YQ", y, Q); S <- summarise_groups(Z, y, lab)
    bz <- sum(S$m * rowSums(sweep(S$Zbar, 2, zc)^2)); by <- sum(S$m * (S$ybar - mean(y))^2)
    data.frame(Q = Q, between_markers = bz / tot_z, between_phenotype = by / tot_y)
  }))
}
rows <- list()
for (si in seq_along(SIM$SCEN)) for (rep in seq_len(SIM$NREP) - 1) {
  d <- make_sim_rep(si, rep); x <- split_info(d$Z[d$tr, ], d$y[d$tr])
  rows[[length(rows) + 1]] <- cbind(data = "Simulation", scenario = SIM$SCEN[si], rep = rep, x); rm(d); gc(verbose = FALSE)
}
D <- readRDS(file.path("..", "data", "processed", "real_pooled.rds"))
for (rep in 0:4) {
  pt <- readRDS(file.path("..", "outputs", "results", sprintf("real_partition_%d.rds", rep)))
  rows[[length(rows) + 1]] <- cbind(data = "Real data", scenario = "real_wheat", rep = rep, split_info(D$Z[pt$tr, ], D$y[pt$tr]))
}
out <- do.call(rbind, rows); print(aggregate(cbind(between_markers, between_phenotype) ~ data + Q, out, mean))
write.csv(out, file.path("..", "outputs", "results", "information_split.csv"), row.names = FALSE)
