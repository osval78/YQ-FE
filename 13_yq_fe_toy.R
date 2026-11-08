###############################################################################
# 13_yq_fe_toy.R -- toy example of the exact two-stage YQ (YQ-FE): data of
# Table 2 (n = 12, p = 4, Q = 3 phenotype-quantile groups of 4 individuals),
# delta fixed at 1 so that every number can be checked by hand.
# Output: outputs/tables/yq_fe_toy.txt
###############################################################################
source("00_functions.R")
sink(file.path("..", "outputs", "tables", "yq_fe_toy.txt"), split = TRUE)
Z <- rbind(c(-0.01, 1.05, 0.74, 0.72), c(1.62, -1.21, -0.63, -1.32), c(-0.11, 1.00, -0.02, 0.50),
           c(-1.91, 0.15, -0.91, 1.78), c(0.89, 0.95, -0.06, 0.61), c(0.66, -0.34, -0.50, -0.11),
           c(-0.61, -0.59, -0.28, -0.73), c(0.77, -1.60, 0.82, -0.63), c(-0.55, -1.35, -0.14, -0.25),
           c(0.19, -0.53, 0.09, 1.82), c(0.41, -0.57, 0.95, -0.13), c(0.59, 0.61, -0.39, -1.93))
y <- c(-0.62, 2.37, -1.01, -1.82, 0.99, 0.16, 0.68, 3.27, 1.12, 1.46, 1.97, 1.81)
n <- 12; p <- 4; delta <- 1
lab <- make_groups("YQ", y, 3); cat("groups:", lab, "\n")
# one-stage: mixed-model equations [1'1 1'Z; Z'1 Z'Z + p delta I] [mu; b] = [1'y; Z'y]
mme <- function(XtX, Xty) solve(XtX, Xty)
X1 <- cbind(1, Z); L <- diag(c(0, rep(p * delta, p)))
th1 <- mme(crossprod(X1) + L, crossprod(X1, y)); cat("\none-stage (mu, b):", round(th1, 5), "\n")
# exact two-stage via Stage 1 summaries
S1 <- yq_fe_stage1(Z, y, lab)
cat("\nStage 1: group sizes", S1$m, "\n group means of y:", round(S1$ybar, 4), "\n group means of markers:\n"); print(round(S1$Zbar, 4))
cat("\n within-group SSCP W:\n"); print(round(S1$W, 4)); cat(" sW:", round(S1$sW, 4), "  tW:", round(S1$tW, 4), "\n")
m <- S1$m
XtX2 <- rbind(c(sum(m), drop(crossprod(m, S1$Zbar))), cbind(drop(crossprod(S1$Zbar, m)), crossprod(S1$Zbar * sqrt(m)) + S1$W))
Xty2 <- c(sum(m * S1$ybar), drop(crossprod(S1$Zbar, m * S1$ybar)) + S1$sW)
th2 <- mme(XtX2 + L, Xty2); cat("\nexact two-stage YQ-FE (mu, b):", round(th2, 5), "\n")
cat("max |difference| =", signif(max(abs(th1 - th2)), 3), "\n")
# original YQ: between-group part only
XtX0 <- rbind(c(sum(m), drop(crossprod(m, S1$Zbar))), cbind(drop(crossprod(S1$Zbar, m)), crossprod(S1$Zbar * sqrt(m))))
th0 <- mme(XtX0 + L, Xty2 - c(0, S1$sW)); cat("original YQ (mu, b):", round(th0, 5), "\n")
# data form with Helmert contrasts
C <- as.matrix(within_contrasts(lab)); lab2 <- as.integer(factor(lab))
cat("\nHelmert contrasts of group 1 (rows = members in data order):\n"); print(round(C[lab2 == 1, 1:3], 4))
u <- sqrt(m) * S1$ybar; w <- drop(crossprod(C, y)); Xu <- sqrt(m) * S1$Zbar; Xw <- crossprod(C, Z)
X2 <- rbind(cbind(sqrt(m), Xu), cbind(0, Xw)); y2 <- c(u, w)
th3 <- mme(crossprod(X2) + L, crossprod(X2, y2)); cat("data form (means + contrasts) (mu, b):", round(th3, 5), "\n")
cat("transformed data: u =", round(u, 4), "\n w =", round(w, 4), "\n")
cat("check ||y||^2 =", round(sum(y^2), 4), " ||u||^2 + ||w||^2 =", round(sum(u^2) + sum(w^2), 4), "\n")
sink()
