# TABLES v5 (continuity study). Reads <root>/continuity_<MODE>/rep.csv and sens.csv; writes tables/*.csv. Base R only.
# AMENDMENT 2026-10-17 (after the confirmatory run, before any table was read): replications with an undefined value (benchmark equal to 1 in a start period) are dropped with na.rm; see docs/PROTOCOL_continuity_v5.md section 6.
# Measures: Pi (random-exit benchmark), PiC (continuity-adjusted, 10 classes of uninterrupted duration; primary), PiC5 (5 coarse classes).
# depth 0 = depth-weighted summary over all list positions; depth 5, 10, 20, 50 = fixed top-k cohorts.
MODE <- Sys.getenv("SIM_MODE", unset = "confirm"); ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-persistence")
OUT <- file.path(ROOT, paste0("continuity_", MODE)); TAB <- file.path(OUT, "tables"); dir.create(TAB, showWarnings = FALSE)
D <- read.csv(file.path(OUT, "rep.csv")); L5 <- D[D$lag == 5, ]; key <- function(d) paste(d$beta, d$delta, d$lambda, d$gamma)
FIX <- c(10, 20, 50)                                                                     # fixed depths with hypothesis tests; k = 5 is descriptive (cohort of five members, noisy)
M1 <- 0.02; M2 <- 0.04; M3 <- 0.05; M4 <- 0.03; M5 <- 0.08; M6W <- 0.03; M6F <- 0.05; M8 <- 0.03; M9 <- 0.02; M10 <- 0.05; M11 <- 0.06   # margins fixed in docs/PROTOCOL_continuity_v5.md
cells <- unique(L5[c("beta", "delta", "lambda", "gamma")]); rownames(cells) <- NULL
cellstat <- function(v, dp, d = L5) { d <- d[d$depth == dp, ]; do.call(rbind, lapply(seq_len(nrow(cells)), function(i) { x <- d[[v]][key(d) == key(cells[i, ])]
  x <- x[!is.na(x)]; data.frame(cells[i, ], depth = dp, mean = mean(x), se = sd(x) / sqrt(length(x)), n = length(x)) })) }
pone <- function(m, se, n, margin, dir) pt((if (dir == "gt") (m - margin) else (margin - m)) / se, n - 1)
iut_gt <- function(cs, margin) max(1 - pone(cs$mean, cs$se, cs$n, margin, "gt"))
iut_eq <- function(cs, margin) max(pmax(1 - pone(cs$mean, cs$se, cs$n, -margin, "gt"), 1 - pone(cs$mean, cs$se, cs$n, margin, "lt")))
allcells <- do.call(rbind, lapply(c(0, 5, 10, 20, 50), function(dp) { a <- cellstat("Pi", dp); b <- cellstat("PiC", dp); c5 <- cellstat("PiC5", dp)
  merge(merge(a, b, by = c("beta", "delta", "lambda", "gamma", "depth"), suffixes = c("_Pi", "_PiC")), setNames(c5, c(names(c5)[1:5], "mean_PiC5", "se_PiC5", "n_PiC5")), by = c("beta", "delta", "lambda", "gamma", "depth")) }))
write.csv(allcells, file.path(TAB, "cells_lag5.csv"), row.names = FALSE)
sel0 <- function(d) d$beta == 0 & d$delta == 0 & d$gamma == 0; h2f <- function(d) d$beta == 0 & (d$delta > 0 | d$gamma > 0)
cs <- function(v, dps, sel) do.call(rbind, lapply(dps, function(dp) { a <- cellstat(v, dp); a[sel(a), ] }))
H <- list(); add <- function(id, cells, est, worst, crit, p) H[[length(H) + 1]] <<- data.frame(id = id, cells = cells, est = est, worst = worst, criterion = crit, p = p)
# H1 Propositions 1 and 2 under random exit, depth-weighted
a <- cs("Pi", 0, sel0); b <- cs("PiC", 0, sel0); add("H1", nrow(a), max(abs(c(a$mean, b$mean))), max(abs(c(a$mean, b$mean)) + 1.645 * c(a$se, b$se)), "90% intervals of both measures within +/-0.02", max(iut_eq(a, M1), iut_eq(b, M1)))
# H2 the same at fixed depths
a <- cs("Pi", FIX, sel0); b <- cs("PiC", FIX, sel0); add("H2", nrow(a) + nrow(b), max(abs(c(a$mean, b$mean))), max(abs(c(a$mean, b$mean)) + 1.645 * c(a$se, b$se)), "90% intervals within +/-0.04 at k = 10, 20, 50", max(iut_eq(a, M2), iut_eq(b, M2)))
# H3 false persistence in Pi (beta = 0, delta > 0 or gamma > 0) at depth 0, 10, 20, 50
a <- cs("Pi", c(0, FIX), h2f); add("H3", nrow(a), min(a$mean), min(a$mean - 1.645 * a$se), "Pi above 0.05 in every cell and depth", iut_gt(a, M3))
# H4 adjustment, depth-weighted
a <- cs("PiC", 0, h2f); add("H4", nrow(a), max(abs(a$mean)), max(abs(a$mean) + 1.645 * a$se), "90% interval of PiC within +/-0.03", iut_eq(a, M4))
# H5 adjustment at fixed depths
a <- cs("PiC", FIX, h2f); add("H5", nrow(a), max(abs(a$mean)), max(abs(a$mean) + 1.645 * a$se), "90% interval of PiC within +/-0.08 at k = 10, 20, 50", iut_eq(a, M5))
# H6 detection (beta > 0)
a <- cs("PiC", 0, function(d) d$beta > 0); b <- cs("PiC", FIX, function(d) d$beta > 0)
add("H6", nrow(a) + nrow(b), min(c(a$mean, b$mean)), min(c(a$mean - 1.645 * a$se, b$mean - 1.645 * b$se)), "PiC above 0.03 (depth-weighted) and 0.05 (k = 10, 20, 50) in every cell", max(iut_gt(a, M6W), iut_gt(b, M6F)))
# H7 dose response: PiC rises from beta 0 to 1 and from 1 to 2 within every (delta, lambda, gamma), at depth 0, 10, 20, 50
combos <- unique(L5[c("delta", "lambda", "gamma")]); pv <- c(); es <- c()
for (dp in c(0, FIX)) for (i in seq_len(nrow(combos))) for (bb in 0:1) { f <- function(be) L5$PiC[L5$depth == dp & L5$beta == be & L5$delta == combos$delta[i] & L5$lambda == combos$lambda[i] & L5$gamma == combos$gamma[i]]
  y <- f(c(0, 1)[bb + 1]); x <- f(c(1, 2)[bb + 1]); y <- y[!is.na(y)]; x <- x[!is.na(x)]; d <- mean(x) - mean(y); se <- sqrt(var(x) / length(x) + var(y) / length(y)); es <- c(es, d); pv <- c(pv, 1 - pt(d / se, length(x) - 1)) }
add("H7", length(es), min(es), NA, "PiC rises with beta in every step and depth", max(pv))
# H8 granularity: coarse classes leave more false persistence than fine classes at k = 20 (beta = 0, delta > 0, lambda > 0)
w <- L5[L5$depth == 20 & L5$beta == 0 & L5$delta > 0 & L5$lambda > 0, ]; w$d <- w$PiC5 - w$PiC
cc <- do.call(rbind, lapply(split(w, key(w)), function(z) { d <- z$d[!is.na(z$d)]; data.frame(mean = mean(d), se = sd(d) / sqrt(length(d)), n = length(d)) }))
add("H8", nrow(cc), min(cc$mean), min(cc$mean - 1.645 * cc$se), "PiC5 - PiC above 0.03 at k = 20 in every cell", iut_gt(cc, M8))
# H9 depth profile: rank information is larger in smaller cohorts, PiC(10) - PiC(50) above 0.02, mean within beta = 1 and beta = 2
w <- merge(L5[L5$depth == 10 & L5$beta > 0, c("beta", "delta", "lambda", "gamma", "rep", "PiC")], L5[L5$depth == 50 & L5$beta > 0, c("beta", "delta", "lambda", "gamma", "rep", "PiC")], by = c("beta", "delta", "lambda", "gamma", "rep"), suffixes = c("10", "50")); w$d <- w$PiC10 - w$PiC50
cc <- do.call(rbind, lapply(split(w, w$beta), function(z) { d <- z$d[!is.na(z$d)]; data.frame(mean = mean(d), se = sd(d) / sqrt(length(d)), n = length(d)) }))
add("H9", nrow(cc), min(cc$mean), min(cc$mean - 1.645 * cc$se), "PiC(10) - PiC(50) above 0.02, mean at each beta level", iut_gt(cc, M9))
# H10 lag: Pi(5) - Pi(1) above 0.05 in the H3 cells, depth-weighted
w <- merge(D[D$lag == 5 & D$depth == 0, c("beta", "delta", "lambda", "gamma", "rep", "Pi")], D[D$lag == 1 & D$depth == 0, c("beta", "delta", "lambda", "gamma", "rep", "Pi")], by = c("beta", "delta", "lambda", "gamma", "rep"), suffixes = c("5", "1"))
w$d <- w$Pi5 - w$Pi1; w <- w[h2f(w), ]; cc <- do.call(rbind, lapply(split(w, key(w)), function(z) { d <- z$d[!is.na(z$d)]; data.frame(mean = mean(d), se = sd(d) / sqrt(length(d)), n = length(d)) }))
add("H10", nrow(cc), min(cc$mean), min(cc$mean - 1.645 * cc$se), "Pi(5) - Pi(1) above 0.05 in every cell", iut_gt(cc, M10))
# H11 cost of adjustment: beta > 0 with delta = lambda = gamma = 0, |Pi - PiC| within 0.06, depth-weighted
a <- L5[L5$depth == 0 & L5$beta > 0 & L5$delta == 0 & L5$lambda == 0 & L5$gamma == 0, ]; a$d <- a$Pi - a$PiC
cc <- do.call(rbind, lapply(split(a, a$beta), function(z) { d <- z$d[!is.na(z$d)]; data.frame(mean = mean(d), se = sd(d) / sqrt(length(d)), n = length(d)) }))
add("H11", nrow(cc), max(abs(cc$mean)), max(abs(cc$mean) + 1.645 * cc$se), "90% interval of Pi - PiC within +/-0.06", iut_eq(cc, M11))
Ht <- do.call(rbind, H); Ht$p_holm <- p.adjust(Ht$p, "holm"); Ht$decision <- ifelse(Ht$p_holm < 0.05, ifelse(grepl("within", Ht$criterion), "Equivalent", "Supported"), "Not supported")
write.csv(Ht, file.path(TAB, "hypotheses.csv"), row.names = FALSE); print(Ht, digits = 3)
# canonical scenarios, all lags and depths
CAN <- data.frame(label = c("Random exit", "Duration-dependent exit", "Frailty only", "Continuity link and duration-dependent exit", "Rank effect", "Rank effect, continuity channels", "Strong rank effect", "Strong rank effect, continuity channels"),
  beta = c(0, 0, 0, 0, 1, 1, 2, 2), delta = c(0, 1, 0, 1, 0, 1, 0, 1), lambda = c(0, 0, 0, 1, 0, 1, 0, 1), gamma = c(0, 0, 1, 0, 0, 0, 0, 0))
mse <- function(x) { x <- x[!is.na(x)]; c(mean(x), sd(x) / sqrt(length(x))) }
can <- do.call(rbind, lapply(seq_len(nrow(CAN)), function(i) { d <- D[key(D) == key(CAN[i, ]), ]
  do.call(rbind, lapply(split(d, list(d$lag, d$depth)), function(z) { p <- mse(z$Pi); q <- mse(z$PiC); r <- mse(z$PiC5)
    data.frame(CAN[i, ], lag = z$lag[1], depth = z$depth[1], Pi = p[1], Pi_se = p[2], PiC = q[1], PiC_se = q[2], PiC5 = r[1], PiC5_se = r[2], rho = mean(z$rho), S = mean(z$S), SC = mean(z$SC)) })) }))
rownames(can) <- NULL; write.csv(can, file.path(TAB, "canonical_by_lag.csv"), row.names = FALSE)
if (file.exists(file.path(OUT, "sens.csv"))) { S <- read.csv(file.path(OUT, "sens.csv")); S <- S[S$lag == 5, ]
  sn <- aggregate(cbind(Pi, PiC, PiC5) ~ N + EXIT + beta + delta + lambda + gamma + depth, S, mean, na.action = na.pass); write.csv(sn, file.path(TAB, "sensitivity_lag5.csv"), row.names = FALSE) }
cat("Tables written to", TAB, "\n")
