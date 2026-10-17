# SIMULATION v5: does rank add information beyond continuity (uninterrupted duration)? Continuity-adjusted excess persistence, fixed top-k depths. (Internal variable "ten" is the uninterrupted duration.)
# Run:  SIM_MODE=pilot|explore|confirm  SIM_OUTPUT_ROOT=<folder>  SIM_CORES=<n>  Rscript persistence_continuity_sim.R
# Model. A ranked list of N members. Member i has a permanent component z_i ~ N(0,1) and a transitory AR(1) component w_it (coefficient PHI_W, stationary SD SD_W).
# Score x_it = z_i + lambda * log(duration_it) + w_it. Rank 1 = highest score. r_it = (rank-1)/(N-1) in [0,1] (0 best).
# Each transition exactly X = round(EXIT*N) members leave. They are drawn without replacement with weights exp(beta*r_it + gamma*u_i - delta*log(duration_it)),
# where u_i ~ N(0,1) is an unobserved member frailty independent of z_i. Leavers are replaced by new members (duration 1, fresh z, u, w). A leaver never returns.
#   beta   = effect of current rank on exit (0 = rank carries no information about exit)
#   delta  = duration dependence of exit (exit weight falls with log duration)
#   lambda = link between uninterrupted duration and score
#   gamma  = unobserved persistent frailty (nuisance heterogeneity)
# Measures. Pi = depth-weighted excess persistence against the random-exit benchmark S(h) = (1 - X/N)^h.
#   PiC = same with the continuity-adjusted benchmark S_C: exit probabilities estimated per period within classes of uninterrupted duration (1, 2, 3-4, 5-8, 9+ periods on the list).
SEED <- c(pilot = 1101L, explore = 7931L, confirm = 9523L); MODE <- Sys.getenv("SIM_MODE", unset = "pilot")
ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-persistence"); CORES <- as.integer(Sys.getenv("SIM_CORES", unset = "2"))
REPS <- c(pilot = 20L, explore = 150L, confirm = 500L)[[MODE]]; OUT <- file.path(ROOT, paste0("continuity_", MODE)); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
NN <- 200L; EXIT0 <- 0.10; TN <- 46L; STARTS <- c(21L, 24L, 27L, 30L, 33L, 36L); LAGS <- c(1L, 2L, 3L, 5L, 10L)
PHI_W <- 0.7; SD_W <- 0.7; CLS <- c(0, 1, 2, 4, 8, Inf); CLF <- c(0, 1, 2, 3, 4, 6, 8, 11, 16, 24, Inf);   # coarse (5) and fine (10) classes of uninterrupted duration
 KS <- c(5L, 10L, 20L, 50L)   # fixed top-k depths; k = 0 codes the depth-weighted summary over all list positions
GRID <- expand.grid(beta = c(0, 1, 2), delta = c(0, 0.5, 1), lambda = c(0, 0.5, 1), gamma = c(0, 1))
if (MODE == "pilot") GRID <- GRID[c(1, 27, 28, 54), ]
SENS <- rbind(expand.grid(N = c(50L, 800L), EXIT = c(0.05, 0.20), beta = 0, delta = 0, lambda = 0, gamma = 0), expand.grid(N = c(50L, 800L), EXIT = c(0.05, 0.20), beta = 0, delta = 1, lambda = 1, gamma = 0),
              expand.grid(N = c(50L, 800L), EXIT = c(0.05, 0.20), beta = 1, delta = 0.5, lambda = 0.5, gamma = 1), expand.grid(N = c(50L, 800L), EXIT = c(0.05, 0.20), beta = 2, delta = 1, lambda = 1, gamma = 0))
SENS_REPS <- c(pilot = 5L, explore = 60L, confirm = 200L)[[MODE]]
library(parallel)

one_run <- function(beta, delta, lambda, gamma, N = NN, EXIT = EXIT0) {
  X <- as.integer(round(EXIT * N))
  ten <- matrix(0L, TN, N); rk <- matrix(0L, TN, N); id <- matrix(0L, TN, N)      # per period and slot: uninterrupted duration, rank, member id
  z <- rnorm(N); u <- rnorm(N); w <- rnorm(N, 0, SD_W); tn <- rep(1L, N); ids <- seq_len(N); nid <- N + 1L; a <- sqrt(1 - PHI_W^2) * SD_W
  for (s in seq_len(TN)) {
    x <- z + lambda * log(tn) + w; r <- rank(-x, ties.method = "first"); ten[s, ] <- tn; rk[s, ] <- r; id[s, ] <- ids
    if (s == TN) break
    wt <- exp(beta * (r - 1) / (N - 1) + gamma * u - delta * log(tn))
    leave <- order(runif(N)^(1 / wt), decreasing = TRUE)[seq_len(X)]                # weighted draw without replacement (Efraimidis-Spirakis)
    w <- PHI_W * w + a * rnorm(N); tn <- tn + 1L
    z[leave] <- rnorm(X); u[leave] <- rnorm(X); w[leave] <- rnorm(X, 0, SD_W); tn[leave] <- 1L; ids[leave] <- nid:(nid + X - 1L); nid <- nid + X
  }
  Emat <- function(brk) { nb <- length(brk) - 1L; cl <- matrix(pmin(nb, pmax(1L, as.integer(cut(ten, brk, labels = FALSE)))), TN, N); E <- matrix(0, TN, nb)   # E[s, b]: exit rate at transition s-1 -> s among members of class b at s-1
    for (s in 2:TN) { gone <- !(id[s - 1, ] %in% id[s, ]); for (b in seq_len(nb)) { m <- cl[s - 1, ] == b; E[s, b] <- if (any(m)) mean(gone[m]) else 0 } }; E }
  bench <- function(E, brk, t, h, tt) { Si <- rep(1, N); lo <- brk[-length(brk)] + 1; for (q in 0:(h - 1)) Si <- Si * (1 - E[t + 1L + q, pmin(length(lo), findInterval(tt + q, lo))]); cumsum(Si) / seq_len(N) }
  E1 <- Emat(CLS); E2 <- Emat(CLF)
  nd <- length(KS) + 1L; acc <- matrix(0, length(LAGS) * nd, 7L)
  for (t in STARTS) { o <- order(rk[t, ]); mem <- id[t, o]; tt <- ten[t, o]; k <- seq_len(N)
    for (j in seq_along(LAGS)) { h <- LAGS[j]; ret <- as.numeric(mem %in% id[t + h, ]); rho <- cumsum(ret) / k
      S <- (1 - X / N)^h; pik <- (rho - S) / (1 - S)
      SC <- bench(E1, CLS, t, h, tt); SF <- bench(E2, CLF, t, h, tt)
      wt <- function(Sx) { ok <- (1 - Sx) > 1e-8; if (any(ok)) sum(k[ok] * (rho[ok] - Sx[ok]) / (1 - Sx[ok])) / sum(k[ok]) else NA }
      pic <- (rho - SC) / (1 - SC); pif <- (rho - SF) / (1 - SF)
      acc[(j - 1L) * nd + seq_len(nd), ] <- acc[(j - 1L) * nd + seq_len(nd), ] + rbind(c(sum(k * pik) / sum(k), wt(SC), wt(SF), rho[N], S, SC[N], SF[N]), cbind(pik[KS], pic[KS], pif[KS], rho[KS], S, SC[KS], SF[KS]))
    } }
  out <- data.frame(lag = rep(LAGS, each = nd), depth = rep(c(0L, KS), length(LAGS)), acc / length(STARTS)); names(out)[3:9] <- c("Pi", "PiC5", "PiC", "rho", "S", "SC5", "SC")   # PiC5: 5 coarse classes; PiC: 10 fine classes (primary)
  out
}
set.seed(SEED[[MODE]], kind = "L'Ecuyer-CMRG")
t0 <- Sys.time()
rows <- mclapply(seq_len(nrow(GRID)), function(g) { p <- GRID[g, ]
  do.call(rbind, lapply(seq_len(REPS), function(rep) { m <- one_run(p$beta, p$delta, p$lambda, p$gamma)
    data.frame(p, rep = rep, m, row.names = NULL) })) }, mc.cores = CORES, mc.set.seed = TRUE)
D <- do.call(rbind, rows); write.csv(D, file.path(OUT, "rep.csv"), row.names = FALSE)
rows2 <- mclapply(seq_len(nrow(SENS)), function(g) { p <- SENS[g, ]
  do.call(rbind, lapply(seq_len(SENS_REPS), function(rep) { m <- one_run(p$beta, p$delta, p$lambda, p$gamma, p$N, p$EXIT)
    data.frame(p, rep = rep, m, row.names = NULL) })) }, mc.cores = CORES, mc.set.seed = TRUE)
write.csv(do.call(rbind, rows2), file.path(OUT, "sens.csv"), row.names = FALSE)
cat(MODE, ": ", nrow(GRID), "cells x", REPS, "reps in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
