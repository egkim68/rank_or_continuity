# FIGURES v5 (PNG only; separate from the simulation). Reads tables/canonical_by_lag.csv written by make_tables_continuity.R.
MODE <- Sys.getenv("SIM_MODE", unset = "confirm"); ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-persistence")
OUT_DIR <- file.path(ROOT, paste0("continuity_", MODE)); FIG_DIR <- file.path(OUT_DIR, "figures"); dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
can <- read.csv(file.path(OUT_DIR, "tables", "canonical_by_lag.csv")); colP <- "#b5532a"; colT <- "#1f6f78"; colG <- "#b8932f"
lab <- c("Random exit" = "Random exit", "Duration-dependent exit" = "Duration-dependent\nexit", "Frailty only" = "Frailty only",
         "Continuity link and duration-dependent exit" = "Continuity link and\nduration-dependent exit", "Rank effect" = "Rank effect",
         "Rank effect, continuity channels" = "Rank effect,\ncontinuity channels", "Strong rank effect" = "Strong rank effect",
         "Strong rank effect, continuity channels" = "Strong rank effect,\ncontinuity channels")
# Figure 2: Pi and PiC at lag 5 in eight scenarios, depth-weighted and top 10, 20, 50
c5 <- can[can$lag == 5, ]; labs <- names(lab); lab2 <- c("Random exit","Duration-\ndependent exit","Frailty only","Continuity link,\nduration-\ndependent exit","Rank effect","Rank effect,\ncontinuity\nchannels","Strong rank\neffect","Strong rank\neffect, continuity\nchannels"); names(lab2) <- labs; n <- length(labs); y <- rev(seq_len(n)); dps <- c(0, 10, 20, 50); ttl <- c("Depth-weighted", "Top 10", "Top 20", "Top 50")
png(file.path(FIG_DIR, "Fig_scenarios.png"), width = 5200, height = 4300, res = 1000); layout(rbind(1:4, c(5, 5, 5, 5)), widths = c(2.3, 1, 1, 1), heights = c(8, 1.5))
for (j in seq_along(dps)) { d <- c5[c5$depth == dps[j], ]; d <- d[match(labs, d$label), ]
  par(mar = c(5.6, if (j == 1) 7.2 else 0.5, 2.4, 0.5), mgp = c(3.4, 0.7, 0))
  plot(NA, xlim = c(-0.05, 1), ylim = c(0.5, n + 0.5), yaxt = "n", xlab = "Excess\npersistence", ylab = "", main = ttl[j], cex.main = 1.05, cex.lab = 1, cex.axis = 0.95)
  abline(v = 0, col = "grey50"); abline(h = y, col = "grey90"); if (j == 1) axis(2, y, lab2[labs], las = 1, cex.axis = 0.95, tick = FALSE)
  points(d$Pi, y + 0.17, pch = 17, col = colP, cex = 0.9); points(d$PiC, y - 0.17, pch = 16, col = colT, cex = 0.9) }
par(mar = c(0, 0, 0, 0)); plot.new(); legend("center", c("Excess persistence (random-exit benchmark)", "Continuity-adjusted excess persistence"), pch = c(17, 16), col = c(colP, colT), ncol = 1, bty = "n", cex = 1.05, pt.cex = 0.95); dev.off()
# Figure 3: lag profiles in four scenarios, top 20 cohort
sc <- c("Duration-dependent exit", "Frailty only", "Rank effect", "Rank effect, continuity channels")
png(file.path(FIG_DIR, "Fig_lag_profiles.png"), width = 6400, height = 4900, res = 1320); layout(matrix(c(1, 2, 3, 4, 5, 5), 3, 2, byrow = TRUE), heights = c(6, 6, 1.2)); par(mar = c(3.3, 3.9, 3.7, 1.0), mgp = c(2.1, 0.7, 0))
for (s in sc) { d <- can[can$label == s & can$depth == 20, ]; d <- d[order(d$lag), ]; plot(NA, xlim = c(1, 10), ylim = c(-0.05, 1), xlab = "Lag (periods)", ylab = "Excess persistence", main = sub(", ", ",\n", s), cex.main = 1.0, cex.lab = 1.0, cex.axis = 1.0)
  abline(h = 0, col = "grey50", lwd = 1); lines(d$lag, d$Pi, col = colP, lwd = 0.55, lty = 2); points(d$lag, d$Pi, pch = 17, col = colP, cex = 0.8); lines(d$lag, d$PiC, col = colT, lwd = 0.55); points(d$lag, d$PiC, pch = 16, col = colT, cex = 0.8) }
par(mar = c(0, 0, 0, 0)); plot.new(); legend("center", c("Excess persistence", "Continuity-adjusted excess persistence"), pch = c(17, 16), lty = c(2, 1), lwd = 0.55, col = c(colP, colT), ncol = 2, bty = "n", cex = 0.95, pt.cex = 1.1, seg.len = 2.5); dev.off()
# Figure 4: profile over fixed depth k at lag 5, with 5 coarse and 10 fine classes
sc4 <- c("Duration-dependent exit", "Continuity link and duration-dependent exit", "Rank effect", "Rank effect, continuity channels"); kk <- c(5, 10, 20, 50)
png(file.path(FIG_DIR, "Fig_depth_profiles.png"), width = 4800, height = 3500, res = 880); layout(matrix(c(1, 2, 3, 4, 5, 5), 3, 2, byrow = TRUE), heights = c(5, 5, 1.2)); par(mar = c(4.2, 4.4, 2.6, 0.5), mgp = c(2.6, 0.8, 0))
for (s in sc4) { d <- can[can$label == s & can$lag == 5 & can$depth %in% kk, ]; d <- d[order(d$depth), ]; plot(NA, xlim = c(5, 50), ylim = c(-0.05, 1), log = "x", xaxt = "n", xlab = "Top-k cohort (k)", ylab = "Excess persistence", main = lab[[s]], cex.main = 1.0, cex.lab = 1.0, cex.axis = 1.0)
  axis(1, kk, cex.axis = 1.0); abline(h = 0, col = "grey50", lwd = 1)
  lines(d$depth, d$Pi, col = colP, lwd = 0.8, lty = 2); points(d$depth, d$Pi, pch = 17, col = colP, cex = 0.8)
  lines(d$depth, d$PiC5, col = colG, lwd = 0.8, lty = 3); points(d$depth, d$PiC5, pch = 15, col = colG, cex = 0.8)
  lines(d$depth, d$PiC, col = colT, lwd = 0.8); points(d$depth, d$PiC, pch = 16, col = colT, cex = 0.8) }
par(mar = c(0, 0, 0, 0)); plot.new(); legend("center", c("Excess persistence", "Adjusted, 5 classes", "Adjusted, 10 classes"), pch = c(17, 15, 16), lty = c(2, 3, 1), lwd = 0.8, col = c(colP, colG, colT), ncol = 3, bty = "n", cex = 0.95, pt.cex = 1.4, seg.len = 2.5); dev.off()
cat("Figures written to", FIG_DIR, "\n")
