# Protocol v5: does rank add information beyond continuity?

Written 2026-10-17 before the confirmatory run. Scripts: `persistence_continuity_sim.R`, `make_tables_continuity.R` (SHA-256 checksums in `FREEZE_continuity.txt`). This protocol replaces protocol v4 (the earlier study used 5 duration classes, depth-weighted values only, and also analyzed SJR lists; none of that analysis is used here).

## 1. Question
Excess persistence compares the continued membership of the top-k members of a list with a benchmark that treats exits as random. A value above zero is usually read as information in rank. It can also arise when members differ in how likely they are to stay, because the benchmark gives every member the same exit probability. The study asks when that happens, at which list depths (fixed top-k cohorts and the depth-weighted summary), and whether a benchmark that conditions on uninterrupted duration (consecutive periods on the list) isolates the information in rank.

## 2. Measures (lag h, start period t)
* Pi: excess persistence against the random-exit benchmark S(h) = prod_s (1 - X_s / N_(s-1)).
* PiC (primary): continuity-adjusted excess persistence. Exit rates are estimated for each period within 10 classes of uninterrupted duration (1, 2, 3, 4, 5-6, 7-8, 9-11, 12-16, 17-24, 25+ periods at the end of the previous period).
* PiC5: the same with 5 coarse classes (1, 2, 3-4, 5-8, 9+), reported as a granularity check.
* Depth: fixed top-k cohorts k = 5, 10, 20, 50 (N = 200 members, so k = 2.5%, 5%, 10%, 25% of the list), and the depth-weighted summary over all positions (weights proportional to k), coded depth 0. Hypotheses are tested at depth 0 and k = 10, 20, 50; k = 5 is descriptive because a cohort of five members gives noisy estimates.

## 3. Simulation
Ranked list of N = 200 members, 46 periods, 10% of members leave at each transition (exactly 20), replaced by new members. Score = permanent component z + lambda * log(duration) + AR(1) transitory component (coefficient 0.7, SD 0.7). Exit weights exp(beta * r + gamma * u - delta * log(duration)), r the relative rank (0 best), u an unobserved frailty. No re-entry. Start periods 21, 24, 27, 30, 33, 36 (averaged); lags 1, 2, 3, 5, 10; lag 5 is primary.
* beta (rank effect) in {0, 1, 2}; delta (duration dependence of exit) in {0, 0.5, 1}; lambda (duration-score link) in {0, 0.5, 1}; gamma (frailty) in {0, 1}: 54 cells, 500 replications each. A sensitivity grid varies N (50, 800) and the exit rate (5%, 20%) for four cells, 200 replications.
* Seeds: pilot 1101, exploration 7931, confirmation 9523.

## 4. Hypotheses (lag 5)
Tests are intersection-union tests over cells (p is the largest cell p-value, one-sided or equivalence), with Holm adjustment across the eleven rows, alpha = .05. Equivalence uses 90% intervals.
* H1 (Propositions 1 and 2, depth-weighted). With beta = delta = gamma = 0 (any lambda), Pi and PiC lie within +/-0.02 of zero.
* H2 (the same at fixed depths). With beta = delta = gamma = 0, Pi and PiC lie within +/-0.04 of zero at k = 10, 20, 50.
* H3 (false persistence). With beta = 0 and delta > 0 or gamma > 0 (15 cells), Pi exceeds 0.05 at depth 0 and at k = 10, 20, 50.
* H4 (adjustment, depth-weighted). In the same cells, PiC lies within +/-0.03 of zero.
* H5 (adjustment at fixed depths). In the same cells, PiC lies within +/-0.08 of zero at k = 10, 20, 50.
* H6 (detection). With beta > 0 (36 cells), PiC exceeds 0.03 at depth 0 and 0.05 at k = 10, 20, 50.
* H7 (dose response). Within every combination of delta, lambda, and gamma, PiC rises from beta 0 to 1 and from 1 to 2 at depth 0 and at k = 10, 20, 50.
* H8 (granularity). With beta = 0, delta > 0, lambda > 0 (8 cells), PiC5 - PiC exceeds 0.03 at k = 20.
* H9 (depth profile). With beta > 0, the mean of PiC(10) - PiC(50) exceeds 0.02 at each beta level.
* H10 (lag). In the H3 cells, Pi(5) - Pi(1) exceeds 0.05 at depth 0.
* H11 (cost of adjustment). With beta > 0 and delta = lambda = gamma = 0, Pi - PiC lies within +/-0.06 at depth 0.

## 5. Disclosures
* The grid values, the 10 fine classes, the margins, and the eleven hypotheses were set after a pilot and two exploratory runs (seed 7931, 150 replications per cell). In the first exploratory run the 5 coarse classes left a residual false persistence of up to 0.20 at fixed depths, which led to adopting 10 classes as the primary adjustment and keeping the 5-class version as H8. The number of start periods was raised from four to six to reduce noise at fixed depths.
* The margins were chosen with the exploratory estimates in view (smallest H3 estimate 0.122; largest H4 estimate 0.005; largest H5 estimate 0.031; smallest H6 estimate 0.044 at depth 0 and 0.089 at fixed depths; smallest H8 estimate 0.044; smallest H9 estimate 0.054; smallest H10 estimate 0.061; largest H11 estimate 0.044).
* An earlier study in this project (protocol v4, 4 start periods, 5 classes, depth-weighted values, SJR application) and protocols v1 to v3 used different designs; their data are not used here. The confirmatory run tests the margins on new data, but the design is not fully independent of the exploratory runs.

## 6. Amendment after the confirmatory run (2026-10-17)
The confirmatory run was completed and the first run of the tables script returned missing values for H6, H7, and H9. The cause is replications in which PiC (or PiC5) is undefined, because a benchmark equals 1 in a start period (all members of the cohort are in classes with no recorded exits), which gives 0/0. Such values are dropped, as stated for the earlier protocol. They occur mainly at lag 1 and at small k, and almost never at lag 5 (counts reported in the manuscript). The tables script was amended to drop them (`na.rm`); no hypothesis, margin, or test was changed. The checksum of the first version of `make_tables_continuity.R` was 75bc1c9256add0e60f4fe680dca8881a72c7d1d03c6116bdaf880f5ce4e63552; `FREEZE_continuity.txt` records the amended version.
