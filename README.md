# Rank or Continuity? Continuity-Adjusted Excess Persistence in Ranked Lists

Simulation code for the paper. R, base packages only. Run from this folder:

    SIM_MODE=confirm SIM_OUTPUT_ROOT=out SIM_CORES=2 Rscript persistence_continuity_sim.R
    SIM_MODE=confirm SIM_OUTPUT_ROOT=out Rscript make_tables_continuity.R
    SIM_MODE=confirm SIM_OUTPUT_ROOT=out Rscript make_figures_continuity.R

* `docs/PROTOCOL_continuity_v5.md`: frozen protocol (hypotheses, margins, amendment).
* `docs/FREEZE_continuity.txt`: SHA-256 checksums of the simulation script, tables script, and protocol.
* Seeds: pilot 1101, exploration 7931, confirmation 9523. Output folders are created automatically; figures are PNG files.
