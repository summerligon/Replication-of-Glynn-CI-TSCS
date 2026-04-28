# run_all.R
# ---------------------------------------------------------------
# Master replication script for:
#   Blackwell, M. & Glynn, A. (2018).
#   "How to Make Causal Inferences with Time-Series Cross-Sectional
#    Data under Selection on Observables."
#   American Political Science Review, 112(4), 1067-1082.
#
# Run this from the repo root. All outputs (figures, tables) are
# written to the output/ directory.
#
# Estimated runtimes (single core):
#   burgoon.R        ~5  min
#   swank-steinmo.R  ~3  sec
#   causal-tscs-sim.R ~10 min (faster with parallel cores)
# ---------------------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# If running via Rscript from command line, use instead:
# setwd(getSrcDirectory(function(x) x))

dir.create("output", showWarnings = FALSE)

message("=== [1/3] Running burgoon.R ===")
source("code/burgoon.R")

message("=== [2/3] Running swank-steinmo.R ===")
source("code/swank-steinmo.R")

message("=== [3/3] Running causal-tscs-sim.R ===")
source("code/causal-tscs-sim.R")

message("=== All done! Outputs written to output/ ===")
