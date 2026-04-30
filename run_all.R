repo <- "/workspaces/Replication-of-Glynn-CI-TSCS"
setwd(repo)
dir.create(file.path(repo, "output"), showWarnings = FALSE)
message("=== [1/3] burgoon.R ===")
source(file.path(repo, "code", "burgoon.R"))
message("=== [2/3] swank-steinmo.R ===")
source(file.path(repo, "code", "swank-steinmo.R"))
message("=== [3/3] causal-tscs-sim.R ===")
source(file.path(repo, "code", "causal-tscs-sim.R"))
message("=== Done! Check output/ for figures and tables ===")

