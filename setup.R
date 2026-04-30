# setup.R
# Run this once to organize all replication files into the correct structure
# ---------------------------------------------------------------

repo <- "/workspaces/Replication-of-Glynn-CI-TSCS"
setwd(repo)

# 1. Create folder structure
message("Creating folders...")
dir.create(file.path(repo, "code"),   showWarnings = FALSE)
dir.create(file.path(repo, "data"),   showWarnings = FALSE)
dir.create(file.path(repo, "output"), showWarnings = FALSE)

# 2. Copy R scripts to code/
r_scripts <- c("burgoon.R", "causal-tscs-sim.R", "panel-utils.R", "swank-steinmo.R")
for (f in r_scripts) {
  from <- file.path(repo, "dataverse_files", "code", "code", f)
  to   <- file.path(repo, "code", f)
  if (file.exists(from)) {
    file.copy(from, to)
    message("Copied ", f, " -> code/")
  } else {
    message("WARNING: ", f, " not found at ", from)
  }
}

# 3. Copy data files to data/
data_files <- c("burgoon.csv", "swank.dta")
for (f in data_files) {
  from <- file.path(repo, "dataverse_files", "data", "data", f)
  to   <- file.path(repo, "data", f)
  if (file.exists(from)) {
    file.copy(from, to)
    message("Copied ", f, " -> data/")
  } else {
    message("WARNING: ", f, " not found at ", from)
  }
}

# 4. Patch file paths in scripts
message("\nPatching file paths in scripts...")

patch_file <- function(filepath, replacements) {
  if (!file.exists(filepath)) {
    message("SKIP (not found): ", filepath); return()
  }
  lines <- readLines(filepath)
  for (r in replacements) lines <- gsub(r$from, r$to, lines, fixed = TRUE)
  writeLines(lines, filepath)
  message("Patched: ", basename(filepath))
}

code_dir <- file.path(repo, "code")
data_dir <- file.path(repo, "data")
out_dir  <- file.path(repo, "output")
source_line <- paste0('source("', code_dir, '/panel-utils.R")')

patch_file(file.path(code_dir, "burgoon.R"), list(
  list(from='source("panel-utils.R")',        to=source_line),
  list(from='read.csv("burgoon.csv"',         to=paste0('read.csv("', data_dir, '/burgoon.csv"')),
  list(from='cairo_pdf(filename =',           to='pdf(file ='),
  list(from='cairo_pdf(',                     to='pdf('),
  list(from=', family = "Minion Pro"',        to=''),
  list(from=', family = "Work Sans Regular"', to=''),
  list(from='"fig6-burgoon.pdf"',             to=paste0('"', out_dir, '/fig6-burgoon.pdf"')),
  list(from='"fig7-burgoon-iptw.pdf"',        to=paste0('"', out_dir, '/fig7-burgoon-iptw.pdf"'))
))

patch_file(file.path(code_dir, "swank-steinmo.R"), list(
  list(from='source("panel-utils.R")',        to=source_line),
  list(from='read.dta("swank.dta")',          to=paste0('foreign::read.dta("', data_dir, '/swank.dta")')),
  list(from='cairo_pdf(filename =',           to='pdf(file ='),
  list(from='cairo_pdf(',                     to='pdf('),
  list(from=', family = "Minion Pro"',        to=''),
  list(from=', family = "Work Sans Regular"', to='')
))

patch_file(file.path(code_dir, "causal-tscs-sim.R"), list(
  list(from='source("panel-utils.R")',        to=source_line),
  list(from='cairo_pdf(filename =',           to='pdf(file ='),
  list(from='cairo_pdf(',                     to='pdf('),
  list(from=', family = "Minion Pro"',        to=''),
  list(from=', family = "Work Sans Regular"', to='')
))

# 5. Create run_all.R
message("\nCreating run_all.R...")
writeLines(paste0(
  'repo <- "', repo, '"
setwd(repo)
dir.create(file.path(repo, "output"), showWarnings = FALSE)
message("=== [1/3] burgoon.R ===")
source(file.path(repo, "code", "burgoon.R"))
message("=== [2/3] swank-steinmo.R ===")
source(file.path(repo, "code", "swank-steinmo.R"))
message("=== [3/3] causal-tscs-sim.R ===")
source(file.path(repo, "code", "causal-tscs-sim.R"))
message("=== Done! Check output/ for figures and tables ===")
'), file.path(repo, "run_all.R"))
message("Created run_all.R")

# 6. Install missing packages
message("\nChecking packages...")
pkgs <- c("foreign","reshape","plyr","plm","geepack","mgcv","MASS","sandwich",
          "lmtest","survey","rms","aod","msm","mvtnorm","boot","stargazer",
          "xtable","ggplot2","coefplot","RColorBrewer","splines","here")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  message("Installing: ", paste(to_install, collapse=", "))
  install.packages(to_install, dependencies=TRUE,
                   repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest")
} else {
  message("All packages already installed!")
}

message("\n=== Setup complete! File structure: ===")
print(list.files(repo, recursive=TRUE, pattern="\\.R$|\\.csv$|\\.dta$"))