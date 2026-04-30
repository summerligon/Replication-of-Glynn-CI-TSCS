repo <- "/workspaces/Replication-of-Glynn-CI-TSCS"
code_dir <- file.path(repo, "code")
data_dir <- file.path(repo, "data")
out_dir  <- file.path(repo, "output")
source_line <- paste0('source("', code_dir, '/panel-utils.R")')

patch_file <- function(filepath, replacements) {
  lines <- readLines(filepath)
  for (r in replacements) lines <- gsub(r$from, r$to, lines, fixed=TRUE)
  writeLines(lines, filepath)
  message("Patched: ", basename(filepath))
}

patch_file(file.path(code_dir,"causal-tscs-sim.R"), list(
  list(from='source("panel-utils.R")', to=source_line),
  list(from='cairo_pdf(filename =', to='pdf(file ='),
  list(from='cairo_pdf(', to='pdf('),
  list(from=', family = "Minion Pro"', to=''),
  list(from=', family = "Work Sans Regular"', to='')
))