# Work around a SIGBUS crash when renv's compiled extension is dyn.load()'d
# on this machine (Apple Silicon); falls back to renv's pure-R implementation.
Sys.setenv(RENV_EXT_ENABLED = "FALSE")

source("renv/activate.R")

if (interactive() && Sys.getenv("TERM_PROGRAM") == "vscode") {
  options(vsc.rstudioapi = TRUE)
  if (requireNamespace("httpgd", quietly = TRUE)) {
    options(vsc.use_httpgd = TRUE)
  }
}
