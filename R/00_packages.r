# ---- centralized package management ----
pkgs <- c(
    "tidyverse", "readxl", "plotly", "RColorBrewer",
    "truncnorm", "gt", "pbapply", "highcharter",
    "echarts4r", "bezier", "here", "furrr", "progressr", "pbapply",
    "parallel", "profvis"
  )

install.packages(setdiff(pkgs, rownames(installed.packages())))
invisible(lapply(pkgs, library, character.only = TRUE))
