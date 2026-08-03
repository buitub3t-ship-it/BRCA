source("config.R")

info <- capture.output(sessionInfo())
writeLines(info, file.path(RESULTS_DIR, "05_sessionInfo.txt"))
message("Saved results/05_sessionInfo.txt")
