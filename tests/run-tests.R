# The whole suite. `R CMD check` runs every .R file directly under tests/, so this is the only
# one there — everything else lives in cases/ and is sourced from here, in order.

library(mzLibR)

here <- if (file.exists("cases/harness.R")) "cases" else file.path("tests", "cases")
source(file.path(here, "harness.R"))

for (case in sort(list.files(here, pattern = "^test-.*\\.R$", full.names = TRUE))) {
  cat("== ", basename(case), "\n", sep = "")
  source(case)
}

mz_report()
