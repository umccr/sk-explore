# Required packages
library(webr)
library(ggplot2)
library(viridis)
library(hrbrthemes)

# Data
pi_col <- data.frame(
  outcome = c("WGTS reporting", "Failed DNA QC", "Consented to other study",
              "Failed WGS QC", "Lack of Patient Consent", "Low purity", "Failed DNA extraction",
              "Lack of clinician buy-in"),
  category = c("Successsful", "Unsuccessful", "Other",
               "Unsuccessful", "Other", "Unsuccessful", "Unsuccessful",
               "Other"),
  value = c(45, 5, 4,
            3, 3, 1, 1,
            1)
)

# Data update
pi_col <- data.frame(
  outcome = c("Failed sample QC (not sequenced)", "Failed WGS QC (sequenced)", "WGTS reporting"),
  category = c("Successsful", "Unsuccessful", "Unsuccessful", "Unsuccessful"),
  value = c(5, 4, 45)
)

# Pie Donut chart
PieDonut(pi_col, aes(category, outcome, count=value), title = "Samples overview", ratioByGroup = FALSE,
         showRatioDonut = FALSE, maxx=1.5, labelposition = 1, pieAlpha = 0.6, r1 = 0.9)
