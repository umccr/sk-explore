# Required packages
library(webr)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(plotly)
library(dplyr)

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
  category = c("UnSuccessful", "UnSuccessful", "Successful"),
  value = c(5, 4, 45),
  stringsAsFactors = FALSE
)

# Pie Donut chart using webr
PieDonut(pi_col, aes(category, outcome, count=value), color = "white", title = "Samples overview", ratioByGroup = FALSE,
         showRatioDonut = FALSE, maxx=1.5, labelposition = 1, pieAlpha = 0.6, r1 = 0.9, addDonutLabel = FALSE)


# Pie Donut using plotly
# Data preparation
# Inner layer data (categories)
inner_data <- pi_col %>%
  group_by(category) %>%
  summarize(value = sum(value)) %>%
  mutate(
    labels = category, # Labels for the inner layer
    parents = "Category" # Parent for all categories
  )

# Outermost center layer
center_data <- data.frame(
  labels = "Category",
  parents = "",
  value = sum(pi_col$value)
)

# Outer layer data (outcomes)
outer_data <- pi_col %>%
  mutate(
    labels = outcome,
    parents = category
  )

# Combine all layers
sunburst_data <- bind_rows(center_data, inner_data, outer_data)

# Custom colors
custom_colors <- c(
  "Category" = "white",
  "UnSuccessful" = "#FF6666",
  "Successful" = "#66CC99",
  "Failed sample QC (not sequenced)" = "#FF9999",
  "Failed WGS QC (sequenced)" = "#FF9999",
  "WGTS reporting" = "#66CC99"
)


sunburst_data <- sunburst_data %>%
  mutate(color = custom_colors[labels])

plot_ly() %>%
  config(
    toImageButtonOptions = list(
      format = "svg",
      filename = "myplot",
      width = 600,
      height = 700
    )
  )

fig <- plot_ly(
  sunburst_data,
  labels = ~labels,
  parents = ~parents,
  values = ~value,
  type = 'sunburst',
  branchvalues = 'total',
  marker = list(colors = ~color), # Apply custom colors
  textinfo = 'label+percent entry', # Show labels and percentages
  rotation = 90
)

fig <- fig %>%
  layout(title = "Samples Overview") %>%
  config(
    toImageButtonOptions = list(
      format = "svg",
      filename = "fig",
      width = 600,
      height = 600
    )
  )
# Print with labels
fig


fig_nolabel <- plot_ly(
  sunburst_data,
  labels = ~labels,
  parents = ~parents,
  values = ~value,
  type = 'sunburst',
  branchvalues = 'total',
  marker = list(colors = ~color),         # Apply custom colors
  textinfo = 'none',
  rotation = 153
)

fig_nolabel <- fig_nolabel %>%
  config(
    toImageButtonOptions = list(
      format = "svg",
      filename = "fig_nolabel",
      width = 600,
      height = 600
    )
  )

# Print without labels
fig_nolabel
