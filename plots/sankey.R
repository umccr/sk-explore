# Required package
library(plotly)

# Define nodes and links data
nodes <- data.frame(name = c(
  "Start of the Study",
  "WGTS Report Successful (44)",
  "Low Purity (1)",
  "Failed WGS QC (3)",
  "Failed DNA QC (5)",
  "Failed DNA Extraction (1)",
  "Consented to Other Study (4)",
  "Lack of Patient Consent (3)",
  "Lack of Clinician Buy-In (1)"
))

links <- data.frame(
  source = c(0, 0, 0, 0, 0, 0, 0, 0, 0),
  target = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  value = c(44, 1, 3, 5, 1, 4, 3, 1, 0)
)

# Sorting links by value in descending order
links <- links[order(-links$value), ]

# Create the Sankey plot
sankey_plot <- plot_ly(
  type = "sankey",
  orientation = "h",
  node = list(
    pad = 15,
    thickness = 20,
    line = list(color = "black", width = 0.5),
    label = nodes$name
  ),
  link = list(
    source = links$source,
    target = links$target,
    value = links$value
  )
)

# Customize the layout - can be further explored
layout <- list(
  title = "Sortable Sankey Diagram with Values",
  font = list(size = 12)
)

# Create and display the final plot
sankey_plot <- sankey_plot %>% layout(layout)




