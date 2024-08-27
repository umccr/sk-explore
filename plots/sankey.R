# Required packages
library(plotly)
library(networkD3)

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

#links <- data.frame(
#  source = c(0, 0, 0, 0, 0, 0, 0, 0, 0),
#  target = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
#  value = c(44, 1, 3, 5, 1, 4, 3, 1, 0)
#)

links <- data.frame(
  source = c(1,2,3,4,5,2),
  target = c(1,1,1,1,1,2),
  value =  c(8,4,1,1,2)
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
sankey_plot <- sankey_plot %>% plotly::layout(layout)
htmlwidgets::saveWidget(sankey_plot, file = "/Users/kanwals/UMCCR/research/projects/column-pi/sankey_plot.html")

# Plot diagnostic dilemma categories
# Define the data as text and read csv
csv_data <- "
source,    target,    value,  group
One_1,     One_2,     8,    mono
One_1,     One_2,     3,    pleo
One_1,     One_2,     7,    rcc
One_1,     One_2,     2,    swi
Two_1,     One_2,     4,    mono
Two_1,     Two_2,     1,    mono
Two_1,     One_2,     4,    pleo
Two_1,     One_2,     4,    swi
Two_1,     Two_2,     1,    pleo
Three_1,     One_2,   1,    mono
Three_1,     One_2,   4,    pleo
Three_1,     Two_2,   1,  pleo
Four_1,   One_2,   1,   mono
Four_1,   One_2,   2,   pleo
Four_1,   Three_2,   1,   pleo
Five_1,   One_2,    2,    mono
"
csv_data_update <- "
source,    target,    value,  group
One_1,     One_2,     8,    mono
Two_1,     One_2,     4,    mono
Two_1,     Two_2,     1,    mono
Three_1,     One_2,   1,    mono
Four_1,   One_2,   2,   mono
Five_1,   One_2,    2,    mono
One_1,     One_2,     3,    pleo
Two_1,     One_2,     4,    pleo
Two_1,     Two_2,     1,    pleo
Three_1,     One_2,   4,    pleo
Three_1,     Two_2,   1,  pleo
Four_1,   One_2,   1,   pleo
Four_1,   Three_2,   1,   pleo
Four_1,   Four_2, 1,  pleo
One_1,     One_2,     5,    rcc
One_1,     One_2,     2,    swi
Two_1,     One_2,     2,    swi


"
links <- read.csv(text = csv_data_update, header = TRUE, strip.white = TRUE)

# Extract unique node IDs from source and target columns
id <- unique(c(as.character(links$source), as.character(links$target)))

# Create labels by removing the suffix
label <- sub("_[0-9]", "", id)

# Create nodes data frame
nodes <- data.frame(id = id, label = label)

# Add a 'group' column to each node. Here I decided to put all of them in the same group to make them grey
# reference https://r-graph-gallery.com/322-custom-colours-in-sankey-diagram.html
nodes$group <- as.factor(c("my_unique_group"))

# Update links data frame to add id sourcw and target columns
links$IDsource <- match(links$source, nodes$id) - 1
links$IDtarget <- match(links$target, nodes$id) - 1

# Define color scale for groups
color_scale <- 'd3.scaleOrdinal()
                .domain(["mono", "pleo", "rcc", "swi", "my_unique_group"])
                .range(["#69b3a2", "#b3697a", "#b37d69", "#a269b3", "grey"])'

# Create Sankey network diagram
sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "IDsource",
  Target = "IDtarget",
  Value = "value",
  NodeID = "label",
  fontSize = 16,
  colourScale = color_scale,
  LinkGroup = "group",
  NodeGroup = "group"
)

