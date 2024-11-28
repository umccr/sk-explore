# Required packages
library(plotly)
library(networkD3)

csv_malignant <- "
source,    target,    value,  group
Malignant_1,     Malignant_2,     12, mono
Malignant_1,     Malignant_2,     14, pleo
Malignant_1,     Malignant_2,     5, rcc
Malignant_1,     Malignant_2,     3, swi
Uncertain,     Malignant_2,     4,  mono
Uncertain,     Benign,     2, rcc
Uncertain,     Intermediate,   1, mono
Uncertain,   Uncertain_2,   1, swi
Uncertain,   Uncertain_2,   1, mono
Uncertain,   Uncertain-neoplastic,    1, pleo
"

links <- read.csv(text = csv_malignant, header = TRUE, strip.white = TRUE)

# Extract unique node IDs from source and target columns
id <- unique(c(as.character(links$source), as.character(links$target)))

# Create labels by removing the suffix
label <- sub("_[0-9]", "", id)

# Create nodes data frame
nodes <- data.frame(id = id, label = label)

# Assign group based on node label
nodes$group <- ifelse(grepl("Malignant", nodes$label), "Malignant",
                      ifelse(grepl("Intermediate", nodes$label), "Intermediate",
                             ifelse(grepl("Benign", nodes$label), "Benign",
                                    ifelse(grepl("Uncertain", nodes$label), "Uncertain", "Other"))))

# Update links data frame to add id source and target columns
links$IDsource <- match(links$source, nodes$id) - 1
links$IDtarget <- match(links$target, nodes$id) - 1

# Define color scale for groups
color_scale <- 'd3.scaleOrdinal()
.domain(["Malignant", "Intermediate", "Benign", "Uncertain", "mono", "pleo", "rcc", "swi", "Other"])
.range(["#333333", "#858585", "#d9d9d9", "#FFFFFF", "#000584", "#B6008C", "#E52918", "#FFB400", "#808080"])'

# Create Sankey network diagram
p <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "IDsource",
  Target = "IDtarget",
  Value = "value",
  NodeID = "label",
  fontSize = 0,
  colourScale = color_scale,
  LinkGroup = "group",
  NodeGroup = "group",
  iterations = 0
)
