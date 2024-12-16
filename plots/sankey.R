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

csv_data_update2 <- "
source,    target,    value,  group
Five_1,   One_2,    2,    mono
Four_1,   Four_2, 1,  pleo
Four_1,   Three_2,   1,   pleo
Four_1,   One_2,   2,   mono
Four_1,   One_2,   1,   pleo
Three_1,     One_2,   1,    mono
Three_1,     One_2,   5,    pleo
Two_1,     One_2,     4,    mono
Two_1,     Two_2,     1,    mono
Two_1,     One_2,     4,    pleo
Two_1,     Two_2,     1,    pleo
Two_1,     One_2,     2,    swi
One_1,     One_2,     8,    mono
One_1,     One_2,     3,    pleo
One_1,     One_2,     7,    rcc
One_1,     One_2,     2,    swi
"

# Malignant potential data
csv_malignant <- "
source,    target,    value,  group
Malignant_1,     Malignant_2,     12, mono
Malignant_1,     Malignant_2,     14, pleo
Malignant_1,     Malignant_2,     5, rcc
Malignant_1,     Malignant_2,     3, swi
Uncertain,     Malignant,     4,  mono
Uncertain,     Benign,     2, rcc
Uncertain,     Intermediate,   1, mono
Uncertain_1,   Uncertain_2,   1, swi
Uncertain_1,   Uncertain_2,   1, mono
Uncertain,   Uncertain-neoplastic,    1, pleo
"

MakeSankey <- function(input_csv) {
  links <- read.csv(text = input_csv, header = TRUE, strip.white = TRUE)

  # Extract unique node IDs from source and target columns
  id <- unique(c(as.character(links$source), as.character(links$target)))

  # Create labels by removing the suffix
  label <- sub("_[0-9]", "", id)
  #label <- c("Five", "Four", "Three", "Two", "One", "Four", "Three", "Two", "One")

  # Create nodes data frame
  nodes <- data.frame(id = id, label = label)

  # Add a 'group' column to each node. And define grey color scale
  # reference https://r-graph-gallery.com/322-custom-colours-in-sankey-diagram.html
  #nodes$group <- as.factor(c("my_unique_group"))
  nodes$group <- factor(nodes$label, levels = c("Five", "Four", "Three", "Two", "One"))

  # Update links data frame to add id source and target columns
  links$IDsource <- match(links$source, nodes$id) - 1
  links$IDtarget <- match(links$target, nodes$id) - 1

  # Define color scale for groups
  # color_scale <- 'd3.scaleOrdinal()
  #               .domain(["mono", "pleo", "rcc", "swi", "my_unique_group"])
  #               .range(["#0072b2", "#cc79a7", "#d55e00", "#e69f00", "grey"])'

  color_scale <- 'd3.scaleOrdinal()
                .domain(["mono", "pleo", "rcc", "swi", "Five", "Four", "Three", "Two", "One"])
                .range(["#0072b2", "#000584", "#E52918", "#FFB400", "#333333", "#5c5c5c", "#858585", "#adadad", "#d9d9d9"])'

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

  # JS <-
  #   '
  #   function(el, x, data){
  #   var svg = d3.select("svg")
  #   // Get the width and height of the SVG container
  #   var svgWidth = svg.node().getBoundingClientRect().width;
  #   var svgHeight = svg.node().getBoundingClientRect().height;
  #
  #   // Calculate the x and y coordinates for the legend position
  #   var legendX = svgWidth - 10; // Adjust this value to position the legend horizontally
  #   var legendY = svgHeight / 2; // Adjust this value to position the legend vertically
  #
  #   // Handmade legend
  #   svg.append("circle").attr("cx", legendX).attr("cy", legendY - 20).attr("r", 6).style("fill", "#5485AB")
  #   svg.append("circle").attr("cx", legendX).attr("cy", legendY).attr("r", 6).style("fill", "#BA4682")
  #   svg.append("text").attr("x", legendX + 10).attr("y", legendY - 20).text("variable M").style("font-size", "15px").attr("alignment-baseline","middle")
  #   svg.append("text").attr("x", legendX + 10).attr("y", legendY).text("variable W").style("font-size", "15px").attr("alignment-baseline","middle")
  #   }
  #   '
  # p <- htmlwidgets::onRender(p,JS)

  return(p)
}


# Plot malignant potential
malignant_lineages <- MakeSankey(csv_data_update2)
htmlwidgets::saveWidget(malignant_lineages, "/Users/kanwals/UMCCR/research/projects/column-pi/sankey_plot-test.html")
webshot2::webshot(url = "/Users/kanwals/UMCCR/research/projects/column-pi/sankey_plot-test.html", file = "/Users/kanwals/UMCCR/research/projects/column-pi/sankey_plot-test.jpg")
