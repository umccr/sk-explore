library(magrittr)
library(reshape2)
library(plotly)
library(tidyr)
library(here)

# Path to combined expression files
file_path_organoid <- "/Users/kanwals/UMCCR/research/projects/brinapant/Avner RNAseq/Avner_Birinapant_Organoid_filtered_wo_rep_CPM_TMM_filtered_log_annot.txt"
file_path_tissue <- "/Users/kanwals/UMCCR/research/projects/brinapant/Avner RNAseq/Avner_Birinapant_Tissue_WGS_matched_filtered_wo_rep_CPM_TMM_filtered_log_annot.txt"


# Read the content of the text files
data_organoid <- read.table(file_path_organoid, header = TRUE, sep = "\t")
data_tissue <- read.table(file_path_tissue, header = TRUE, sep = "\t")

# Extract expression for gene of interest (MLKL) in this case from tissue and organoid
gene_expression_organoid <- data_organoid[data_organoid$X == "REV3L", ] %>%
     subset(select=-X)
gene_expression_tissue <- data_tissue[data_tissue$X == "REV3L", ] %>%
     subset(select=-X)

# Find the unique column names from both data frames
all_columns <- unique(c(names(gene_expression_organoid), names(gene_expression_tissue)))

# Add missing columns to gene_expression_organoid with NA values
missing_columns_organoid <- base::setdiff(all_columns, names(gene_expression_organoid))
gene_expression_organoid[missing_columns_organoid] <- NA

# Add missing columns to gene_expression_tissue with NA values
missing_columns_tissue <- base::setdiff(all_columns, names(gene_expression_tissue))
gene_expression_tissue[missing_columns_tissue] <- NA

# Add a unique identifier to each data frame
gene_expression_organoid$source <- "organoid"
gene_expression_tissue$source <- "tissue"

# Combine the data frames and convert to long format using gather
combined_data <- rbind(gene_expression_organoid, gene_expression_tissue)
combined_data_long <- gather(combined_data, key = "variable", value = "value", -source)

# Convert the values column to numeric (if not already)
combined_data_long$value <- as.numeric(combined_data_long$value)

# Plot results
p <- plot_ly(data = combined_data_long, x = ~variable, y = ~value, color = ~source, type = "bar") %>%
  plotly::layout(title = "Expression of REV3L",
         yaxis = list(title = "Expression Value"),
         xaxis = list(title = "Sample Type"))

htmlwidgets::saveWidget(as_widget(p), "/Users/kanwals/UMCCR/research/projects/brinapant/Avner RNAseq/expression_REV3L.html")
