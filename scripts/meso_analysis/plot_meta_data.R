# Load required libraries
library(readxl)
library(VariantAnnotation)
library(maftools)
library(data.table)
library(ComplexHeatmap)
library(circlize)

# Set working directory
setwd("~/UMCCR/research/projects/mesogenome/data")

# Load and preprocess clinical data
load_and_preprocess_clinical_data <- function() {
  df <- read_excel("clinical.xlsx")
  
  df$Age <- df$`Age at Diagnosis`
  df$"Sample purity" <- cut(
    as.numeric(df$`Sample purity`),
    breaks = c(0, 0.3, 0.5, 0.7, Inf),
    labels = c("<30%", "30-50%", "50-70%", ">70%"),
    right = FALSE 
  )
  
  df$Survival_yr <- as.numeric(df$`Survival from Diagnosis Including Alive (months)`) / 12
  df$"Survival (year)" <- cut(
    df$Survival_yr,
    breaks = c(0, 1, 2, 3, Inf),
    labels = c("<1", "1-2", "2-3", "3+"),
    right = FALSE
  )
  
  clinical_vars <- c(
    "ID", "Age", "Gender", "Asbestos Exposure", "Smoking history",
    "Histology", "Meso type", "Sample purity", "Prior treatment", "Survival (year)"
  )
  cat_df <- df[, clinical_vars]
  
  cat_df$`Prior treatment`[cat_df$`Prior treatment` == "Post-Tx"] <- "No"
  cat_df$`Prior treatment`[cat_df$`Prior treatment` == "Pre-Tx"] <- "Yes"
  
  return(cat_df)
}

# Process VCF files
process_vcf_files <- function() {
  vcf_files <- list.files(path = "somatic-variants", pattern = "*.vcf.gz", full.names = TRUE)
  samples <- sapply(vcf_files, function(f) strsplit(basename(f), "__")[[1]][1])
  
  count_mutations <- function(vcf_file) {
    vcf <- readVcf(vcf_file, genome = "hg38")  
    ref_len <- nchar(ref(vcf))
    alt_list <- alt(vcf)
    alt_len <- sapply(alt_list, function(alt) if (length(alt) == 0) 0 else max(nchar(as.character(alt))))
    
    counts <- c(
      "Single nucleotide" = sum(ref_len == 1 & alt_len == 1),
      "Dinucleotide" = sum(ref_len == 2 & alt_len == 2),
      "Trinucleotide" = sum(ref_len == 3 & alt_len == 3),
      "Insertion" = sum(alt_len > ref_len),
      "Deletion" = sum(ref_len > alt_len)
    )
    
    return(counts)
  }
  
  mutation_data <- lapply(vcf_files, count_mutations)
  mutation_df <- do.call(rbind, mutation_data)
  rownames(mutation_df) <- samples
  return(as.data.frame(mutation_df))
}

# Process structural variant files
process_sv_files <- function() {
  folder_path <- "manta-sv"
  tsv_files <- list.files(path = folder_path, pattern = "*.tsv", full.names = TRUE)
  samples <- sapply(tsv_files, function(f) strsplit(basename(f), "__")[[1]][1])
  
  count_sv_types <- function(tsv_file) {
    data <- fread(tsv_file, sep = "\t", header = TRUE)
    sv_types <- data$svtype
    return(table(sv_types))
  }
  
  sv_data <- lapply(tsv_files, count_sv_types)
  all_sv_types <- unique(unlist(lapply(sv_data, names)))
  
  sv_df <- data.frame(matrix(0, nrow = length(sv_data), ncol = length(all_sv_types)))
  colnames(sv_df) <- all_sv_types
  rownames(sv_df) <- samples
  
  for (i in 1:length(sv_data)) {
    for (sv_type in names(sv_data[[i]])) {
      sv_df[i, sv_type] <- sv_data[[i]][sv_type]
    }
  }
  
  return(sv_df)
}

# Prepare data for heatmap
prepare_heatmap_data <- function(cat_df, mutation_df, sv_df) {
  cat_df <- as.data.frame(cat_df)
  rownames(cat_df) <- cat_df$ID
  cat_df <- cat_df[, -1]  # Remove the ID
  
  common_samples <- intersect(rownames(cat_df), rownames(mutation_df))
  common_samples <- common_samples[common_samples != "SBJ05849"]
  
  cat_df <- cat_df[common_samples, , drop = FALSE]
  mutation_df <- mutation_df[common_samples, , drop = FALSE]
  sv_df <- sv_df[common_samples, , drop = FALSE]
  
  age_vac <- cat_df$Age
  cat_df <- cat_df[, !colnames(cat_df) %in% "Age"]
  
  dummy_mat <- matrix(age_vac, nrow = 1, ncol = nrow(cat_df))
  colnames(dummy_mat) <- rownames(cat_df)
  
  return(list(cat_df = cat_df, mutation_df = mutation_df, sv_df = sv_df, age_vac = age_vac, dummy_mat = dummy_mat))
}

# Define color schemes
define_color_schemes <- function() {
  custom_colors <- list(
    "Gender" = c("M" = "#3B4992", "F" = "#EE4B2B"),
    "Asbestos Exposure" = c("Yes" = "#F28E2B", "No" = "#E15759", "Unknown" = "#BAB0AC"),
    "Smoking history" = c("Yes" = "#76B7B2", "No" = "#59A14F", "Unknown" = "#BAB0AC"),
    "Histology" = c("Epithelioid" = "#AF7AA1", "Sarcomatoid" = "#4E79A7", "Biphasic" = "#D37295", "WDPMT" = "#76448A"),
    "Meso type" = c("Pleural" = "#F2C14E", "Peritoneal" = "#8C564B", "Other" = "#BAB0AC"),
    "Sample purity" = c("<30%" = "#A0CBE8", "30-50%" = "#FFBF86", "50-70%" = "#FF8C94", ">70%" = "#E15759"),
    "Prior treatment" = c("Yes" = "#FF9F1C", "No" = "#2F4B26"),
    "Survival (year)" = c("<1" = "#E15759", "1-2" = "#F28E2B", "2-3" = "#59A14F", "3+" = "#76B7B2")
  )
  
  mutation_colors <- c(
    "Single nucleotide" = "#0077BB", "Dinucleotide" = "#009988", "Trinucleotide" = "#EE7733",
    "Insertion" = "#CC3311", "Deletion" = "#882255"
  )
  
  sv_colors <- c("BND" = "#2F4B26", "DEL" = "#D62728", "DUP" = "#9467BD", "INS" = "#E377C2")
  
  return(list(custom_colors = custom_colors, mutation_colors = mutation_colors, sv_colors = sv_colors))
}

# Create heatmap
create_heatmap <- function(heatmap_data, color_schemes) {
  mutation_matrix <- as.matrix(heatmap_data$mutation_df)
  sv_matrix <- as.matrix(heatmap_data$sv_df)
  
  bottom_ha <- HeatmapAnnotation(
    mutation = anno_barplot(mutation_matrix, baseline = 0, which = "column", border = FALSE,
                            bar_width = 0.92, beside = FALSE, gp = gpar(fill = color_schemes$mutation_colors, col = NA),
                            axis = TRUE, add_numbers = FALSE, height = unit(6, "cm"),
                            axis_param = list(side = "left", gp = gpar(col = "black", lwd = 1)), add_x_axis = TRUE),
    sv = anno_barplot(sv_matrix, baseline = 0, which = "column", border = FALSE,
                      bar_width = 0.92, beside = FALSE, gp = gpar(fill = color_schemes$sv_colors, col = NA),
                      axis = TRUE, add_numbers = FALSE, height = unit(6, "cm"),
                      axis_param = list(side = "left", gp = gpar(col = "black", lwd = 1)), add_x_axis = TRUE),
    gap = unit(5, "mm")
  )
  
  mutation_legend <- Legend(labels = names(color_schemes$mutation_colors), title = "Somatic Variants", 
                            legend_gp = gpar(fill = color_schemes$mutation_colors), direction = "horizontal")
  
  sv_legend <- Legend(labels = names(color_schemes$sv_colors), title = "Structural Variants", 
                      legend_gp = gpar(fill = color_schemes$sv_colors), direction = "horizontal")
  
  ha <- HeatmapAnnotation(df = heatmap_data$cat_df, col = color_schemes$custom_colors,
                          gp = gpar(col = "white", lwd = 0.3), border = FALSE, gap = unit(0.3, "points"))
  
  color_fun <- colorRamp2(c(min(heatmap_data$age_vac, na.rm = TRUE), max(heatmap_data$age_vac, na.rm = TRUE)), 
                          c("#a6cee3", "#1f78b4"))
  
  ht <- Heatmap(heatmap_data$dummy_mat,
                top_annotation = ha,
                bottom_annotation = bottom_ha,
                col = color_fun,
                name = "Age",
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 10),
                heatmap_legend_param = list(title = "Age",
                                            at = round(seq(min(heatmap_data$age_vac), max(heatmap_data$age_vac), length.out = 5), 1),
                                            direction = "horizontal"),
                height = unit(0.5, "cm")
  )
  
  return(list(ht = ht, mutation_legend = mutation_legend, sv_legend = sv_legend))
}

# Main execution
main <- function() {
  cat_df <- load_and_preprocess_clinical_data()
  mutation_df <- process_vcf_files()
  sv_df <- process_sv_files()
  
  heatmap_data <- prepare_heatmap_data(cat_df, mutation_df, sv_df)
  color_schemes <- define_color_schemes()
  
  heatmap_objects <- create_heatmap(heatmap_data, color_schemes)
  
  # Draw and save the heatmap
  pdf("heatmap_output_age_no_SBJ05849.pdf", width = 12, height = 14)
  draw(heatmap_objects$ht, 
       annotation_legend_list = list(mutation = heatmap_objects$mutation_legend, 
                                     sv = heatmap_objects$sv_legend), 
       annotation_legend_side = "top")
  dev.off()
}

# Run the main function
main()