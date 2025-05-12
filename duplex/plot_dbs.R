# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggpubr)

# Set working directory
setwd("~/UMCCR/research/projects/CX5461/")

# Function to read and process sample data
read_sample_data <- function(file_path, sheet_name) {
  read_excel(file_path, sheet = sheet_name) %>%
    mutate(Sample = as.character(Sample))
}

# Function to process sample groups
process_sample_groups <- function(file_path) {
  sample_groups_raw <- read_excel(file_path, sheet = "samplegroup", col_names = FALSE)

  major_groups <- as.character(sample_groups_raw[1, ])
  sub_groups <- as.character(sample_groups_raw[2, ])
  sub_group_desc <- as.character(sample_groups_raw[3, ])

  map_df(seq_len(ncol(sample_groups_raw)), ~{
    samples_in_group <- na.omit(as.character(sample_groups_raw[4:nrow(sample_groups_raw), .x]))
    if (length(samples_in_group) > 0) {
      tibble(
        Sample = samples_in_group,
        MajorGroup = major_groups[.x],
        SubGroup = sub_groups[.x],
        SubGroupDescription = sub_group_desc[.x]
      )
    }
  }) %>%
    filter(SubGroup != "Control")
}

# Function to read and process mutation data
read_mutation_data <- function(sample_name, subgroup_name) {
  file_path <- file.path("./output_default_png_excel", sample_name, "mutational_counts.xlsx")

  if (file.exists(file_path)) {
    read_excel(file_path, sheet = "DBS_counts") %>%
      set_names(c("Context", "Count")) %>%
      mutate(
        Type = paste0(sub(">.*", "", gsub("_", ">", Context)), ">NN"),
        Sample = sample_name,
        SubGroup = subgroup_name
      ) %>%
      group_by(Type, Sample, SubGroup) %>%
      summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")
  } else {
    warning(paste("File not found:", file_path))
    NULL
  }
}

# Function to create plot
create_plot <- function(plot_data, sg) {
  sample_order <- get_sample_order(sg)

  gap_data <- tibble(
    Sample = sample_order[grepl("^gap", sample_order)],
    Type = NA_character_,
    Count = NA_real_,
    Proportion = NA_real_,
    SubGroup = sg
  )

  plot_data <- bind_rows(plot_data, gap_data) %>%
    mutate(Sample = factor(Sample, levels = sample_order)) %>%
    arrange(Sample)

  ggplot(plot_data, aes(x = Sample, y = Proportion, fill = Type)) +
    geom_col(data = ~ filter(., !grepl("^gap", Sample)), width = 0.7, color = "black", linewidth = 0.2) +
    geom_vline(xintercept = which(grepl("^gap", sample_order)), linetype = "dashed", color = "black", linewidth = 0.5) +
    labs(
      title = paste0("Single Base Substitution Proportion - ", sg),
      x = "Sample",
      y = "Double Base Substitution Proportion",
      fill = "Mutation Type"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 1.5, color = "black", fill = NA),
      axis.line = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(
      values = c(
        "AC>NN" = "#1f77b4", "AT>NN" = "#ff7f0e", "CC>NN" = "#2ca02c",
        "CG>NN" = "#d62728", "CT>NN" = "#9467bd", "GC>NN" = "#8c564b",
        "TA>NN" = "#e377c2", "TC>NN" = "#7f7f7f", "TG>NN" = "#bcbd22",
        "TT>NN" = "#17becf"
      ),
      na.value = "white"
    ) +
    scale_x_discrete(drop = FALSE, labels = ~ ifelse(grepl("^gap", .), "",
                                                     ifelse(. == "Control", "Control",
                                                            ifelse(grepl("Base$", ., ignore.case = TRUE), "Baseline",
                                                                   sub(".*BMN|.*BMTTu|.*BMTu|.*PBMC|.*SkN|.*C4D15|.*SkTu", "", .)))))
}

# Function to get sample order - updated according to May discussion
get_sample_order <- function(sg) {
  switch(sg,
         "PMC-07-CLL-RT" = c("07BMNBase", "07BMNC1D1", "gap1", "07BMTuBase", "07BMTuC1D2", "gap2", "07PBMCBase", "07PBMCEOT"),
         "PMC-09-MM" = c("09BMNBase", "09BMNEOT", "gap1", "09BMTuBase", "09BMTuEOT", "gap2", "09PBMCBase", "09PBMCEOT"),
         "PMC-14-MM" = c("14BMTTuBase", "14BMTTuC1D2", "14BMTTuC4D15EOT", "gap1", "14PBMCBase", "14PBMCEOT"),
         "PMC-01-CTCL" = c("01SkNBase", "01SkNC1D9", "gap1", "01SkTuBase", "01SkTuC1D9", "gap2", "01PBMCBase", "01PBMCEOT"),
         character(0)
  )
}

# Main execution - assumes signature analysis has been performed prior to running this to have sample name, snv_count,
# indel_count, dnv_count values
main <- function() {
  excel_file <- "mutation_counts_summary.xlsx"

  sample_counts <- read_sample_data(excel_file, "samplecount")
  sample_group_info <- process_sample_groups(excel_file)

  all_data <- map2_dfr(sample_group_info$Sample, sample_group_info$SubGroup, read_mutation_data)

  all_data_df <- all_data %>%
    group_by(Sample) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()

  subgroups <- unique(all_data_df$SubGroup)

  for (sg in subgroups) {
    plot_data <- all_data_df %>% filter(SubGroup == sg)
    p <- create_plot(plot_data, sg)
    print(p)
    # ggsave(paste0("./figure/dbs/", sg, "dbs_substitution_gap_plot.pdf"), plot = p, width = 6, height = 6)
  }
}

# Call the main function
main()
