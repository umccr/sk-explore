# Required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggpubr)

# --- Data Loading and Pre-processing ---

load_and_preprocess_data <- function(excel_file) {
  # Load sample counts
  sample_counts <- read_excel(excel_file, sheet = "samplecount") %>%
    mutate(Sample = as.character(Sample))

  # Load and process sample groups
  sample_groups_raw <- read_excel(excel_file, sheet = "samplegroup", col_names = FALSE)

  major_groups <- as.character(sample_groups_raw[1, ])
  sub_groups <- as.character(sample_groups_raw[2, ])
  sub_group_desc <- as.character(sample_groups_raw[3, ])

  group_df_list <- map(seq_len(ncol(sample_groups_raw)), function(i) {
    samples_in_group <- na.omit(as.character(sample_groups_raw[4:nrow(sample_groups_raw), i]))
    if (length(samples_in_group) > 0) {
      data.frame(
        Sample = samples_in_group,
        MajorGroup = major_groups[i],
        SubGroup = sub_groups[i],
        SubGroupDescription = sub_group_desc[i],
        stringsAsFactors = FALSE
      )
    }
  })

  sample_group_info <- bind_rows(group_df_list) %>%
    mutate(Sample = as.character(Sample))

  # Load sample types
  sample_types <- read_excel(excel_file, sheet = "sampletypes")

  # Merge data
  merged_data <- left_join(sample_counts, sample_group_info, by = "Sample") %>%
    mutate(
      MajorGroup = factor(MajorGroup),
      SubGroup = factor(SubGroup)
    )

  return(list(
    merged_data = merged_data,
    sample_types = sample_types
  ))
}

# --- Plotting Function ---

create_mutation_plot_single_axis <- function(data, y_var, title) {
  # Define custom sample order - currently based on serial sampling of patients
  custom_sample_order <- c(
    "07BMNBase", "07BMNC1D1", "07BMTuBase", "07BMTuC1D2", "07PBMCBase", "07PBMCEOT",
    "gap1",
    "09BMNBase", "09BMNEOT", "09BMTuBase", "09BMTuEOT", "09PBMCBase", "09PBMCEOT",
    "gap2",
    "14BMTTuBase", "14BMTTuC1D2", "14BMTTuC4D15EOT", "14PBMCBase", "14PBMCEOT",
    "gap3",
    "01SkNBase", "01SkNC1D9", "01SkTuBase", "01SkTuC1D9", "01PBMCBase", "01PBMCEOT",
    "gap4",
    "Control"
  )
  custom_sample_order <- unique(custom_sample_order[!is.na(custom_sample_order)])

  # Data prep
  data_filtered <- data %>%
    filter(!is.na(Sample) & !is.na(.data[[y_var]])) %>%
    filter(Sample %in% custom_sample_order) %>%
    mutate(
      Sample = as.character(Sample),
      SampleType = case_when(
        grepl("Base$", Sample, ignore.case = TRUE) ~ "Baseline",
        Sample == "Control" ~ "Control",
        TRUE ~ "Other"
      ),
      SampleType = factor(SampleType, levels = c("Baseline", "Control", "Other")),
      Sample = factor(Sample, levels = custom_sample_order)
    ) %>%
    arrange(Sample)

  # Add gap samples
  gap_samples <- data.frame(
    Sample = c("gap1", "gap2", "gap3", "gap4"),
    SampleType = factor("Gap", levels = c("Baseline", "Control", "Other", "Gap")),
    y_value = 0
  )
  names(gap_samples)[3] <- y_var

  data_filtered <- bind_rows(data_filtered, gap_samples) %>%
    mutate(Sample = factor(Sample, levels = custom_sample_order)) %>%
    arrange(Sample)

  # Determine Y-axis limits
  y_upper_limit <- max(data_filtered[[y_var]], na.rm = TRUE) * 1.3

  # Define vertical line positions
  vline_positions <- c(7.0, 14.0, 20.0, 27.0)
  vline_data <- data.frame(xintercept = vline_positions)

  # Create plot
  p <- ggplot(data_filtered, aes(x = Sample, y = .data[[y_var]], fill = SampleType)) +
    geom_col(alpha = 1, width = 0.65, color = "black", linewidth = 0.1) +
    geom_vline(data = vline_data, aes(xintercept = xintercept),
               linetype = "dashed", color = "black", linewidth = 0.5) +
    labs(title = title, x = "Treatment Timepoint", y = "Mutation Frequency", fill = "Sample Type") +
    scale_fill_manual(
      name = "Sample Type",
      values = c("Baseline" = "#7534b3", "Control" = "#2b8cbe", "Other" = "#fc8d62", "Gap" = "transparent"),
      drop = FALSE
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    ) +
    scale_y_continuous(
      limits = c(0, y_upper_limit),
      expand = expansion(mult = c(0, 0)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    geom_text(aes(label = round(.data[[y_var]], 0)), vjust = -0.5, size = 2.5) +
    scale_x_discrete(labels = function(x) {
      ifelse(grepl("^gap", x), "",
             ifelse(x == "Control", "Control",
                    ifelse(grepl("Base$", x, ignore.case = TRUE), "Baseline",
                           sub(".*BMN|.*BMTu|.*PBMC|.*SkN|.*SkTu", "", x)
                    )
             )
      )
    })

  return(p)
}

# --- Main Execution ---

main <- function() {
  excel_file <- "~/UMCCR/research/projects/CX5461/mutation_counts_summary.xlsx"
  data <- load_and_preprocess_data(excel_file)

  # Generate plots
  plot_types <- list(
    snv = list(var = "snv_count", title = "SNV Counts by Sample"),
    indel = list(var = "indel_count", title = "Indel Counts by Sample"),
    dbs = list(var = "dnv_count", title = "DBS Counts by Sample")
  )

  for (plot_type in names(plot_types)) {
    plot <- create_mutation_plot_single_axis(
      data$merged_data,
      plot_types[[plot_type]]$var,
      plot_types[[plot_type]]$title
    )
    print(plot)
    ggsave(paste0(plot_type, "_counts_sample_single_axis_anno.pdf"), plot = plot, width = 8, height = 4)
  }
}

# Run the main function
main()
