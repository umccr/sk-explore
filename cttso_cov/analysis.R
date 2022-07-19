library(here)
library(dplyr)
library(tidyverse)
library(stringr)

list_of_files <- list.files(path = "~/UMCCR/data/projects/cttso-qc/MDX220130_L2200592",
                            recursive = TRUE,
                            pattern = "\\.bed.gz$",
                            full.names = TRUE)

df <- list_of_files %>%
  set_names() %>%
  map_df(read.table,
         .id = "sample_id",
         col.names = c("chr", "region_start", "region_end", "mean_coverage" )) %>%
  mutate(sample_id = tools::file_path_sans_ext(basename(sample_id))) %>%
  mutate(sample_id = gsub(".regions.bed", "", sample_id)) %>%
  dplyr::select(sample_id, region_start, region_end, mean_coverage)


# display input data in a table
datatable(df, rownames = FALSE, filter="top", options = list(pageLength = 5, scrollX=T), caption = "Mean coverage of TERT promoter region and region of interest")

#read input
coverage_beds <-list.files("~/UMCCR/data/projects/cttso-qc/mosdepth_output/",
           full.names = T)

data <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(data) = c("sample", "mean_coverage_ratio")
file <- data.frame()

for (i in 1:length(coverage_beds)){
  print(coverage_beds[i])
  file <- read.table(coverage_beds[i],col.names = c("chr", "region_start", "region_end", "mean_coverage" )) %>%
    mutate(sample_name = basename(coverage_beds[i])) %>%
    mutate(sample_name = gsub(".regions.bed.gz", "", sample_name)) %>%
    mutate(ratio = mean_coverage[1] / mean_coverage[2]) %>%
    distinct(.keep_all = FALSE)


  sample_name = unique(file$sample_name)
  ratio_mean_cov = unique(file$ratio)

  data[nrow(data) + 1, ] <- c(sample_name, round(ratio_mean_cov, 2))
}

# sort ratio
data <- data[order(data$mean_coverage_ratio),]
# convert ratio values to numeric and omit nan values
data$mean_coverage_ratio <- as.numeric(as.character(data$mean_coverage_ratio))
data <- na.omit(data)
# plot results https://r-graph-gallery.com/267-reorder-a-variable-in-ggplot2.html
data %>%
  arrange(mean_coverage_ratio) %>%
  mutate(sample=factor(sample, levels=sample)) %>% # This trick update the factor levels
  ggplot( aes(x=mean_coverage_ratio, y=sample)) +
  geom_segment( aes(xend=0, yend=sample)) +
  geom_point( size=4, color="orange") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("")
