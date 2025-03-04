library(dplyr)
library(tools)
library(stringr) 

# input data
setwd("~/UMCCR/research/projects/mesogenome/data/purple-cnv")
input_folder <- "~/UMCCR/research/projects/mesogenome/data/purple-cnv"
output_file <- "~/UMCCR/research/projects/mesogenome/data/purple-cnv/combined_segments.seg" 
file_list <- list.files(input_folder, pattern = "*.tsv", full.names = TRUE)

# combine cnv outputs to run with gistic
combined_data <- data.frame()
for (file in file_list) {
  sample_data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sample_id <- tools::file_path_sans_ext(basename(file))
  
  sample_data <- sample_data %>%
    mutate(Sample = sample_id,
           Chrom = chromosome,
           Start = start,
           Stop = end,
           Marker = if_else(
             str_detect(bafCount, "\\([0-9]+\\)"),
             as.numeric(str_extract(bafCount, "(?<=\\()[0-9]+(?=\\))")),
             as.numeric(bafCount)  
           ),
           Seg.CN = log2(copyNumber) - 1) %>%  
    dplyr::select(Sample, Chrom, Start, Stop, Marker, Seg.CN)
  
  combined_data <- bind_rows(combined_data, sample_data)
}
combined_data$Sample <- gsub("\\.purple.cnv.somatic", "", combined_data$Sample)

# write output
write.table(combined_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
unique(combined_data$Sample)