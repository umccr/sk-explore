ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(dplyr)
library(patchwork)
library(kableExtra)
library(deconstructSigs)

# read data
sample.mut.ref <- read.table("./data/sample.mut.filter.ref.tsv",
                             header = FALSE, 
                             sep = "\t",
                             stringsAsFactors = FALSE)
sample.mut.ref$V1 <- 1

# convert to deconstructSigs input
sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref, 
                                sample.id = "V1", 
                                chr = "V2", 
                                pos = "V3", 
                                ref = "V4", 
                                alt = "V5",
                                bsg =BSgenome.Hsapiens.UCSC.hg38)
col_names <- colnames(sigs.input)
third_letter <- gsub("^.\\[([A-Z])>[A-Z]\\].$", "\\1", col_names)
sorted_indices <- order(third_letter)
sigs.input <- sigs.input[, sorted_indices]
head(sigs.input)

# read reference
cosmic_v33 <- read.table(
  "./reference/COSMIC_v3.3.1_SBS_GRCh38.txt",
  header = TRUE, 
  sep = "\t",
  row.names = 1 
)
cosmic_v33_t <- as.data.frame(t(cosmic_v33))
col_names <- colnames(cosmic_v33_t)
mutation_types <- gsub("^.\\[([A-Z]>[A-Z])\\].$", "\\1", col_names)
sorted_indices <- order(mutation_types)
cosmic_v33_t<- cosmic_v33_t[, sorted_indices]
head(cosmic_v33_t)

# call signatures
sample_1 = deconstructSigs::whichSignatures(tumor.ref = sigs.input, 
                           signatures.ref = cosmic_v33_t, 
                           sample.id = 1, 
                           contexts.needed = TRUE,
                           signature.cutoff = 0.1,
                           tri.counts.method = "genome")

# plot output
deconstructSigs::plotSignatures(sample_1, sig.type = "SBS")
