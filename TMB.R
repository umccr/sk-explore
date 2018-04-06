#Required packages
library(vcfR)
library(dplyr)

#Input file
vcf <-  read.vcfR("../../Data/ensemble-pon-pass.vcf")

# Have a named vector with chr + chr_length
# You can extract the chr_length for chr, and then bin based on that
# In a named vector, you can access each element by name such as chr_lengths[chr_you_want]

chr_lengths <- c(249250621L, 243199373L, 198022430L, 191154276L,
                 180915260L, 171115067L, 159138663L, 146364022L,
                 141213431L, 135534747L, 135006516L, 133851895L,
                 115169878L, 107349540L, 102531392L,  90354753L,
                 81195210L,  78077248L,   59128983L,  63025520L,
                 48129895L,  51304566L,  155270560L,  59373566L, 16569L)
chr_names <-  c(1:22, "X", "Y", "MT")
names(chr_lengths) <- chr_names
bin_size <- 40000000

#Create a vcf with chromosome number and position coloumns, from the priginal vcf. 
#Get chromosome and filter vcf for that chromosome

filter_chr <- function(x, vcf) {
  vcf_chr_pos <- data.frame(chr=getCHROM(vcf), pos=getPOS(vcf), stringsAsFactors = FALSE)
  return(filter(vcf_chr_pos, chr == x))
}

#Prepare bins for the specific chromosome

bin_chr <- function(chr_lengths, chr_name, bin_size) {
  return(seq(from = 0, to = chr_lengths[chr_name], by = bin_size))
}

#Count mutations in a chromosome - cbind combines vector, matrix or data-frame by coloumns 

count_mut_in_chr <- function(chr_pos, bin_vec) {
  x <- cbind(chr_pos$pos, findInterval(chr_pos$pos, bin_vec))
  mean(table(x[, 2]))
}

#Apply same function to all chromosomes

results <- vector("numeric", length = length(chr_names))
for (i in 1:length(chr_names)) {
  results[i] <- count_mut_in_chr(filter_chr(chr_names[i], vcf), bin_chr(chr_lengths, chr_names[i], bin_size))
}








