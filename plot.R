#Required packages
library(vcfR)
library(dplyr)

#Input file

vcf <-  read.vcfR("../../Data/p013-tumor-ensemble-annotated.vcf.gz")
chrom <- vcfR::getCHROM(vcf)
x <- table(chrom)

# How can I sort not alphanumerically
# Use gtools::mixedorder to get the ordered indices (or mixedsort to sort them)
sorted_ind <- gtools::mixedorder(names(x))
x <- x[sorted_ind]
chr_fac <- factor(names(x), levels = gtools::mixedsort(names(x)))
df <- data.frame(chr = chr_fac, mut_num = c(unname(x)))
# target_size_mb = 40
# 
# sample_calls_tmb <- sample_calls %>% dplyr::filter(stringr::str_detect(CONSEQUENCE,"^(stop_gained|stop_lost|start_lost|frameshift_variant|missense_variant|inframe_deletion|inframe_insertion|synonymous)"))
# 
# tmb_estimate <- round(as.numeric(nrow(sample_calls_tmb)/ target_size_mb), digits = 2)
                      


# chrom | n mut

df

ggplot(df, aes(x = chr, y = mut_num)) +
  geom_point(colour = "purple") +
  theme_bw()
plot(rnorm(100))        
