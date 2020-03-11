library(vcfR)
library(dplyr)
library(here)

#read vcf input
bcbio_2 <- read.vcfR(here("/UMCCR/data/projects/varlap/0000.vcf"), verbose = FALSE)

#convert vcf object to tidy object - the function returns a list with three components: fix, gt, and meta 
bcbio_tidy <- vcfR2tidy(bcbio_2)

#extract fix component, select, rename and mutate columns of interest
bcbio_fix_df <- bcbio_tidy$fix %>%
  dplyr::select(CHROM, POS, REF, ALT) %>%
  dplyr::rename(chrom=CHROM, pos=POS, ref=REF, alt=ALT) %>%
  dplyr::mutate(class="bcbio")

#save output
write.csv(bcbio_fix_df, here("/UMCCR/data/projects/varlap/bcbio_fix.csv"), row.names = FALSE, quote = FALSE)
