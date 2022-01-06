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



maf <- read.table("~/Documents/UMCCR/data/projects/ctDNA/SFRC01022/sample-only-variants/vcf2maf/0000.maf", sep="\t", as.is=TRUE, header=TRUE, quote = "")

bcbio_fix_df <- maf %>%
  dplyr::select(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  dplyr::rename(chrom=Chromosome, pos=Start_Position, ref=Reference_Allele, alt=Tumor_Seq_Allele2) %>%
  dplyr::mutate(class="bcbio")

write.csv(bcbio_fix_df,"~/Documents/UMCCR/data/projects/ctDNA/SFRC01022/sample-only-variants/vcf2maf/0000.csv", row.names = FALSE, quote = FALSE)
