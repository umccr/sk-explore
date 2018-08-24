---
title: "Tumour Mutational Burden (TMB)"
author: "Sehrish Kanwal"
date: "Fri 2018-Aug-24"
output:
  html_document:
    keep_md: true
editor_options:
  chunk_output_type: console
---



Tumour mutational burden (TMB) "measures the quantity of mutations found in a tumor" (https://www.focr.org/tmb). Also, it is defined as "a quantitative measure of the total number of mutations per coding area of a tumor genome." (https://www.genengnews.com/gen-exclusives/the-promise-of-tumor-mutational-burden-for-cancer-immunotherapy/77900833).
This type of biomarker is currently under study to evaluate whether it may help predict the likelihood a patient's response to immunotherapy in a range of advanced cancers. 

Tumors that have higher levels of TMB are believed to express more neoantigens – a type of cancer-specific antigen – that may allow for a more robust immune response and therefore a more durable response to immunotherapy.

Required packages


```r
library(vcfR, warn.conflicts=F, quietly=T)
library(dplyr, warn.conflicts=F, quietly=T)
## Warning: package 'dplyr' was built under R version 3.5.1
library(ggplot2, warn.conflicts=F, quietly=T)
library(reticulate, warn.conflicts=F, quietly=T)
```


```r
vcf <-  read.vcfR("~/Documents/UMCCR/data/ensemble-pon-pass.vcf") 
## Scanning file to determine attributes.
## File attributes:
##   meta lines: 343
##   header_line: 344
##   variant count: 20361
##   column count: 11
## Meta line 343 read in.
## All meta lines processed.
## gt matrix initialized.
## Character matrix gt created.
##   Character matrix gt rows: 20361
##   Character matrix gt cols: 11
##   skip: 0
##   nrows: 20361
##   row_num: 0
## Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant 10000Processed variant 11000Processed variant 12000Processed variant 13000Processed variant 14000Processed variant 15000Processed variant 16000Processed variant 17000Processed variant 18000Processed variant 19000Processed variant 20000Processed variant: 20361
## All variants processed

#Calculating mutations per megabase
fix <- getFIX(vcf)
vcf_number_rows <- nrow(fix)
mutations_megabase <- ceiling(vcf_number_rows/3200)

#Summarizing annotations for variants in the vcf
ann <- vcfR::extract.info(vcf, "ANN")
region_ann <- sapply(ann, function(x){
  y <- strsplit(x, "\\|")[[1]][2]
})
region_ann <- unname(region_ann)
table(region_ann)
## region_ann
##                                                          3_prime_UTR_variant 
##                                                                          104 
##                               5_prime_UTR_premature_start_codon_gain_variant 
##                                                                            8 
##                                                          5_prime_UTR_variant 
##                                                                           36 
##                                                conservative_inframe_deletion 
##                                                                            1 
##                                                      downstream_gene_variant 
##                                                                         1311 
##                                                           frameshift_variant 
##                                                                            9 
## frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant 
##                                                                            1 
##                                                      initiator_codon_variant 
##                                                                            1 
##                                                            intergenic_region 
##                                                                        10332 
##                                                               intron_variant 
##                                                                         6642 
##                                                             missense_variant 
##                                                                          101 
##                                       missense_variant&splice_region_variant 
##                                                                            5 
##                                           non_coding_transcript_exon_variant 
##                                                                           53 
##                                       splice_acceptor_variant&intron_variant 
##                                                                            2 
##                                          splice_donor_variant&intron_variant 
##                                                                            4 
##                                                        splice_region_variant 
##                                                                            1 
##                                         splice_region_variant&intron_variant 
##                                                                           17 
##                                     splice_region_variant&synonymous_variant 
##                                                                            1 
##                                                                   start_lost 
##                                                                            1 
##                                                                  stop_gained 
##                                                                            4 
##                                                                    stop_lost 
##                                                                            1 
##                                                           synonymous_variant 
##                                                                           47 
##                                                        upstream_gene_variant 
##                                                                         1679

#Calculating mutations per megabase in coding region
coding_variants = 0
coding_variants <- ifelse(region_ann %in% c("frameshift_variant", "missense_variant", "missense_variant&splice_region_variant"),
       TRUE,
       FALSE)
coding_variants <- table(coding_variants)
mutations_megabase_coding <- ceiling(as.vector(coding_variants[2])/40) #40MB is the estimated size of coding region in human genome - as used by PCGR as well. We can use 36MB if we go with exact calculations, as only 1.2% of the total genome is considered coding (total genome * percent protein coding = 3,000,000,000 * 0.012 = 36,000,000 ~36MB)
```

* The _total number of mutations_ in the vcf are **20361** and 
* Number of mutations per megabase are **7**.
* The _total number of mutations in the coding region_ are **115**
* Number of mutations per megabase in the coding region are **3**







