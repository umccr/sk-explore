---
title: "Tumour Mutational Burden (TMB)"
author: "Sehrish Kanwal"
date: "Thu 2018-Aug-30"
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
library(vcfR)
library(knitr)
```


```r
vcf <-  read.vcfR("~/Documents/UMCCR/data/vcfs/MH17T001P004-somatic-ensemble-pon_hardfiltered.vcf.gz", verbose =  FALSE) 

#Calculating mutations per megabase
fix <- getFIX(vcf)
vcf_number_rows <- nrow(fix)
mutations_megabase <- round(vcf_number_rows/3200, digits = 2)

#Summarizing annotations for variants in the vcf
ann <- vcfR::extract.info(vcf, "ANN")
region_ann <- sapply(ann, function(x){
  y <- strsplit(x, "\\|")[[1]][2]
})
variant_annotation <- unname(region_ann)

#Creating a nice table output for annotations summary
region_ann_df <- data.frame(table(variant_annotation))
kable(region_ann_df, caption = "Table summarizing all annotations in the vcf and the total number of variants suppporting these annotations", format = "markdown")
```



|variant_annotation                                       | Freq|
|:--------------------------------------------------------|----:|
|3_prime_UTR_variant                                      |    7|
|5_prime_UTR_variant                                      |    1|
|downstream_gene_variant                                  |   67|
|frameshift_variant                                       |    1|
|intergenic_region                                        |  222|
|intron_variant                                           |  203|
|missense_variant                                         |    2|
|sequence_feature                                         |    3|
|splice_region_variant&intron_variant                     |    3|
|splice_region_variant&non_coding_transcript_exon_variant |    1|
|TF_binding_site_variant                                  |    1|
|upstream_gene_variant                                    |  107|

```r

#Calculating mutations per megabase in coding region
coding_variants = 0
coding_variants <- region_ann %in% c("frameshift_variant", "missense_variant", "missense_variant&splice_region_variant")
coding_variants <- table(coding_variants)
mutations_megabase_coding <- round(as.vector(coding_variants[2])/40, digits = 2) 
#40MB is the estimated size of coding region in human genome - as used by PCGR as well. 
#We can use 36MB if we go with exact calculations, as only 1.2% of the total genome is considered coding. 
#total genome * percent protein coding = 3,000,000,000 * 0.012 = 36,000,000 ~36MB
```

* The _total number of mutations_ in the vcf are **618** and 
* Number of mutations per megabase are **0.19**.
* The _total number of mutations in the coding region_ are **3**
* Number of mutations per megabase in the coding region are **0.08**







