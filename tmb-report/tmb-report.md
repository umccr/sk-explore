---
title: "Tumour Mutational Burden (TMB)"
author: "Sehrish Kanwal"
date: "Thu 2018-Nov-15"
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
#library(vcfR)
library(knitr)
library(kableExtra)
```


```r
#vcf <-  read.vcfR("~/Documents/UMCCR/data/vcfs/ensemble-with_chr_prefix.vcf", verbose =  FALSE) 
vcf <- "~/Documents/UMCCR/data/vcfs/ensemble-pon-pass.vcf"
#bcftools <- "bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n'"

#lets just read the annotaions from vcf file
bcftools <- "bcftools query -f '%INFO/ANN\n'"

#vcf <- system(paste(bcftools, vcf_file), intern = TRUE) %>%
  #tibble::tibble(all_cols = .) %>%
  #tidyr::separate(col = .data$all_cols,
                  #into = c("CHROM", "POS", "ANN"),
                  #sep = "\t", convert = TRUE)

#run bcftools command on the vcf and convert output to dataframe
vcf_ann <- system(paste(bcftools, vcf), intern = TRUE) %>%
  tibble::tibble(ANN = .)

#Calculating mutations per megabase
#fix <- getFIX(vcf)
vcf_number_rows <- nrow(vcf_ann)
mutations_megabase <- round(vcf_number_rows/3200, digits = 2)

#Summarizing annotations for variants in the vcf
#ann <- vcfR::extract.info(vcf, "ANN")
#ann <- vcf$ANN
region_ann <- sapply(vcf_ann$ANN, function(x){
  y <- strsplit(x, "\\|")[[1]][2]
})
variant_annotation <- unname(region_ann)

#Creating a nice table output for annotations summary
region_ann_df <- data.frame(table(variant_annotation))
kable(region_ann_df, caption = "Table summarizing all annotations in the vcf and the total number of variants suppporting these annotations")  %>%
  kable_styling(font_size = 12, "striped", "bordered")
```

<table class="table table-striped" style="font-size: 12px; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">Table summarizing all annotations in the vcf and the total number of variants suppporting these annotations</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> variant_annotation </th>
   <th style="text-align:right;"> Freq </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 3_prime_UTR_variant </td>
   <td style="text-align:right;"> 104 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5_prime_UTR_premature_start_codon_gain_variant </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5_prime_UTR_variant </td>
   <td style="text-align:right;"> 36 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> conservative_inframe_deletion </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> downstream_gene_variant </td>
   <td style="text-align:right;"> 1311 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> frameshift_variant </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> frameshift_variant&amp;splice_donor_variant&amp;splice_region_variant&amp;intron_variant </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> initiator_codon_variant </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> intergenic_region </td>
   <td style="text-align:right;"> 10332 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> intron_variant </td>
   <td style="text-align:right;"> 6642 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> missense_variant </td>
   <td style="text-align:right;"> 101 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> missense_variant&amp;splice_region_variant </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> non_coding_transcript_exon_variant </td>
   <td style="text-align:right;"> 53 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> splice_acceptor_variant&amp;intron_variant </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> splice_donor_variant&amp;intron_variant </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> splice_region_variant </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> splice_region_variant&amp;intron_variant </td>
   <td style="text-align:right;"> 17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> splice_region_variant&amp;synonymous_variant </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> start_lost </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stop_gained </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stop_lost </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> synonymous_variant </td>
   <td style="text-align:right;"> 47 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> upstream_gene_variant </td>
   <td style="text-align:right;"> 1679 </td>
  </tr>
</tbody>
</table>

```r
  
#Calculating mutations per megabase in coding region
coding_variants = 0
coding_variants <- region_ann %in% c("frameshift_variant", "missense_variant", "missense_variant&splice_region_variant")
coding_variants <- table(coding_variants)
mutations_megabase_coding <- round(as.vector(coding_variants[2])/40, digits = 2) 
#40MB is the estimated size of coding region in human genome - as used by PCGR as well. 
#We can use 36MB if we go with exact calculations, as only 1.2% of the total genome is considered coding. 
#total genome * percent protein coding = 3,000,000,000 * 0.012 = 36,000,000 ~36MB

#Displaying results in a table
region <- c("Wholegenome", "Coding")
total_mutations <- c(vcf_number_rows, as.vector(coding_variants[2]))
mutations_mb <- c(mutations_megabase, mutations_megabase_coding)

result_display <- data.frame(region, total_mutations, mutations_mb) 
#kable(result_display,  caption = "Table summarizing somatic burden result")

kable(result_display, caption = "Table summarizing somatic burden result") %>%
  kable_styling(font_size = 12, "striped", "bordered")
```

<table class="table table-striped" style="font-size: 12px; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">Table summarizing somatic burden result</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> region </th>
   <th style="text-align:right;"> total_mutations </th>
   <th style="text-align:right;"> mutations_mb </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Wholegenome </td>
   <td style="text-align:right;"> 20361 </td>
   <td style="text-align:right;"> 6.36 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Coding </td>
   <td style="text-align:right;"> 115 </td>
   <td style="text-align:right;"> 2.88 </td>
  </tr>
</tbody>
</table>

* The _total number of mutations_ in the vcf are **20361** and 
* Number of mutations per megabase are **6.36**.
* The _total number of mutations in the coding region_ are **115**
* Number of mutations per megabase in the coding region are **2.88**

## Addendum


```r
devtools::session_info()
## ─ Session info ──────────────────────────────────────────────────────────
##  setting  value                       
##  version  R version 3.5.0 (2018-04-23)
##  os       macOS  10.14.1              
##  system   x86_64, darwin15.6.0        
##  ui       X11                         
##  language (EN)                        
##  collate  en_AU.UTF-8                 
##  ctype    en_AU.UTF-8                 
##  tz       Australia/Melbourne         
##  date     2018-11-15                  
## 
## ─ Packages ──────────────────────────────────────────────────────────────
##  package     * version date       lib source        
##  assertthat    0.2.0   2017-04-11 [1] CRAN (R 3.5.0)
##  backports     1.1.2   2017-12-13 [1] CRAN (R 3.5.0)
##  base64enc     0.1-3   2015-07-28 [1] CRAN (R 3.5.0)
##  callr         3.0.0   2018-08-24 [1] CRAN (R 3.5.0)
##  cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.0)
##  colorspace    1.3-2   2016-12-14 [1] CRAN (R 3.5.0)
##  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)
##  desc          1.2.0   2018-05-01 [1] CRAN (R 3.5.0)
##  devtools      2.0.1   2018-10-26 [1] CRAN (R 3.5.0)
##  digest        0.6.18  2018-10-10 [1] CRAN (R 3.5.0)
##  evaluate      0.12    2018-10-09 [1] CRAN (R 3.5.0)
##  fs            1.2.6   2018-08-23 [1] CRAN (R 3.5.0)
##  glue          1.3.0   2018-07-17 [1] CRAN (R 3.5.0)
##  highr         0.7     2018-06-09 [1] CRAN (R 3.5.0)
##  hms           0.4.2   2018-03-10 [1] CRAN (R 3.5.0)
##  htmltools     0.3.6   2017-04-28 [1] CRAN (R 3.5.0)
##  httr          1.3.1   2017-08-20 [1] CRAN (R 3.5.0)
##  kableExtra  * 0.9.0   2018-05-21 [1] CRAN (R 3.5.0)
##  knitr       * 1.20    2018-02-20 [1] CRAN (R 3.5.0)
##  magrittr      1.5     2014-11-22 [1] CRAN (R 3.5.0)
##  memoise       1.1.0   2017-04-21 [1] CRAN (R 3.5.0)
##  munsell       0.5.0   2018-06-12 [1] CRAN (R 3.5.0)
##  pillar        1.3.0   2018-07-14 [1] CRAN (R 3.5.0)
##  pkgbuild      1.0.2   2018-10-16 [1] CRAN (R 3.5.0)
##  pkgconfig     2.0.2   2018-08-16 [1] CRAN (R 3.5.0)
##  pkgload       1.0.2   2018-10-29 [1] CRAN (R 3.5.0)
##  prettyunits   1.0.2   2015-07-13 [1] CRAN (R 3.5.0)
##  processx      3.2.0   2018-08-16 [1] CRAN (R 3.5.0)
##  ps            1.2.1   2018-11-06 [1] CRAN (R 3.5.0)
##  R6            2.3.0   2018-10-04 [1] CRAN (R 3.5.0)
##  Rcpp          1.0.0   2018-11-07 [1] CRAN (R 3.5.0)
##  readr         1.1.1   2017-05-16 [1] CRAN (R 3.5.0)
##  remotes       2.0.2   2018-10-30 [1] CRAN (R 3.5.0)
##  rlang         0.3.0.1 2018-10-25 [1] CRAN (R 3.5.0)
##  rmarkdown     1.10    2018-06-11 [1] CRAN (R 3.5.0)
##  rprojroot     1.3-2   2018-01-03 [1] CRAN (R 3.5.0)
##  rstudioapi    0.8     2018-10-02 [1] CRAN (R 3.5.0)
##  rvest         0.3.2   2016-06-17 [1] CRAN (R 3.5.0)
##  scales        1.0.0   2018-08-09 [1] CRAN (R 3.5.0)
##  sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.5.0)
##  stringi       1.2.4   2018-07-20 [1] CRAN (R 3.5.0)
##  stringr       1.3.1   2018-05-10 [1] CRAN (R 3.5.0)
##  testthat      2.0.1   2018-10-13 [1] CRAN (R 3.5.0)
##  tibble        1.4.2   2018-01-22 [1] CRAN (R 3.5.0)
##  usethis       1.4.0   2018-08-14 [1] CRAN (R 3.5.0)
##  viridisLite   0.3.0   2018-02-01 [1] CRAN (R 3.5.0)
##  withr         2.1.2   2018-03-15 [1] CRAN (R 3.5.0)
##  xml2          1.2.0   2018-01-24 [1] CRAN (R 3.5.0)
##  yaml          2.2.0   2018-07-25 [1] CRAN (R 3.5.0)
## 
## [1] /Library/Frameworks/R.framework/Versions/3.5/Resources/library
```
