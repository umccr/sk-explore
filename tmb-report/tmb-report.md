---
title: "Tumour Mutational Burden (TMB)"
author: "Sehrish Kanwal"
date: "Fri 2018-Oct-26"
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
library(kableExtra)
```


```r
vcf <-  read.vcfR("~/Documents/UMCCR/data/vcfs/ensemble-with_chr_prefix.vcf", verbose =  FALSE) 

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
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> downstream_gene_variant </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> intergenic_region </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> intron_variant </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> missense_variant </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> upstream_gene_variant </td>
   <td style="text-align:right;"> 7 </td>
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
   <td style="text-align:right;"> 38 </td>
   <td style="text-align:right;"> 0.01 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Coding </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.10 </td>
  </tr>
</tbody>
</table>

* The _total number of mutations_ in the vcf are **38** and 
* Number of mutations per megabase are **0.01**.
* The _total number of mutations in the coding region_ are **4**
* Number of mutations per megabase in the coding region are **0.1**

## Addendum


```r
devtools::session_info()
## Session info -------------------------------------------------------------
##  setting  value                       
##  version  R version 3.5.0 (2018-04-23)
##  system   x86_64, darwin15.6.0        
##  ui       X11                         
##  language (EN)                        
##  collate  en_AU.UTF-8                 
##  tz       Australia/Melbourne         
##  date     2018-10-26
## Packages -----------------------------------------------------------------
##  package     * version date       source        
##  ape           5.2     2018-09-24 CRAN (R 3.5.0)
##  backports     1.1.2   2017-12-13 CRAN (R 3.5.0)
##  base        * 3.5.0   2018-04-24 local         
##  cluster       2.0.7-1 2018-04-13 CRAN (R 3.5.0)
##  colorspace    1.3-2   2016-12-14 CRAN (R 3.5.0)
##  compiler      3.5.0   2018-04-24 local         
##  crayon        1.3.4   2017-09-16 CRAN (R 3.5.0)
##  datasets    * 3.5.0   2018-04-24 local         
##  devtools      1.13.6  2018-06-27 CRAN (R 3.5.0)
##  digest        0.6.17  2018-09-12 CRAN (R 3.5.0)
##  evaluate      0.11    2018-07-17 CRAN (R 3.5.0)
##  graphics    * 3.5.0   2018-04-24 local         
##  grDevices   * 3.5.0   2018-04-24 local         
##  grid          3.5.0   2018-04-24 local         
##  highr         0.7     2018-06-09 CRAN (R 3.5.0)
##  hms           0.4.2   2018-03-10 CRAN (R 3.5.0)
##  htmltools     0.3.6   2017-04-28 CRAN (R 3.5.0)
##  httr          1.3.1   2017-08-20 CRAN (R 3.5.0)
##  kableExtra  * 0.9.0   2018-05-21 CRAN (R 3.5.0)
##  knitr       * 1.20    2018-02-20 CRAN (R 3.5.0)
##  lattice       0.20-35 2017-03-25 CRAN (R 3.5.0)
##  magrittr      1.5     2014-11-22 CRAN (R 3.5.0)
##  MASS          7.3-50  2018-04-30 CRAN (R 3.5.0)
##  Matrix        1.2-14  2018-04-13 CRAN (R 3.5.0)
##  memoise       1.1.0   2017-04-21 CRAN (R 3.5.0)
##  memuse        4.0-0   2017-11-10 CRAN (R 3.5.0)
##  methods     * 3.5.0   2018-04-24 local         
##  mgcv          1.8-24  2018-06-18 CRAN (R 3.5.0)
##  munsell       0.5.0   2018-06-12 CRAN (R 3.5.0)
##  nlme          3.1-137 2018-04-07 CRAN (R 3.5.0)
##  parallel      3.5.0   2018-04-24 local         
##  permute       0.9-4   2016-09-09 CRAN (R 3.5.0)
##  pillar        1.3.0   2018-07-14 CRAN (R 3.5.0)
##  pinfsc50      1.1.0   2016-12-02 CRAN (R 3.5.0)
##  pkgconfig     2.0.2   2018-08-16 CRAN (R 3.5.0)
##  R6            2.3.0   2018-10-04 CRAN (R 3.5.0)
##  Rcpp          0.12.19 2018-10-01 CRAN (R 3.5.0)
##  readr         1.1.1   2017-05-16 CRAN (R 3.5.0)
##  rlang         0.2.2   2018-08-16 CRAN (R 3.5.0)
##  rmarkdown     1.10    2018-06-11 CRAN (R 3.5.0)
##  rprojroot     1.3-2   2018-01-03 CRAN (R 3.5.0)
##  rstudioapi    0.8     2018-10-02 CRAN (R 3.5.0)
##  rvest         0.3.2   2016-06-17 CRAN (R 3.5.0)
##  scales        1.0.0   2018-08-09 CRAN (R 3.5.0)
##  stats       * 3.5.0   2018-04-24 local         
##  stringi       1.2.4   2018-07-20 CRAN (R 3.5.0)
##  stringr       1.3.1   2018-05-10 CRAN (R 3.5.0)
##  tibble        1.4.2   2018-01-22 CRAN (R 3.5.0)
##  tools         3.5.0   2018-04-24 local         
##  utils       * 3.5.0   2018-04-24 local         
##  vcfR        * 1.8.0   2018-04-17 CRAN (R 3.5.0)
##  vegan         2.5-2   2018-05-17 CRAN (R 3.5.0)
##  viridisLite   0.3.0   2018-02-01 CRAN (R 3.5.0)
##  withr         2.1.2   2018-03-15 CRAN (R 3.5.0)
##  xml2          1.2.0   2018-01-24 CRAN (R 3.5.0)
##  yaml          2.2.0   2018-07-25 CRAN (R 3.5.0)
```
