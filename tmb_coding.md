---
title: "Tumour Mutational Burden (TMB)"
author: "Sehrish Kanwal"
date: "Tue 2018-Aug-21"
output: 
  html_document: 
    keep_md: yes
editor_options: 
  chunk_output_type: console
---


Tumour mutational burden (TMB) "measures the quantity of mutations found in a tumor" (https://www.focr.org/tmb). Also, it is defined as "a quantitative measure of the total number of mutations per coding area of a tumor genome." (https://www.genengnews.com/gen-exclusives/the-promise-of-tumor-mutational-burden-for-cancer-immunotherapy/77900833).
This type of biomarker is currently under study to evaluate whether it may help predict the likelihood a patient's response to immunotherapy in a range of advanced cancers.

As multiple definition indicate, currently, there is a lack of standardization for TMB calculation and reporting. However, it is imperative to use this information to further inform clinically-meaningful results, which will support informed decision-making for patients. Most importantly, it provides a quantitative measure that can be used to better inform treatment decisions.

The report describes the work done for analysing and understanding TMB on one of the ensemble somatic variant calls (`/data/cephfs/punim0010/data/Results/Patients/CUP_SC932/final/umccrised/cup_tissue/somatic/ensemble-pon-pass.vcf.gz`). Many thanks to Peter for providing specific pointers and useful explanation on how to tackle various issues in R. The current script is an attempt to filter variants specifically in the coding region.

Required packages


```r
library(vcfR, warn.conflicts=F, quietly=T)
library(dplyr, warn.conflicts=F, quietly=T)
## Warning: package 'dplyr' was built under R version 3.5.1
library(ggplot2, warn.conflicts=F, quietly=T)
library(reticulate, warn.conflicts=F, quietly=T)
```

Reading in the the input file


```r
vcf <-  read.vcfR("~/Documents/UMCCR/data/ensemble-pon-pass.vcf") # 20,361 variants in total
#head(vcf)
```

Have a named vector with chr + chr_length.
You can extract the chr_length for chr, and then bin based on that.
In a named vector, you can access each element by name such as chr_lengths[chr_you_want].


```r
chr_lengths <- c(249250621L, 243199373L, 198022430L, 191154276L,
                 180915260L, 171115067L, 159138663L, 146364022L,
                 141213431L, 135534747L, 135006516L, 133851895L,
                 115169878L, 107349540L, 102531392L,  90354753L,
                 81195210L,  78077248L,   59128983L,  63025520L,
                 48129895L,  51304566L,  155270560L,  59373566L)
chr_names <- c(1:22, "X", "Y")
names(chr_lengths) <- chr_names
```

Declare bin size/target size (in bases) - it can be variable. This size will be used to create bins for all the chromosomes in the vcf and count how many mutations are in each bin.
*NOTE:*
Binning the chromosomes has no biological significance/basis. Just trying to play with different R fucnctions. 


```r
bin_size <- 40000000
```



```r
chrom <- vcfR::getCHROM(vcf)
x <- table(chrom)
x
## chrom
##    1   10   11   12   13   14   15   16   17   18   19    2   20   21   22 
## 1142  828  875  712 1046  689  448  562  274  585  290 1434  284  368  130 
##    3    4    5    6    7    8    9    X    Y 
## 1155 1543 1502 1014 1130 1187  706 2455    2
# How can I sort not alphanumerically
# Use gtools::mixedorder to get the ordered indices (or mixedsort to sort them)
sorted_ind <- gtools::mixedorder(names(x))
sorted_ind
##  [1]  1 12 16 17 18 19 20 21 22  2  3  4  5  6  7  8  9 10 11 13 14 15 23
## [24] 24
x <- x[sorted_ind]
x
## chrom
##    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
## 1142 1434 1155 1543 1502 1014 1130 1187  706  828  875  712 1046  689  448 
##   16   17   18   19   20   21   22    X    Y 
##  562  274  585  290  284  368  130 2455    2
chr_fac <- factor(names(x), levels = gtools::mixedsort(names(x)))
chr_fac
##  [1] 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 X 
## [24] Y 
## 24 Levels: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ... Y
```

**Calculating mutations per megabase**


```r
fix <- getFIX(vcf)
vcf_number_rows <- nrow(fix)
mutations_megabase <- ceiling(vcf_number_rows/3200) 
```

* The _total number of mutations_ in the vcf are **20361** and 
* Number of mutations per megabase_ are **7**.

**Counting mutations for different annotations**


```r
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
```

Get chromosome and filter vcf for that chromosome


```r
# takes a vcf and a chromosome
# outputs df with pos + chr for that chromosome
filter_chr <- function(x, vcf) {
  vcf_chr_pos_ann <- data.frame(chr=getCHROM(vcf), pos=getPOS(vcf), ann=vcfR::extract.info(vcf, "ANN"), stringsAsFactors = FALSE)
  return(filter(vcf_chr_pos_ann, chr == x))
}
```

Prepare bins for the specific chromosome


```r
bin_chr <- function(chr_lengths, chr_name, bin_size) {
  return(seq(from = 0, to = chr_lengths[chr_name], by = bin_size))
}
```

Count mutations in a chromosome, using specific annotations. For now focussing on using multiple annotations that are more or less related to coding sequence, as the vcf I am testing this code on has no variant with the annotation "coding_sequence_variant", i.e. a sequence variant that changes coding sequence (http://sequenceontology.org/browser/current_svn/term/SO:0001580).


```r
count_mut_per_bin <- function(chr_pos_ann, bin_vec) {
  #extract only annotations form vcf
  region_ann <- sapply(chr_pos_ann$ann, function(x){
     y <- strsplit(x, "\\|")[[1]][2]
     })
  region_ann <- unname(region_ann)
  dfs_list <- list()
  for (i in 1:length(region_ann)){
    #Filter annotations of interest 
    if (region_ann[i] == "frameshift_variant" | region_ann[i] == "missense_variant" | 
        region_ann[i] == "missense_variant&splice_region_variant"){
        df <- data.frame(chr_pos_ann[i,2], findInterval(chr_pos_ann[i,2], bin_vec))
        dfs_list[[i]] <- df
    }
  }
  
#check if the annotation used for filtering exists in the vcf - i.e. there are rows in the combined df
if( length(dfs_list) != 0){
  #Bind rows together
  combined <- bind_rows(dfs_list)
  tab <- table(combined[, 2])
  df <- data.frame(bin_num = names(tab),
                   num_mut_in_bin = c(unname(tab)),
                   stringsAsFactors = FALSE)
} else {
  df <- data.frame(bin_num = 0,
                   num_mut_in_bin = 0,
                   stringsAsFactors = FALSE)
}
df
}
```

Apply same function to all chromosomes. This gives a numerical vector containing the number of variants for each chromosome 


```r
results <- vector("numeric", length = length(chr_names))
results3 <- vector("numeric", length = length(chr_names))
for (i in 1:length(chr_names)) {
  #check if there exists any mutation with the specified annotation, in the vcf. 
  if(count_mut_per_bin(filter_chr(chr_names[i], vcf), bin_chr(chr_lengths, chr_names[i], bin_size))[1, 1] != 0){
    res <- count_mut_per_bin(filter_chr(chr_names[i], vcf), bin_chr(chr_lengths, chr_names[i], bin_size))
    results[i] <- mean(res[, 2])
    results3[i] <- sum(res[, 2])
  } else{
    results[i] <- 0
    results3[i] <- 0
  }
}
print(results)
##  [1] 2.800000 1.750000 1.500000 1.000000 1.333333 1.666667 2.000000
##  [8] 1.333333 1.333333 1.666667 1.666667 2.000000 2.000000 1.500000
## [15] 1.000000 2.333333 3.000000 1.000000 5.000000 1.000000 1.000000
## [22] 1.500000 2.750000 0.000000
print(results3)
##  [1] 14  7  6  2  4  5  8  4  4  5  5  6  4  3  1  7  3  1 10  1  1  3 11
## [24]  0
```


* Calculate number of mutations per bin for all chromosomes:


```r
results2 <- vector("list", length = length(chr_names))
for (i in 1:length(chr_names)) {
  #check if there exists any mutation with the specified annotation, in the vcf. 
  if(count_mut_per_bin(filter_chr(chr_names[i], vcf), bin_chr(chr_lengths, chr_names[i], bin_size))[1, 1] != 0){
    results2[[i]] <- count_mut_per_bin(filter_chr(chr_names[i], vcf), bin_chr(chr_lengths, chr_names[i], bin_size))
  }
}
```

* Turn a list that contains X data.frames to a single data.frame


```r
#check if results2 list does not have just all NULLs
if (all(sapply(results2, function(x) is.null(x))) != TRUE){
  names(results2) <- paste("chr", chr_names, sep = "")
  str(results2)
  #binding1 <- do.call("rbind", results2) # this sucks
  binding2 <- dplyr::bind_rows(results2, .id = "chromosome")
}
```

Plotting and comparing mutations (filtered using specific annotations), across bins in all chromosomes


```r
if (exists("binding2")){
  df <- binding2 
  df_new <- dplyr::mutate(df, chr_factor = factor(chromosome, levels = gtools::mixedsort(unique(chromosome))))
  head(df)
  bp <- ggplot(df_new, aes(x = bin_num, y = num_mut_in_bin)) +
    geom_bar(stat = "identity", fill="blue", colour = "pink") +
    facet_wrap(~chr_factor) +
    theme_bw(base_size = 10)
  bp + guides(fill=FALSE)
} else {
  print("No mutations present with the speicifed annotation")
}
```

![](tmb_coding_files/figure-html/mut_per_bin_per_chrom-1.png)<!-- -->

Plotting the total number of variants (filtered using specific annotations) across all chromosomes


```r
df <- data.frame(chr = chr_fac, mut_num = results3)
df
##    chr mut_num
## 1    1      14
## 2    2       7
## 3    3       6
## 4    4       2
## 5    5       4
## 6    6       5
## 7    7       8
## 8    8       4
## 9    9       4
## 10  10       5
## 11  11       5
## 12  12       6
## 13  13       4
## 14  14       3
## 15  15       1
## 16  16       7
## 17  17       3
## 18  18       1
## 19  19      10
## 20  20       1
## 21  21       1
## 22  22       3
## 23   X      11
## 24   Y       0
ggplot(df, aes(x = chr, y = mut_num)) +
  geom_point(colour = "purple") +
  theme_bw()
```

![](tmb_coding_files/figure-html/total_variants_per_chrom-1.png)<!-- -->

Plotting the mean of variants (filtered using specific annotations) across all bins in a chromosome i.e. results 


```r
df <- data.frame(chr = chr_fac, mut_num <- results)
ggplot(df, aes(x = chr, y = mut_num)) +
  geom_point(colour = "red") +
  theme_bw()
```

![](tmb_coding_files/figure-html/mean_of_varinats_per_bin_per_chrom-1.png)<!-- -->





