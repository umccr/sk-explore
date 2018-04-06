---
title: "Tumour Mutational Burden"
author: "Sehrish Kanwal"
date: "Thu 2018-Apr-05"
output: 
  html_document: 
    keep_md: yes
editor_options: 
  chunk_output_type: console
---



Required packages


```r
library(vcfR)
library(dplyr)
library(ggplot2)
```

Input file


```r
vcf <-  read.vcfR("../Data/ensemble-pon-pass.vcf") # 20,361 variants in total
## Scanning file to determine attributes.
## File attributes:
##   meta lines: 343
##   header line: 344
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

Declare bin size (in bases)


```r
bin_size <- 40000000
```

Get chromosome and filter vcf for that chromosome


```r
# takes a vcf and a chromosome
# outputs df with pos + chr for that chromosome
filter_chr <- function(x, vcf) {
  vcf_chr_pos <- data.frame(chr=getCHROM(vcf), pos=getPOS(vcf), stringsAsFactors = FALSE)
  return(filter(vcf_chr_pos, chr == x))
}
```

Prepare bins for the specific chromosome


```r
bin_chr <- function(chr_lengths, chr_name, bin_size) {
  return(seq(from = 0, to = chr_lengths[chr_name], by = bin_size))
}
```

Count mutations in a chromosome


```r
count_mut_in_chr <- function(chr_pos, bin_vec) {
  x <- cbind(chr_pos$pos, findInterval(chr_pos$pos, bin_vec))
  mean(table(x[, 2]))
}
```

Apply same function to all chromosomes. This gives a numerical vector containing the number of variants for each chromosome 


```r
results <- vector("numeric", length = length(chr_names))
for (i in 1:length(chr_names)) {
  results[i] <- count_mut_in_chr(filter_chr(chr_names[i], vcf), bin_chr(chr_lengths, chr_names[i], bin_size))
}
print(results)
##  [1] 163.14286 204.85714 231.00000 308.60000 300.40000 202.80000 282.50000
##  [8] 296.75000 176.50000 207.00000 218.75000 178.00000 348.66667 229.66667
## [15] 149.33333 187.33333  91.33333 292.50000 145.00000 142.00000 184.00000
## [22]  65.00000 613.75000   2.00000
```

Plotting the total number of variants across all chromosomes


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
df <- data.frame(chr = chr_fac, mut_num = c(unname(x)))
df
##    chr mut_num
## 1    1    1142
## 2    2    1434
## 3    3    1155
## 4    4    1543
## 5    5    1502
## 6    6    1014
## 7    7    1130
## 8    8    1187
## 9    9     706
## 10  10     828
## 11  11     875
## 12  12     712
## 13  13    1046
## 14  14     689
## 15  15     448
## 16  16     562
## 17  17     274
## 18  18     585
## 19  19     290
## 20  20     284
## 21  21     368
## 22  22     130
## 23   X    2455
## 24   Y       2
ggplot(df, aes(x = chr, y = mut_num)) +
  geom_point(colour = "purple") +
  theme_bw()
```

![](tmb_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

Plotting the mean of variants across all bins in a chromosome i.e. results 


```r
df <- data.frame(chr = chr_fac, mut_num <- results)
ggplot(df, aes(x = chr, y = mut_num)) +
  geom_point(colour = "red") +
  theme_bw()
```

![](tmb_files/figure-html/unnamed-chunk-10-1.png)<!-- -->




















