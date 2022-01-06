---
title: "Coverage assessment using mosdepth"
author: "Sehrish Kanwal"
date: "Tue 2019-Feb-12"
output:
  html_document: 
    keep_md: true
editor_options:
  chunk_output_type: console
---





## Required packages


```r
library(knitr)
library(dplyr)
library(kableExtra)
```


```r
lf <- function(...) {
  data.frame(fname = list.files(...)) %>% 
    knitr::kable(row.names = TRUE)
}
```

Currently focussing on `regions.bed` file that contains information about coverage of hotspot regions and `quantized.bed` file that contains information about callability windows.

### Hotspot regions


```r
#prepare input data location
mosdepth <- "/Users/kanwals/Documents/UMCCR/data/mosdepth"

#read in sample directory containing mosdepth outputs
mosdepth_output <- file.path(mosdepth, "CCR180136_WH18F001P025_fakeY_changedthreshold/")
#call lf function to list files in the directory
lf(mosdepth_output)
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> fname </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> CCR180136_WH18F001P025.mosdepth.global.dist.txt </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> CCR180136_WH18F001P025.mosdepth.region.dist.txt </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> CCR180136_WH18F001P025.quantized.bed.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:left;"> CCR180136_WH18F001P025.quantized.bed.gz.csi </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> CCR180136_WH18F001P025.regions.bed.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:left;"> CCR180136_WH18F001P025.regions.bed.gz.csi </td>
  </tr>
</tbody>
</table>

```r

#read regions.bed file - first grep the appropriate bed file from mosdepth sample output folder and then create a path for that name
regions_bed_name <- list.files(mosdepth_output)[grep(".regions.bed.gz$",list.files(mosdepth_output))]

#prepare file path for the regions.bed file
regions_bed_path <- file.path(mosdepth_output, regions_bed_name)
regions_bed <- as.data.frame(read.table(regions_bed_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
head(regions_bed)
##   V1      V2      V3  V4 V5
## 1  1 8073509 8073510   C 82
## 2  1 9777666 9777667   C 86
## 3  1 9777666 9777667   C 86
## 4  1 9780851 9780852   G 70
## 5  1 9780851 9780854 GAG 70
## 6  1 9787030 9787031   G 71

#read hotspots file - to include the column that shows the change from ref. nucleotide in the regions.bed file
hotspots <- as.data.frame(read.table("~/Documents/UMCCR/data/mosdepth/hotspots/hotspot_fake.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE))

#merge region_bed file with hotspots column that has info about the possible nucleotide change.Also avoid adding an additional Row.names column by setting it to NULL, The output of this has chr, start, end, ref, cov and alt. Rearrange this output and exclude end (V3) column 
regions_bed_hotspot <- transform(merge(regions_bed, hotspots$V5, by = "row.names", sort = FALSE), Row.names=NULL) %>%
  dplyr::select(V1, V2, V4, y, V5)
#assign column names  
colnames(regions_bed_hotspot) <- c("chrom","pos", "ref", "alt", "cov")
head(regions_bed_hotspot)
##   chrom     pos ref alt cov
## 1     1 8073509   C   A  82
## 2     1 9777666   C   A  86
## 3     1 9777666   C   G  86
## 4     1 9780851   G   A  70
## 5     1 9780851 GAG AAA  70
## 6     1 9787030   G   A  71

#add callability column to regions.bed file - based on the coverage value
call = vector()
call <- sapply(regions_bed_hotspot$cov, function(x){
  if (x == 0)
  {call = c(call, "0")} 
  else if (x >= 1 && x <= 15) 
  {call = c(call, "1:15")}
  else if (x >= 15 && x <= 500) 
  {call = c(call, "15:500")}
  else if (x >= 500) 
  {call = c(call, "500:inf")}
})
regions_bed_new <- dplyr::mutate(regions_bed_hotspot, callability = call)

#present output in a table format - only showing first few entries in the output because the .md file size exceeds github file limit 
head(regions_bed_new) %>%
  kable(caption = "Coverage value and the corresponding callability level for hotspot regions") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  scroll_box(height = "400px")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-y: scroll; height:400px; "><table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Coverage value and the corresponding callability level for hotspot regions</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> chrom </th>
   <th style="text-align:right;"> pos </th>
   <th style="text-align:left;"> ref </th>
   <th style="text-align:left;"> alt </th>
   <th style="text-align:right;"> cov </th>
   <th style="text-align:left;"> callability </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 8073509 </td>
   <td style="text-align:left;"> C </td>
   <td style="text-align:left;"> A </td>
   <td style="text-align:right;"> 82 </td>
   <td style="text-align:left;"> 15:500 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 9777666 </td>
   <td style="text-align:left;"> C </td>
   <td style="text-align:left;"> A </td>
   <td style="text-align:right;"> 86 </td>
   <td style="text-align:left;"> 15:500 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 9777666 </td>
   <td style="text-align:left;"> C </td>
   <td style="text-align:left;"> G </td>
   <td style="text-align:right;"> 86 </td>
   <td style="text-align:left;"> 15:500 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 9780851 </td>
   <td style="text-align:left;"> G </td>
   <td style="text-align:left;"> A </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:left;"> 15:500 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 9780851 </td>
   <td style="text-align:left;"> GAG </td>
   <td style="text-align:left;"> AAA </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:left;"> 15:500 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 9787030 </td>
   <td style="text-align:left;"> G </td>
   <td style="text-align:left;"> A </td>
   <td style="text-align:right;"> 71 </td>
   <td style="text-align:left;"> 15:500 </td>
  </tr>
</tbody>
</table></div>

```r

#only output loci with zero/low coverage
low_coverage_hotspot <- data.frame()
for (row in 1:nrow(regions_bed_new)){
  if(regions_bed_new[row,6] %in% c("0","1:15")){
    #low_coverage_hotspot <- rbind[low_coverage_hotspot,regions_bed_new[row,]] This throws a stupid error
    low_coverage_hotspot <- dplyr::bind_rows(low_coverage_hotspot,regions_bed_new[row,])
  }
}

#Number of hotspots with zero/low coverage
nrow(low_coverage_hotspot)
## [1] 1

#present output in a table format
head(low_coverage_hotspot, 10) %>%
  kable(caption = "Coverage for zero and low callability hotspot regions - only displaying the first few values") %>%
  kable_styling(bootstrap_options = "striped", full_width = F) %>%
  scroll_box(height = "400px")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-y: scroll; height:400px; "><table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Coverage for zero and low callability hotspot regions - only displaying the first few values</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> chrom </th>
   <th style="text-align:right;"> pos </th>
   <th style="text-align:left;"> ref </th>
   <th style="text-align:left;"> alt </th>
   <th style="text-align:right;"> cov </th>
   <th style="text-align:left;"> callability </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Y </td>
   <td style="text-align:right;"> 1271365 </td>
   <td style="text-align:left;"> ATG </td>
   <td style="text-align:left;"> GTA </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
  </tr>
</tbody>
</table></div>

## Callability levels

Here, analysing the quantized output produced by mosdepth to see if it adds any extra information.


```r
#read quatized.bed file - first grep the appropriate bed file from folder and then create a path for that name
quantized_bed_name <- list.files(mosdepth_output)[grep(".quantized.bed.gz$",list.files(mosdepth_output))]
quantized_bed_path <- file.path(mosdepth_output, quantized_bed_name)
quantized_bed <- as.data.frame(read.table(quantized_bed_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
#assign column names
colnames(quantized_bed) <- c("Chr", "Start", "End", "Callability")
#number of values
nrow(quantized_bed)
## [1] 345456

#Extract no coverage regions from quantized.bed file
no_coverage <- quantized_bed[quantized_bed$Callability == "0", ]
nrow(no_coverage)
## [1] 0
if(nrow(no_coverage) != 0)
  head(no_coverage) %>%
    kable(caption = "No coverage regions - only showing first few values") %>%
    kable_styling(bootstrap_options = "striped", full_width = F)
  
#Extract low coverage regions from quantized.bed file
low_coverage <- quantized_bed[quantized_bed$Callability == "1:15", ]
nrow(low_coverage)
## [1] 144467
if(nrow(low_coverage) != 0)
  head(low_coverage) %>%
    kable(caption = "Low coverage regions- only showing first few values") %>%
    kable_styling(bootstrap_options = "striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Low coverage regions- only showing first few values</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Chr </th>
   <th style="text-align:right;"> Start </th>
   <th style="text-align:right;"> End </th>
   <th style="text-align:left;"> Callability </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 9994 </td>
   <td style="text-align:right;"> 9999 </td>
   <td style="text-align:left;"> 1:15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 14 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 10621 </td>
   <td style="text-align:right;"> 10637 </td>
   <td style="text-align:left;"> 1:15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 10672 </td>
   <td style="text-align:right;"> 10789 </td>
   <td style="text-align:left;"> 1:15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 18 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 10820 </td>
   <td style="text-align:right;"> 10856 </td>
   <td style="text-align:left;"> 1:15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 28 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 18847 </td>
   <td style="text-align:right;"> 18853 </td>
   <td style="text-align:left;"> 1:15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 36 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 26993 </td>
   <td style="text-align:right;"> 27054 </td>
   <td style="text-align:left;"> 1:15 </td>
  </tr>
</tbody>
</table>

```r

#Extract high coverage regions from quantized.bed file
high_coverage <- quantized_bed[quantized_bed$Callability == "500:inf", ]
nrow(high_coverage)
## [1] 25318
if(nrow(high_coverage) != 0)
  head(high_coverage) %>%
    kable(caption = "High coverage regions- only showing first few values") %>%
    kable_styling(bootstrap_options = "striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>High coverage regions- only showing first few values</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Chr </th>
   <th style="text-align:right;"> Start </th>
   <th style="text-align:right;"> End </th>
   <th style="text-align:left;"> Callability </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 10010 </td>
   <td style="text-align:right;"> 10230 </td>
   <td style="text-align:left;"> 500:inf </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 10231 </td>
   <td style="text-align:right;"> 10238 </td>
   <td style="text-align:left;"> 500:inf </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 10241 </td>
   <td style="text-align:right;"> 10244 </td>
   <td style="text-align:left;"> 500:inf </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 10257 </td>
   <td style="text-align:right;"> 10311 </td>
   <td style="text-align:left;"> 500:inf </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 10353 </td>
   <td style="text-align:right;"> 10460 </td>
   <td style="text-align:left;"> 500:inf </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 13506 </td>
   <td style="text-align:right;"> 13516 </td>
   <td style="text-align:left;"> 500:inf </td>
  </tr>
</tbody>
</table>



