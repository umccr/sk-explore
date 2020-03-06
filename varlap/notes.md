## vcf manipulation

Documenting useful vcf manipulation commands I have been playing with recently for varlap. Aiming to add class information in individual vcfs and then merging these together as one final vcf with required `Class` information.

#### Update vcf to keep a specific sample

```
$ vcftools --remove-indv CCR180135_WH18B002P025 --vcf 0002.vcf --recode --out ./both_P025.vcf
```

#### Update samplename in a vcf file

```
$ bcftools reheader -s sample.txt -o ./dragen_renamedSample.vcf dragen.vcf
```

where `sample.txt` specifies new sample names, one name per line, in the same order as they appear in the VCF file. 


#### Compress vcf

```
$ bcftools view both_P025.vcf -Oz -o both_P025.vcf.gz
```

#### Concatenate vcfs

```
$ bcftools concat -a ./bcbio_P025.vcf.gz ./dragen_P025.vcf.gz ./both_P025.vcf.gz -o ./result.vcf
```

