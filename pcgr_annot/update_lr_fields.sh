#!/bin/bash

# Check if an input file is provided
if [ -z "$1" ]; then
echo "Usage: $0 <input_file.vcf.gz>"
exit 1
fi

# Get the input file name without the extension
input_file_name="${1%.*}"
filename=$(basename -- "$input_file_name")
input_file_name="${filename%.*}"
echo "$input_file_name"

# Extract FORMAT/VAF and FORMAT/DP
bcftools query -f '%CHROM\t%POS\t[ %VAF]\n' "$1" | bgzip -c > "${input_file_name}_annot-vaf.txt.gz"
bcftools query -f '%CHROM\t%POS\t[ %DP]\n' "$1" | bgzip -c > "${input_file_name}_annot-dp.txt.gz"

# Index the file
tabix -s1 -b2 -e2 "${input_file_name}_annot-vaf.txt.gz"
tabix -s1 -b2 -e2 "${input_file_name}_annot-dp.txt.gz"

# Create header line
echo -e '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Fractions">' >> hdr-vaf.txt
echo -e '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">' >> hdr-dp.txt
echo -e '##FORMAT=<ID=VDP,Number=1,Type=Integer,Description="Read Depth">' >> hdr-dp.txt

# Transfer the annotation and header
bcftools annotate -a "${input_file_name}_annot-vaf.txt.gz" -h hdr-vaf.txt -c CHROM,POS,INFO/VAF "$1" > "${input_file_name}_annot.vcf"
bcftools annotate -a "${input_file_name}_annot-dp.txt.gz" -h hdr-dp.txt -c CHROM,POS,INFO/DP "${input_file_name}_annot.vcf" > "${input_file_name}_annot_2.vcf"

# Update INFO DP and VAF field
bcftools annotate -c INFO/VDP:=INFO/DP "${input_file_name}_annot_2.vcf" > "${input_file_name}_annot_VDP.vcf"

# Zip and index the final file
bgzip -c "${input_file_name}_annot_VDP.vcf" > "${input_file_name}_annot_VDP.vcf.gz"
tabix -p vcf "${input_file_name}_annot_VDP.vcf.gz"

echo "Output file: ${input_file_name}_annot_VDP.vcf.gz"
