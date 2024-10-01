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

# Extract FORMAT/VF
bcftools query -f '%CHROM\t%POS\t[ %VF]\n' "$1" | bgzip -c > "${input_file_name}_annot-test.txt.gz"

# Index the file
tabix -s1 -b2 -e2 "${input_file_name}_annot-test.txt.gz"

# Create header line
echo -e '##INFO=<ID=VF,Number=1,Type=Float,Description="Variant Frequency">' >> hdr.txt

# Transfer the annotation and header
bcftools annotate -a "${input_file_name}_annot-test.txt.gz" -h hdr.txt -c CHROM,POS,INFO/VF "$1" > "${input_file_name}_annot.vcf"

# Update INFO DP field
bcftools annotate -c INFO/VDP:=INFO/DP "${input_file_name}_annot.vcf" > "${input_file_name}_annot_VDP.vcf"

# Extract pass only variants
bcftools view -f PASS "${input_file_name}_annot_VDP.vcf" > "${input_file_name}_annot_VDP_final.vcf"

# Zip and index the final file
bgzip -c "${input_file_name}_annot_VDP_final.vcf" > "${input_file_name}_annot_VDP_final.vcf.gz"
tabix -p vcf "${input_file_name}_annot_VDP_final.vcf.gz"

echo "Output file: ${input_file_name}_annot_VDP_final.vcf.gz"
