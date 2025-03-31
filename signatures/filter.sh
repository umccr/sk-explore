#!/bin/bash
set -euo pipefail
module load bcftools
module load bedtools

DATA_DIR="data_duplex"
OUTPUT_DIR="output_filtered"

# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Loop over each sample directory
for sample_dir in "${DATA_DIR}"/*; do
    sample_name=$(basename "${sample_dir}")
    echo "Processing sample: ${sample_name}"

    # Assume the input VCF file is named: <sample_name>.1.consensus.variant-calls.genome.vcf.gz
    INPUT_VCF="${sample_dir}/${sample_name}.1.consensus.variant-calls.genome.vcf.gz"
    if [ ! -f "${INPUT_VCF}" ]; then
        echo "Input file ${INPUT_VCF} not found, skipping ${sample_name}"
        continue
    fi

    # Create a temporary filter directory within the sample directory
    FILTER_DIR="${sample_dir}/filter"
    mkdir -p "${FILTER_DIR}"

    #############################
    # 1. Initial filtering of low-frequency variants (VAF < 0.01)
    #############################

    # (1) Filter SNVs using HIAF as the reference
    SNV_VCF="${FILTER_DIR}/SNV_mutations.vcf.gz"
    bcftools view -i 'INFO/TYPE=="SNV" && INFO/HIAF < 0.01' "${INPUT_VCF}" -Oz -o "${SNV_VCF}"

    # (2) Filter Indels/Complex variants using ADJAF as the reference
    INDEL_VCF="${FILTER_DIR}/indel_mutations.vcf.gz"
    bcftools view -i '(INFO/TYPE=="Insertion" || INFO/TYPE=="Deletion" || INFO/TYPE=="Complex") && INFO/ADJAF < 0.01' "${INPUT_VCF}" -Oz -o "${INDEL_VCF}"

    #############################
    # 2. Exclude SNVs that fall within clonal MNV regions
    #############################

    # (1) Extract clonal MNVs (non-SNV variants with ADJAF or HIAF >= 0.01)
    CLONAL_MNV_VCF="${FILTER_DIR}/clonal_MNV.vcf.gz"
    bcftools view -i 'INFO/TYPE!="SNV" && (INFO/ADJAF>=0.01 || INFO/HIAF>=0.01)' "${INPUT_VCF}" -Oz -o "${CLONAL_MNV_VCF}"

    # (2) Convert clonal MNV VCF to BED format
    CLONAL_MNV_BED="${FILTER_DIR}/clonal_MNV.bed"
    bcftools query -f '%CHROM\t%POS0\t%END\n' "${CLONAL_MNV_VCF}" > "${CLONAL_MNV_BED}"

    # (3) Convert the initial filtered SNVs to BED format
    SNV_BED="${FILTER_DIR}/SNV_mutations.bed"
    bcftools query -f '%CHROM\t%POS0\t%POS\n' "${SNV_VCF}" > "${SNV_BED}"

    # (4) Use bedtools to exclude SNVs that overlap clonal MNV regions
    SNV_FILTERED_BED="${FILTER_DIR}/SNV_filtered.bed"
    bedtools intersect -v -a "${SNV_BED}" -b "${CLONAL_MNV_BED}" > "${SNV_FILTERED_BED}"

    # (5) Extract the final retained SNVs from the SNV VCF using the filtered BED file
    FINAL_SNV_VCF="${FILTER_DIR}/final_SNV_mutations.vcf.gz"
    bcftools view -T "${SNV_FILTERED_BED}" "${SNV_VCF}" -Oz -o "${FINAL_SNV_VCF}"
    
    # Create an index for FINAL_SNV_VCF (force creation)
    bcftools index -f "${FINAL_SNV_VCF}"

    #############################
    # 3. Final integration
    #############################

    # Create an index for INDEL_VCF if not already generated
    bcftools index -f "${INDEL_VCF}"

    FINAL_REPORT_VCF="${FILTER_DIR}/final_report.vcf.gz"
    # Concatenate FINAL_SNV_VCF and INDEL_VCF using the -a option to allow non-contiguous chromosome blocks
    bcftools concat -a "${FINAL_SNV_VCF}" "${INDEL_VCF}" -Oz -o "${FINAL_REPORT_VCF}"

    # Create a TBI index for the final report VCF (force TBI format using -t)
    bcftools index -f -t "${FINAL_REPORT_VCF}"

    # Create a sample output directory in OUTPUT_DIR and copy the final report and its index,
    # appending _filter to the file name.
    SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${sample_name}"
    mkdir -p "${SAMPLE_OUTPUT_DIR}"
    cp "${FINAL_REPORT_VCF}" "${SAMPLE_OUTPUT_DIR}/${sample_name}_filter.vcf.gz"
    cp "${FINAL_REPORT_VCF}.tbi" "${SAMPLE_OUTPUT_DIR}/${sample_name}_filter.vcf.gz.tbi"

    echo "Sample ${sample_name} processed. Final output: ${SAMPLE_OUTPUT_DIR}/${sample_name}_filter.vcf.gz and its index file."
done

echo "All samples processed."