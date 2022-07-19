#!/usr/bin/env bash

# Set to fail
set -euo pipefail

# Globals
WORKFLOW_ID="wfl.3cd456b8d86a4f59b523627d3b6cb2a1"
WORKFLOW_VERSION="0.3.1"
PORTAL_TOKEN= "" # generate via https://data.umccr.org/

readarray -t cttso_bam_files <<< "$( \
  curl --silent --fail --location \
    --header "Authorization: Bearer $PORTAL_TOKEN" \
    --url "https://api.data.prod.umccr.org/workflows?end_status=Succeeded&type_name=tso_ctdna_tumor_only&rowsPerPage=1000" | \
  jq --raw-output \
    '
      .results |
      map(
        .output | fromjson |
        .output_results_dir_by_sample[0] |
        .listing |
        map(
          select(
            (.nameext == ".bam") and
            (.nameroot | startswith("evidence") | not ) and
            (.basename | endswith("cleaned.stitched.bam") | not )
          )
        ) |
        .[0].location
      ) |
      .[]
    ' \
)"

for bam_file in "${cttso_bam_files[@]}"; do
  # get basename without suffix
  bam_basename="$(basename "${bam_file}" .bam)"

  #set variables
  task_name="mosdepth-${bam_basename}"
  by_bed="gds://development/temp/sehrish/mosdepth/regions.bed"
  output_directory="gds://development/temp/sehrish/cttso/mosdepth_output/${bam_basename}"

  # collect bam path
  bam_path="$( \
    python -c \
      "from urllib.parse import urlparse; print(urlparse(\"${bam_file}\").path)" \
  )"

  # use the portal gds api to collect the internal file ids for the bam file and its respective index
  file_objs_json_str="$( \
    curl \
      --fail --silent --location \
      --header "Authorization: Bearer $PORTAL_TOKEN" \
      --url "https://api.data.prod.umccr.org/gds?search=${bam_path}" | \
    jq --raw-output '.results' \
  )"

  # from file_objs_json_str list, collect the bam file id and the index id (both are returns when searching the bam path
  bam_file_id="$( \
    jq --raw-output \
      '
        . |
        map(
             select(.path | endswith(".bam")) |
             .id
           ) |
        .[]
      ' <<< "${file_objs_json_str}" \
  )"

  bam_index_file_id="$( \
    jq --raw-output \
      '
        . |
        map(
             select(.path | endswith(".bam.bai")) |
             .id
           ) |
        .[]
      ' <<< "${file_objs_json_str}" \
  )"

  # use the portal's presigned url endpoint to collect the presigned urls
  bam_presigned_url="$( \
    curl \
      --fail --silent --location \
      --header "Authorization: Bearer $PORTAL_TOKEN" \
      --url "https://api.data.prod.umccr.org/gds/${bam_file_id}/presign" | \
    jq --raw-output \
      '.signed_url' \
    )"

  bam_index_presigned_url="$( \
    curl \
      --fail --silent --location \
      --header "Authorization: Bearer $PORTAL_TOKEN" \
      --url "https://api.data.prod.umccr.org/gds/${bam_index_file_id}/presign" | \
    jq --raw-output \
      '.signed_url' \
    )"

  # prepare input json
  input_json_str="$( \
    jq --null-input --raw-output --compact-output \
      --arg output_directory "${output_directory}" \
      --arg bam_presigned_url "${bam_presigned_url}" \
      --arg bam_index_presigned_url "${bam_index_presigned_url}" \
      --arg task_name "${task_name}" \
      --arg by_bed "${by_bed}" \
      '
        {
          "name": $task_name,
          "input": {
            "bam_sorted": {
              "class": "File",
              "location": $bam_presigned_url,
              "secondaryFiles": [
                {
                  "class": "File",
                  "location": $bam_index_presigned_url
                }
              ]
            },
            "by": {
              "class": "File",
              "location": $by_bed
            },
            "no_per_base": true
          },
          "engineParameters": {
            "outputDirectory": $output_directory
          }
        }
      ' \
  )"

  # launch task
  curl --silent --fail --location \
    --request POST \
    --url "https://aps2.platform.illumina.com/v1/workflows/${WORKFLOW_ID}/versions/${WORKFLOW_VERSION}:launch" \
    --header "Accept: application/json" \
    --header "Content-Type: application/json" \
    --header "Authorization: Bearer ${ICA_ACCESS_TOKEN}" \
    --data "${input_json_str}"
  # ica workflows versions launch
done
