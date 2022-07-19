#!/usr/bin/env bash

# Set to fail
set -euo pipefail

# Globals
WORKFLOW_ID="wfl.3cd456b8d86a4f59b523627d3b6cb2a1"
WORKFLOW_VERSION="0.3.1"
PORTAL_TOKEN="eyJraWQiOiJiN002cGZKRzBIYWZhTmVvZW56ZTUwUCsydHBlK2c2UjF2eUpKQUhOdGo4PSIsImFsZyI6IlJTMjU2In0.eyJhdF9oYXNoIjoiLTZIZzVfamNUbnUweEVua2tSSjdaZyIsInN1YiI6ImYxMTZmMTcwLWExNWItNGE5Mi05M2IzLTAyZDM2MzEwMzRjYyIsImNvZ25pdG86Z3JvdXBzIjpbImFwLXNvdXRoZWFzdC0yX1N1T05KeWtJdF9Hb29nbGUiXSwiZW1haWxfdmVyaWZpZWQiOmZhbHNlLCJpc3MiOiJodHRwczpcL1wvY29nbml0by1pZHAuYXAtc291dGhlYXN0LTIuYW1hem9uYXdzLmNvbVwvYXAtc291dGhlYXN0LTJfU3VPTkp5a0l0IiwiY29nbml0bzp1c2VybmFtZSI6Ikdvb2dsZV8xMDA2OTM5MTM0NjI0NTk5MzA2MDciLCJhdWQiOiI3a2drNzcwNjgzZGxnczZsMmt1Z2Rkc2hnMCIsImlkZW50aXRpZXMiOlt7InVzZXJJZCI6IjEwMDY5MzkxMzQ2MjQ1OTkzMDYwNyIsInByb3ZpZGVyTmFtZSI6Ikdvb2dsZSIsInByb3ZpZGVyVHlwZSI6Ikdvb2dsZSIsImlzc3VlciI6bnVsbCwicHJpbWFyeSI6InRydWUiLCJkYXRlQ3JlYXRlZCI6IjE1OTE5MjUzNTg0MjAifV0sInRva2VuX3VzZSI6ImlkIiwiYXV0aF90aW1lIjoxNjU3NDk5OTE3LCJleHAiOjE2NTc2ODYzNTcsImlhdCI6MTY1NzU5OTk1NywiZW1haWwiOiJzZWhyaXNoLmthbndhbEB1bWNjci5vcmcifQ.WBPn1VHtBXKV6igKZ9gdL1UsQeUmMdY4Tv5fSGM8NP-u6Jh74uWaQ1dIlw2uVBxpM0JYrl062XUqCRk-g-cNwUwMFiNURMIe8RG8kSsgM4gMnqvUgVn7Y_ujlvlC5TC9fVtdv5CoNdJadtz2iAqUat6yTh6_uZlxY_JAuQ7pc5NvwexGxrhSxuoObDwraa0xv_IT8Qs5Cg6G-bDBBf7Mj6I9eAH7sXZf735tq3hMEM_22rDBzZrZlZo8epDy8OrEO4ivq3Rd0v7ih9mNyaKCpNE6SvMiTvLaug2iFmteev7SKVX5wefnS0nPe7FSnkrJh5sWjQKbi_PGAqOWlBHtBw"

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
