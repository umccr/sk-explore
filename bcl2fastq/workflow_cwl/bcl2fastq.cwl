#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

hints:
 DockerRequirement:
   dockerImageId: umccr/pipeline
   dockerPull: umccr/bcl2fastq

requirements:
  ScatterFeatureRequirement: {}
  EnvVarRequirement:
    envDef:
      DEPLOY_ENV: $(inputs.denv)
#  InitialWorkDirRequirement:
#    listing:
#      - $(inputs.samplesheet)
#      - entry: $(inputs.input_folder)
#        writable: true

inputs:
  denv: string
  input_folder: Directory
  samplesheet: File
  config: Directory
  input_folder: Directory
  output_folder: Directory

steps:
  runFolderCheck:
    run: ./tools/runFolderCheck.cwl
    in:
      denv: denv
      input_folder: input_folder
    out:
      - log_out

  sampleSheetCheck:
    run: ./tools/sampleSheetCheck.cwl
    in:
      denv: denv
      samplesheet: samplesheet
      config: config
    out:
      [split_samplesheets]

  bcl2fastq:
    run: ./tools/bcl2fastq.cwl
    scatter: [samplesheet_bcl2fastq]
    in:
      denv: denv
      input_folder: input_folder
      output_folder: output_folder
      samplesheet_bcl2fastq: sampleSheetCheck/split_samplesheets
    out:
      [samplesheets_fastq]

outputs:
  pipeline_split_samplesheets:
    type:
      type: array
      items: File
    outputSource: sampleSheetCheck/split_samplesheets

  pipeline_result:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: bcl2fastq/samplesheets_fastq

