#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

hints:
 DockerRequirement:
   dockerImageId: umccr/pipeline-cwl

requirements:
  EnvVarRequirement:
    envDef:
      DEPLOY_ENV: $(inputs.denv)
  InitialWorkDirRequirement:
    listing:
      - $(inputs.samplesheet)

inputs:
    
  denv: string

  samplesheet:
    type: File
    inputBinding:
      position: 1

  config:
    type: Directory
    inputBinding:
      position: 2

outputs:

  split_samplesheets:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*[!.csv]"

baseCommand: [python, /scripts/samplesheet-check.py]
