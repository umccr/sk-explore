#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

hints:
 DockerRequirement:
   dockerPull: umccr/bcl2fastq

requirements:
  EnvVarRequirement:
    envDef:
      DEPLOY_ENV: $(inputs.denv)
  InitialWorkDirRequirement:
    listing:
      - $(inputs.samplesheet_bcl2fastq)
      - entry: $(inputs.input_folder)
        writable: true

inputs:    
  denv: string

  input_folder:
    type: Directory
    inputBinding:
      position: 1
      prefix: -R

  output_folder:
    type: Directory
    inputBinding:
      position: 3
      prefix: -o
  
  samplesheet_bcl2fastq:
    type: [File, Directory]
    inputBinding:
      prefix: --sample-sheet

arguments: ["--no-lane-splitting"]

outputs:
  log_out:
    type: stdout

  samplesheets_fastq:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.fastq.gz"

stdout: bcl2fastq.log

baseCommand: []
