denv: 'prod'

input_folder:
  class: Directory
  location: /data/bcl/

samplesheet:
  class: File
  path: /data/bcl/SampleSheet.csv

config:
  class: Directory
  location: /home/ssm-user/.config/gspread_pandas

output_folder: 
  class: Directory
  location: /data/bcl/output
