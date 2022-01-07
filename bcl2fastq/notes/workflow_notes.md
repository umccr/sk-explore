<!-- vscode-markdown-toc -->
* 1. [Input data and scripts](#Inputdataandscripts)
* 2. [Workflow steps](#Workflowsteps)
	* 2.1. [1. Run folder check](#Runfoldercheck)
		* 2.1.1. [Run command:](#Runcommand:)
		* 2.1.2. [Output:](#Output:)
	* 2.2. [2. Sample sheet check](#Samplesheetcheck)
		* 2.2.1. [Run command:](#Runcommand:-1)
		* 2.2.2. [Output:](#Output:-1)
	* 2.3. [3. bcl2fastq](#bcl2fastq)
		* 2.3.1. [Run command:](#Runcommand:-1)
		* 2.3.2. [Output:](#Output:-1)
	* 2.4. [4. MultiQC](#MultiQC)
		* 2.4.1. [4.1. InterOp](#InterOp)
		* 2.4.2. [4.2. bcl2fastq](#bcl2fastq-1)
		* 2.4.3. [4.3. Output for both interop and bcl2fastq:](#Outputforbothinteropandbcl2fastq:)
	* 2.5. [5. Run folder checksums](#Runfolderchecksums)
		* 2.5.1. [Run command](#Runcommand)
		* 2.5.2. [Output:](#Output:-1)
	* 2.6. [6. Fastq checksums](#Fastqchecksums)
		* 2.6.1. [Run command:](#Runcommand:-1)
		* 2.6.2. [Output:](#Output:-1)
* 3. [Docker file](#Dockerfile)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

##  1. <a name='Inputdataandscripts'></a>Input data and scripts

The pipeline scripts are in `https://github.com/umccr/infrastructure/tree/master/scripts`

The workflow expects input directory with the following files:

```
├── CopyComplete.txt
├── RTAComplete.txt
├── SampleSheet.csv
└── SequenceComplete.txt
```

##  2. <a name='Workflowsteps'></a>Workflow steps

###  2.1. <a name='Runfoldercheck'></a>1. Run folder check

####  2.1.1. <a name='Runcommand:'></a>Run command:

```
$ export DEPLOY_ENV=dev
$ sh runfolder-check.sh ~/Documents/UMCCR/data/bcl2fastq/
```

####  2.1.2. <a name='Output:'></a>Output:

**Standard output:** written to a log file

`2019-09-20 14:29:32.467077000 runfolder-check.sh: INFO: Invocation with parameters: /Users/kanwals/Documents/UMCCR/data/bcl2fastq/`

###  2.2. <a name='Samplesheetcheck'></a>2. Sample sheet check

####  2.2.1. <a name='Runcommand:-1'></a>Run command:

```
$ python samplesheet-check.py ~/Documents/UMCCR/data/bcl2fastq/SampleSheet.csv
```

####  2.2.2. <a name='Output:-1'></a>Output:

**Standard output:** written to a log file

```
├── SampleSheet.csv.custom.1.truseq
├── SampleSheet.csv.custom.2.10X
```

###  2.3. <a name='bcl2fastq'></a>3. bcl2fastq

####  2.3.1. <a name='Runcommand:-1'></a>Run command:

```
$ sh run-bcl2fastq.sh -R ~/Documents/UMCCR/data/bcl2fastq/ -n bcl2fastq -o ~/Documents/UMCCR/data/bcl2fastq/output/
```

####  2.3.2. <a name='Output:-1'></a>Output:

**Standard output:** written to a log file

```
./output/
├── bcl2fastq_custom.1.truseq.log
└── bcl2fastq_custom.2.10X.log
```

###  2.4. <a name='MultiQC'></a>4. MultiQC

Uses container `docker pull umccr/multiqc`

Note: The `create-multiqc-reports.sh` script sets few default paths at the beginning of the script. 

My data setup and path looks like:

```
fastq_data_base_path="/Users/kanwals/Documents/UMCCR/data/bcl2fastq"
bcl_data_base_path="/Users/kanwals/Documents/UMCCR/data/bcl2fastq"
qc_output_base_path="/Users/kanwals/Documents/UMCCR/data/bcl2fastq/output"
```

####  2.4.1. <a name='InterOp'></a>4.1. InterOp

##### Run command:

```
sh create-multiqc-reports.sh bcl2fastq 190704_A00130_0109_BHLHLWDSXX
```

####  2.4.2. <a name='bcl2fastq-1'></a>4.2. bcl2fastq

##### Run command:

```
sh create-multiqc-reports.sh bcl2fastq 190704_A00130_0109_BHLHLWDSXX
```


####  2.4.3. <a name='Outputforbothinteropandbcl2fastq:'></a>4.3. Output for both interop and bcl2fastq:

**Standard output:** written to a log file

```
├── 190704_A00130_0109_BHLHLWDSXX_bcl2fastq_qc.html
├── 190704_A00130_0109_BHLHLWDSXX_bcl2fastq_qc_data
│   ├── multiqc.log
│   ├── multiqc_bcl2fastq_bylane.txt
│   ├── multiqc_bcl2fastq_bysample.txt
│   ├── multiqc_data.json
│   ├── multiqc_general_stats.txt
│   └── multiqc_sources.txt
├── 190704_A00130_0109_BHLHLWDSXX_interop_qc.html
└── 190704_A00130_0109_BHLHLWDSXX_interop_qc_data
    ├── multiqc.log
    ├── multiqc_data.json
    └── multiqc_sources.txt
```

###  2.5. <a name='Runfolderchecksums'></a>5. Run folder checksums

####  2.5.1. <a name='Runcommand'></a>Run command

```
sh create-checksums.sh runfolder ~/Documents/UMCCR/data/bcl2fastq/190704_A00130_0109_BHLHLWDSXX 
```

####  2.5.2. <a name='Output:-1'></a>Output:

**Standard output:** written to a log file

And a `runfolder.md5sum` file in `~/Documents/UMCCR/data/bcl2fastq/190704_A00130_0109_BHLHLWDSXX`

###  2.6. <a name='Fastqchecksums'></a>6. Fastq checksums

####  2.6.1. <a name='Runcommand:-1'></a>Run command:

```
sh create-checksums.sh bcl2fastq ~/Documents/UMCCR/data/bcl2fastq/output/
```

####  2.6.2. <a name='Output:-1'></a>Output:

**Standard output:** written to a log file

And a `bcl2fastq.md5sum` file in `~/Documents/UMCCR/data/bcl2fastq/output`


##  3. <a name='Dockerfile'></a>Docker file

```
$ /path_to/git/infrastructure/docker/umccr_pipeline $ docker build ../../scripts/umccr_pipeline/ --file ./Dockerfile
```
creates docker image (for my case image id is `8504550ff2fe`).


