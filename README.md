# Contents

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [Aim of the analysis](#aim-of-the-analysis)
* [Workflow](#workflow)

<!-- vim-markdown-toc -->

## Introduction

This report documents work done as part of exploring/executing [mosdepth](https://github.com/brentp/mosdepth), a tool that reports depth calculation for WGS, exome or target sequencing.

## Aim of the analysis

The basic idea is to flag if we are missing any key regions or abnormally highly covering others. 

The overall plan is to include the coverage assessment of hotspot regions (as used by umccrise `/data/cephfs/punim0010/extras/umccrise/genomes/GRCh37/hmf/Hotspot.tsv.gz` as well) in our bookdown report for patient samples.

## Workflow

1. Select a couple of patient samples - one with a poor cellualrity and the other with a decent one (Purple reports purity information
2. Run mosdepth on both samples
3. Evaluate results on targeted regions 
4. Report results in a document for further discussion
5. Work on automation of parameters
6. Include assessment in the bookdown report
