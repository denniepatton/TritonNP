# *Snakemake workflow for TritonNP*

## Description
From *Triton*, the Greek god of waves, and N(ucleosome) P(hasing). This workflow will run TritonNP on matched BAM and (optional) GC-correction files for regions of interest, generating nucleosome phasing metrics and region-level figures as outputs.

This tool has been updated! For the newest version, which now includes fragmentomic analyses, check out: https://github.com/GavinHaLab/Triton  
However, f you are looking to reproduce LuCaP PDX phasing analyses from https://doi.org/10.1158/2159-8290.CD-22-0692 you are in the right place.

## Contact
Robert Patton

Fred Hutchinson Cancer Center

Contact: <rpatton@fredhutch.org>

Date: June 21, 2022

Website: https://github.com/denniepatton

## Requirements
### Software packages or libraries
  - Python 3.7.4
    - Snakemake 5.19.2
    - Pysam 0.15.4
    - Argparse 1.1
    - Numpy 1.17.3
    - Pandas 0.25.3
    - Scipy 1.4.1

### Scripts
  - GenerateFFTFeatures.py (primary tool for extracting phasing features)
  - CombinePhasingFM.py (combines output files into a single feature matrix)
  - GenerateGCSamples.py (optional; found in config/; used for combining sample info with GC-bias data from Griffin)

### Griffin-based GC correction
TritonNP optionally takes BAM-matched GC bias data produced by the Griffin workflow; the workflow with instructions for generating bias files can be
found at https://github.com/GavinHaLab/Griffin (when used in the snakemake as opposed to a stand-alone tool GC bias is required).

## Sample list
Sample names with paths to matching BAM and GC_bias files should be defined in a YAML file. See `PhasingAnalysis/config/samples.yaml` for an example.
This file may also be generated automatically based on a directory full of BAMs and a path to matched GC_bias files using `GenerateGCSamples.py` as a
stand-alone script.
```
samples:
  sample_1:
    bam: /path/to/sample_1.bam
    GC_bias: /path/to/sample_1.GC_bias.txt
```

## Snakefiles
1. `GenerateFFTFeatures.snakefile`

### [config.yaml] (PhasingAnalysis/config/config.yaml)
See below for details on [config.yaml] (`PhasingAnalysis/config/config.yaml`)

## Run the analysis
### 1. Invoking the full snakemake workflow on a cluster using 'slurm'
There is only one file in use for `slurm`:
  `PhasingAnalysis/config/cluster_slurm.yaml` - This file contains memory, runtime, and number of cores for each task.
To invoke the snakemake pipeline for `qsub`:
```
snakemake -s GenerateFFTFeatures.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40
```

### 2. Using the tool as a stand-alone
GenerateFFTFeatures.py may also be run as a stand-alone script to retrieve features of an individual BAM file, with or without GC_bias, on a local
machine. Once dependencies are loaded run `python GenerateFFTFeatures.py -h` (help) to see all available options, which have more detailed explanations
below in the configuration and settings section.

## Configuration and settings
All (default) settings for the workflow are contained in [config.yaml] (`PhasingAnalysis/config/config.yaml`), while cluster settings and sample
paths are found in `PhasingAnalysis/config/cluster_slurm.yaml` and `PhasingAnalysis/config/samples.yaml`, respectively.

### sites_bed
Path to BED file containing regions of interest

### results_dir
Path to the output directory (default: results)

### reference_genome
Path to GRCh38 reference genome for GC calculation (.fa)

### plot_list
Text file with a list of region-names to generate plots for (may be left empty)

### size_range
Tuple of ints for minimum and maximum fragment sizes to use; defaults to (15, 500)

### map_quality
Int defining the minimum read mapping qualilty to keep

## Software License
TritonNP Copyright (c) 2022 Fred Hutchinson Cancer Research Center
All rights reserved.

This program is free software: you can redistribute it and/or modify it under the terms of the BSD-3-Clause-Clear license. No licenses are granted to any
patent rights of the Fred Hutchinson Cancer Research Center.  

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the BSD-3-Clause-Clear license for more details.  

You should have received a copy of the G BSD-3-Clause-Clear license along with this program.
If not, see https://spdx.org/licenses/BSD-3-Clause-Clear.html. 
