# ECD feature pipeline for non-BS WGS data

## Overview

This repo contains functions for generating fragmentomics and image-based features from WGS data. 

The pipeline takes `.BAM` files as inputs. For each `.BAM`file, we first convert the file to a fragment-wise tab-separated file `.frag.tsv`. The `.frag.tsv` file contains information on each fragment, including chromosome, fragment start coordinate, fragment end coordinate, quality score, fragment length, read ID (or query name), ... Run `01_convert_BAM_to_FRAGS.sh`to convert `.BAM` to `.frag.tsv` file. 

We next construct 4bp end motifs and distance to nearest nucleosome for each fragment: run `02_calculate_features_from_frag.sh`.

To construct Copy number aberration (CNA) features on 1M bins and fragment length ratio in the whole genome, run the function `03_generate_binWise_features.R`. 

Finally, run `04_generate_features_fragmentomics_and_image_features.py` to summarize all features. 

See `run.debug.sh` for an example. 

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)


## Installation

Clone the repository and ensure all dependencies are installed:

```bash
git clone https://github.com/gsecddatalab/ecd_wgs_feature_pipeline.git
cd ecd_wgs_feature_pipeline
# Install required packages as needed
```

## Usage

Follow the steps below to run the pipeline:

1. Convert `.BAM` files to `.frag.tsv`:
    ```bash
    bash 01_convert_BAM_to_FRAGS.sh -i ${inputbam} -o ${outputdir} -n 20 -f ${path_to_fa} -q 50 -w 150 -e 151 -r 350 -c "true"
    ```
2. Calculate fragment features:
    ```bash
    bash 02_calculate_features_from_frag.sh -i ${outputdir}/${sampleid}/${sampleid}.frag.tsv -o ${outputdir} -f ${path_to_fa} -r ${path_to_nucleosome_ref} -c "true"
    ```
3. Generate bin-wise features. This function needs to be run in docker image: `tronghieunguyen/ecd_features`
    ```bash
    Rscript 03_generate_binWise_features.R \
    -s ${outputdir}/${sampleid}/short_long_BAM/${sampleid}.markdup.sorted_50_150.short.bam \
    -l ${outputdir}/${sampleid}/short_long_BAM/${sampleid}.markdup.sorted_50_150.long.bam \
    -f ${outputdir}/${sampleid}/${sampleid}.markdup.sorted.bam \
    -o ${outputdir}/${sampleid}/binwise_features;
    ```
4. Summarize all features:
    ```bash
    python 04_generate_features_fragmentomics_and_image_features.py \
    --input ${outputdir}/${sampleid}/${sampleid}.frag_output.tsv \
    --output ${outputdir}/${sampleid}/fragmentomics_features;
    ```

See `run.debug.sh` for a complete example.

## Features

- FLEN: Fragment length distribution (1D vector).
- EM: frequency of 4bp end motifs (1D vector).
- ND: distribution of nearest distance from fragment to nucleosome (1D vector).
- Image features: combination of FLEN; EM; ND features fragment-wise.

