# CEI
A Cytoplasmic Index for Quantifying Immune-Related A-to-I RNA Editing

This repository contains both the **cloud-computing platform** and **analysis code** associated with our academic paper:

> **"TODO"**  
> Authors
> [DOI: 10.1234/example.doi.2025.001](https://doi.org/10.1234/example.doi.2025.001)


## 📑 Table of Contents
- [Overview](#-overview)
- [Repository Structure](#-repository-structure)
- [Getting Started](#-getting-started)
  - [Step 1 - Initialize the Machine](#step-1---initialize-the-machine)
  - [Step 2 - Initialize Resources](#step-2---initialize-resources)
  - [Step 3 - Run the Analysis](#step-3---run-the-analysis)
- [Profiles](#profiles)
- [Parameters Details](#parameters-details)
  - [AWS – SRA](#for-using-aws-to-run-on-sra-see-nextflow-for-fargate-documentation)
  - [AWS – TCGA](#for-using-aws-to-run-on-tcga-see-nextflow-for-fargate-documentation)
  - [GCP](#for-gcp)
- [Docker Images](#docker-images)
  - [Workflow](#Workflow)
  - [Initialization](#Initialization)

---

## 🧬 Overview

This project offers:
- A cloud-native computing platform built with Nextflow, Docker, and GCP/AWS for high-throughput RNA editing analysis
- Resources for quantifying both the **cytoplasmic editing index (CEI)** and **global Alu editing index (AEI)**
- Reproducible workflows for the analysis presented in the paper
- Example datasets and configurations for running your own analysis

---

## 📁 Repository Structure
CEI provides:

- A scalable cloud-native workflow for computing:
  - The **Cytoplasmic Editing Index (CEI)**  
  - The **Global Alu Editing Index ([AEI](https://www.nature.com/articles/s41592-019-0610-9))**
  - Editing within [known editing sites](https://doi.org/10.1038/s41467-022-28841-4)
  - Gene expression using [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html)
- Infrastructure support for:
  - **AWS and GCP** - note that the AWS pipeline is more advanced and we recommend using this version.
  - Parallelized processing of **SRA** (FASTQ) and **TCGA** (BAM) samples.
- Full analysis and plotting scripts for:
  - Benchmarking CEI vs global index
  - Immune-related editing response
  - Orientation-aware Alu element analysis

---

## 🚀 Getting Started

Before running the workflows, you must follow the initialization sequence below.

### Step 1 - Initialize the Machine
Make sure your machine has updated versions of [Nextflow](https://www.nextflow.io/docs/latest/install.html) and [Docker](https://docs.docker.com/engine/install/). Make sure your user account has permission to access the Docker daemon.    
In addition, the following should be available: `wget`, `curl`, `gzip`.     
Per-platform requirements:    
- AWS - download and initialize [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)      
- GCP - download and initialize [gcloud CLI](https://cloud.google.com/sdk/docs/install) (including `gsutil`)    

### Step 2 - Initialize Resources
Run the resource initialization script:
```bash
nohup sh CloudPipeline/Init/init_main.sh <PLATFORM> <REGION> <BUCKET_NAME> <NUM_THREADS>  > init.out 2> init.err &
```
| Parameter  | Position | Description | Allowed Values |
|------------|----------|-------------|----------------|
| `PLATFORM` | 1 | Which platform is used | `GCP` or `AWS` |
| `REGION` | 2 | Region in which the  resources bucket is to be located | Any region of platform (examples: us-central1, us-east1, etc.) |
| `BUCKET_NAME` | 3 | Wanted resources bucket name | Any name complying with platform criteria |
| `NUM_THREADS` | 4 | Number of threads for STAR and salmon index generation | Positive whole numbers |           

Default for number of threads for generation of STAR and Salmon indices is 10. STAR requires at least 64G RAM for this process.

### Step 3 - Run the Analysis
After initialization, launch the analysis workflow.     
Example for AWS - for GCP, use the files within the GCP directory.      
[Profiles](#profiles) and [parameters](#parameters-details) detailed below.        

**1. Update user parameters configuration file:**    
Change the user parameters within ``CloudPipeline/AWS/SRA_pipeline/rna_editing.user_params.config``

**2. Run:**
```bash
nohup ~/nextflow -c CloudPipeline/AWS/SRA_pipeline/rna_editing.config -bg CloudPipeline/AWS/SRA_pipeline/run rna_editing.nf -profile <SE,stranded,RL75,hg38> --run_title <RUN_TITLE> --srrACC_list <SRR_LIST> > log.out 2> log.err &
```
For restricted access data (dbGaP, supported only on AWS)
```bash
nohup ~/nextflow -c CloudPipeline/AWS/SRA_pipeline/rna_editing.config -bg run CloudPipeline/AWS/SRA_pipeline/rna_editing.nf -profile <SE,stranded,RL75,hg38> --run_title <RUN_TITLE> --srrACC_list <SRR_LIST> --NGC_file <NGC_FILE> > log.out 2> log.err &
```
---
    
## Profiles
The following profiles are supported. Any combination of profile options can be used from these categories, but at least one profile must be used from each category. See [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles) to learn more about profiles.       

| Category  | Options | Description | Steps Affected |
|-----------|---------|-------------|----------------|
| Genome | `hg38`, `mm10` | Genome of organism - human or mouse | Alignment, expression & editing quantification |
| Library type | `SE`, `PE` | Sequencing library type - single-end or paired-end | Preprocessing, alignment, expression & editing quantification |
| Sequencing directionality | `stranded`, `unstranded` | Is the sequencing directional? | Editing quantification |
| Sequencing length | `RL50`, `RL75`, `RL100`, `RL125`, `RL150` | Read length | Preprocessing & alignment |        

Note that stranded data can also be run as unstranded, without utilizing the strand information.     
Read length should match the general read length or below, as reads shorten than wanted length by 3bp or more will be filtered out (see protocol at (TODO add article reference)[].

    
---
    
## Parameters Details
The following parameters must be provided, either directly via flags or within a configuration file:

For using AWS to run on SRA (see [Nextflow for Fargate documentation](https://www.nextflow.io/docs/latest/aws.html#aws-fargate)):
| Parameter  | Description | Type   | Configuration File |
|------------|-------------|--------|--------------------|
| `--ecr_region` | AWS region | AWS parameter | User parameters config |
| `--process_queue` | AWS Batch queue for Fargate | AWS parameter | User parameters config |
| `--jobRole` | AWS Batch job role | AWS parameter | User parameters config |
| `--executionRole` | AWS Batch execution role | AWS parameter | User parameters config |
| `--tower_access_token` | Nextflow [Seqera access token](https://www.nextflow.io/docs/latest/wave.html) | Nextflow parameter | User parameters config |
| `--workspace_id` | Nextflow [Seqera workspace ID](https://www.nextflow.io/docs/latest/wave.html) | Nextflow parameter | User parameters config |
| `--bucket_name` | Resources bucket | AWS parameter | User parameters config |
| `--run_title` | Title of current run | Per-run parameter | General config |
| `--srrACC_list` | File with SRA accessions | Per-run parameter | General config |
| `--NGC_file` | NGC file for [restricted dbGaP access](https://www.ncbi.nlm.nih.gov/sra/docs/sra-dbGAP-cloud-download/), not required | Per-run parameter | General config |

For using AWS to run on TCGA (see [Nextflow for Fargate documentation](https://www.nextflow.io/docs/latest/aws.html#aws-fargate)):
| Parameter  | Description | Type   | Configuration File |
|------------|-------------|--------|--------------------|
| `--ecr_region` | AWS region | AWS parameter | User parameters config |
| `--process_spot_queue` | AWS Batch queue for Fargate | AWS parameter | User parameters config |
| `--jobRole` | AWS Batch job role | AWS parameter | User parameters config |
| `--executionRole` | AWS Batch execution role | AWS parameter | User parameters config |
| `--tower_access_token` | Nextflow [Seqera access token](https://www.nextflow.io/docs/latest/wave.html) | Nextflow parameter | User parameters config |
| `--workspace_id` | Nextflow [Seqera workspace ID](https://www.nextflow.io/docs/latest/wave.html) | Nextflow parameter | User parameters config |
| `--resources_bucket_name` | Resources bucket | AWS parameter | User parameters config |
| `--results_bucket_name` | Output bucket | AWS parameter | User parameters config |
| `--GDC_token` | GCD [authorization token](https://docs.gdc.cancer.gov/Data/Data_Security/Data_Security/) | User parameter | User parameters config |
| `--run_title` | Title of current run in general config | Per-run parameter | General config |
| `--gdc_UUID_list` | File with UUID accessions | Per-run parameter | General config |

For GCP:
| Parameter  | Description | Type   | Configuration File |
|------------|-------------|--------|--------------------|
| `--project_name` | GCP project workspace name | GCP parameter | User parameters config |
| `--region` | GCP region | GCP parameter | User parameters config |
| `--bucket_name` | Resources bucket | GCP parameter | User parameters config |
| `--run_title` | Title of current run in general config | Per-run parameter | General config |
| `--srrACC_list` | File with SRA accessions | Per-run parameter | General config |

    
---
    

## Docker Images
### Workflow
| Step  | Image | Source | 
|-------|-------|--------|
| Download SRA FASTQ files | levanonlab/sratoolkit:3.2.1 | [staphb/sratoolkit:3.2.1](https://hub.docker.com/r/staphb/sratoolkit/) |
| Preprocessing | levanonlab/fastp:0.23.4--hadf994f_2 | [quay.io/biocontainers/fastp:0.23.4--hadf994f_2](https://quay.io/repository/biocontainers/fastp?tab=tags&tag=0.23.4--hadf994f_2) |
| Expression quantification | levanonlab/salmon:1.10.2--hecfa306_0 | [quay.io/biocontainers/salmon:1.10.2--hecfa306_0](https://quay.io/repository/biocontainers/salmon?tab=tags&tag=1.10.2--hecfa306_0) |
| Alignment | levanonlab/star:2.7.10b--h9ee0642_0 | [quay.io/biocontainers/star:2.7.10b--h9ee0642_0](https://quay.io/repository/biocontainers/star?tab=tags&tag=2.7.10b--h9ee0642_0) |
| RNA editing index | levanonlab/rna-editing-index-lite:1.0 |   |
| RNA editing index & timing | levanonlab/rna-editing-index-lite:1.0.time |   |
| Per-site editing quantification | levanonlab/cmpileup:1.0 |   |
| Download TCGA BAM files | levanonlab/gdc-client:2.3.0 |  |

### Initialization       
| Step  | Image | Source | 
|-------|-------|--------|
| Resource processing | quay.io/biocontainers/bedtools | [quay.io/biocontainers/bedtools](https://quay.io/repository/biocontainers/bedtools?tab=tags&tag=latest) |
| Resource processing | quay.io/biocontainers/ucsc-bigbedtobed | [quay.io/biocontainers/ucsc-bigbedtobed](https://quay.io/repository/biocontainers/bedtools?tab=tags&tag=latest) |
| Region index generation | levanonlab/rna-editing-index-lite:1.0 | |
| Salmon index generation | levanonlab/salmon:1.10.2--hecfa306_0 | [quay.io/biocontainers/salmon:1.10.2--hecfa306_0](https://quay.io/repository/biocontainers/salmon?tab=tags&tag=1.10.2--hecfa306_0) |
| STAR index generation | levanonlab/star:2.7.10b--h9ee0642_0 | [quay.io/biocontainers/star:2.7.10b--h9ee0642_0](https://quay.io/repository/biocontainers/star?tab=tags&tag=2.7.10b--h9ee0642_0) |

