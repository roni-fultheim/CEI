# CEI
A Cytoplasmic Index for Quantifying Immune-Related A-to-I RNA Editing

This repository contains both the **cloud-computing platform** and **analysis code** associated with our academic paper:

> **"TODO"**  
> Authors
> [DOI: 10.1234/example.doi.2025.001](https://doi.org/10.1234/example.doi.2025.001)

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
Make sure your machine has updated versions of [Nextflow](https://www.nextflow.io/docs/latest/install.html) and [Docker](https://docs.docker.com/engine/install/).    
In addition, the following should be available: `wget`, `curl`, `gzip`.     
Per-platform requirements:    
- AWS - download and initialize [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)      
- GCP - download and initialize [gcloud CLI](https://cloud.google.com/sdk/docs/install) (including `gsutil`)    

### Step 2 - Initialize Resources
Run the resource initialization script:
```bash
sh Init/init_main.sh "AWS"
```

### Step 3 - Run the Analysis
After initialization, launch the analysis workflow.

**1. Update user parameters configuration file:**
Change the user parameters within ``rna_editing.awsFargate.user_params.config``

**2. Run:**
```bash
nohup ~/nextflow -c rna_editing.awsFargate.config -bg run rna_editing.nf -profile <SE,stranded> --run_title <RUN_TITLE> --srrACC_list <SRR_LIST> > log.out 2> log.err &
```
For restricted access data (dbGaP, supported only on AWS)
```bash
nohup ~/nextflow -c rna_editing.awsFargate.config -bg run rna_editing.nf -profile <SE,stranded> --run_title <RUN_TITLE> --srrACC_list <SRR_LIST> --NGC_file <NGC_FILE> > log.out 2> log.err &
```

---

## Required Parameters Details
The following parameters must be provided, either directly via flags or within a configuration file:

For using AWS to run on SRA (see [Nextflow for Fargate documentation](https://www.nextflow.io/docs/latest/aws.html#aws-fargate)):
| Parameter  | Description |
|------------|-------------|
| `--ecr_region` | AWS region |
| `--process_queue` | AWS Batch queue for Fargate |
| `--jobRole` | AWS Batch job role |
| `--executionRole` | AWS Batch execution role |
| `--bucket_name` | Resources bucket |
| `--tower_access_token` | Nextflow [Seqera access token](https://www.nextflow.io/docs/latest/wave.html) |
| `--workspace_id` | Nextflow [Seqera workspace ID](https://www.nextflow.io/docs/latest/wave.html) |
| `--run_title` | Title of current run in general config, per-run |
| `--srrACC_list` | File with SRA accessions - in NF worklow, per-run |
| `--NGC_file` | NGC file for [restricted dbGaP access](https://www.ncbi.nlm.nih.gov/sra/docs/sra-dbGAP-cloud-download/) in general config, per-run, not required |

For using AWS to run on TCGA (see [Nextflow for Fargate documentation](https://www.nextflow.io/docs/latest/aws.html#aws-fargate)):
| Parameter  | Description |
|------------|-------------|
| `--ecr_region` | AWS region |
| `--resources_bucket_name` | Resources bucket |
| `--results_bucket_name` | Output bucket |
| `--process_spot_queue` | Description of parameter 2 |
| `--jobRole` | Description of parameter 2 |
| `--executionRole` | Description of parameter 2 |
| `--tower_access_token` | Description of parameter 2 |
| `--workspace_id` | Description of parameter 2 |
| `--GDC_token` | Description of parameter 2 |
| `--run_title` | Title of current run in general config, per-run |
| `--gdc_UUID_list` | Description of parameter 2 - in general config, per-run |

For GCP:
| Parameter  | Description |
|------------|-------------|
| `--project_name` | GCP project workspace name |
## TODO REGION
| `--bucket_name` | Resources bucket |
| `--run_title` | Title of current run in general config, per-run |
| `--srrACC_list` | File with SRA accessions - in NF worklow, per-run |



