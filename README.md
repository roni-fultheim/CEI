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
- Resources for quantifying both the **cytoplasmic editing index** and **global Alu editing index**
- Reproducible workflows for the analysis presented in the paper
- Example datasets and configurations for running your own analysis

---

## 📁 Repository Structure
CEI provides:

- A scalable cloud-native workflow for computing:
  - The **Cytoplasmic Editing Index (CEI)**  
  - The **Global Alu Editing Index (AEI)**
  - Editing within [known editing sites](https://doi.org/10.1038/s41467-022-28841-4)
  - Gene expression
- Infrastructure support for:
  - **GCP and AWS**
  - Parallelized processing of **SRA** and **TCGA** samples
- Full analysis and plotting scripts for:
  - Benchmarking CEI vs global index
  - Immune-related editing response
  - Orientation-aware Alu element analysis

---

## 🚀 Getting Started

Before running the workflows, you must follow the initialization sequence below.

### Step 1 — Initialize the Machine
Run the machine initialization script:
TODO - check diff for AWS\GCP
```bash
sh CloudPipeline/Init/cloud_cmd.sh
```

### Step 2 — Initialize Resources
Run the resource initialization script:
TODO - check diff for AWS\GCP
```bash
sh Init/init_main.sh
```

### Step 3 — Run the Analysis
After initialization, launch the analysis workflow.

**Option 1 — Pass parameters via command line (example for AWS):**
```bash
cd <YOUR_WORKDIR>
nohup ~/nextflow -C rna_editing.awsFargate.config -bg run rna_editing.nf -profile SE,stranded --run_title <RUN_TITLE> --ecr_region <REGION> --ecr_user_id <ID> --bucket_namee <BUCKET> > log.out 2> log.err &

```

**Option 2 — Use a configuration file:**
Change the user parameters within ``rna_editing.awsFargate.user_params.config``, then run
```bash
nohup ~/nextflow -C rna_editing.awsFargate.config -bg run rna_editing.nf -profile SE,stranded --run_title <RUN_TITLE> > log.out 2> log.err &
```

---

### Required Parameters
The following parameters must be provided, either directly via `--param` flags or within a configuration file:

For using AWS to run on SRA:
| Parameter  | Description |
|------------|-------------|
| `--ecr_region` | Description of parameter 1 |
| `--ecr_user_id` | Description of parameter 2 |
| `--bucket_name` | Description of parameter 2 |
| `--process_queue` | Description of parameter 2 |
| `--jobRole` | Description of parameter 2 |
| `--executionRole` | Description of parameter 2 |
| `--tower_access_token` | Description of parameter 2 |
| `--workspace_id` | Description of parameter 2 |
| `--run_title` | Description of parameter 2 - in general config, per-run |
| `--NGC_file` | Description of parameter 2 - in general config, per-run, not required |

For using AWS to run on TCGA:
| Parameter  | Description |
|------------|-------------|
| `--ecr_region` | Description of parameter 1 |
| `--ecr_user_id` | Description of parameter 2 |
| `--resources_bucket_name` | Description of parameter 2 |
| `--results_bucket_name` | Description of parameter 2 |
| `--process_spot_queue` | Description of parameter 2 |
| `--jobRole` | Description of parameter 2 |
| `--executionRole` | Description of parameter 2 |
| `--tower_access_token` | Description of parameter 2 |
| `--workspace_id` | Description of parameter 2 |
| `--GDC_token` | Description of parameter 2 |
| `--run_title` | Description of parameter 2 - in general config, per-run |
| `--gdc_UUID_list` | Description of parameter 2 - in general config, per-run |


For GCP:
| Parameter  | Description |
|------------|-------------|
| `--ecr_region` | Description of parameter 1 |
| `--ecr_user_id` | Description of parameter 2 |
| `--bucket_name` | Description of parameter 2 |
| `--process_queue` | Description of parameter 2 |
| `--jobRole` | Description of parameter 2 |
| `--executionRole` | Description of parameter 2 |
| `--tower_access_token` | Description of parameter 2 |
| `--workspace_id` | Description of parameter 2 |

---

### 1. Cloud Setup
See instructions in `CloudPipeline/README.md` for:
- Setting up GCP/AWS environments
- Running the `SRA_pipeline` and `TCGA_pipeline`

