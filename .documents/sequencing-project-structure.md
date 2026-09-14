# Project structure — sequencing pipelines

Single Python + R project to manage multiple sequencing types (RNA-seq, WGS/WES, metagenomics, ChIP-seq), with a shared core and one folder per pipeline.

## Repository tree

```
project/
├── .gitignore
├── README.md
│
├── shared/
│   ├── envs/
│   │   ├── base.yml              # common tools: fastqc, fastp, samtools, multiqc
│   │   └── base.renv.lock
│   ├── lib/
│   │   ├── python/
│   │   │   ├── __init__.py
│   │   │   ├── fastq_io.py        # read/write/validate FASTQ, read counting
│   │   │   ├── alignment_utils.py # wrappers around samtools, BAM stats parsing
│   │   │   ├── reference_utils.py # download/manage genomes and annotations
│   │   │   ├── qc_utils.py        # FastQC/MultiQC report parsing, QC thresholds
│   │   │   └── plotting.py        # plots reused across pipelines
│   │   └── r/
│   │       ├── data_loading.R     # load counts/variants/tables
│   │       ├── stats_helpers.R    # common statistical test wrappers
│   │       └── plotting_theme.R   # shared ggplot theme, color palette
│   └── reference/                 # shared genomes/annotations (NOT in git)
│
├── pipelines/
│   ├── rna-seq/
│   │   ├── config.yml
│   │   ├── Snakefile
│   │   ├── envs/
│   │   │   └── rna-seq.yml        # extends shared/envs/base.yml
│   │   └── scripts/
│   │       ├── python/
│   │       └── r/
│   │
│   ├── wgs/
│   │   ├── config.yml
│   │   ├── Snakefile
│   │   ├── envs/
│   │   │   └── wgs.yml
│   │   └── scripts/
│   │       ├── python/
│   │       └── r/
│   │
│   ├── metagenomics/
│   │   ├── config.yml
│   │   ├── Snakefile
│   │   ├── envs/
│   │   │   └── metagenomics.yml
│   │   └── scripts/
│   │       ├── python/
│   │       └── r/
│   │
│   └── chip-seq/
│       ├── config.yml
│       ├── Snakefile
│       ├── envs/
│       │   └── chip-seq.yml
│       └── scripts/
│           ├── python/
│           └── r/
│
├── data/                          # NOT in git (.gitignore)
│   ├── raw/
│   │   ├── rna-seq/
│   │   ├── wgs/
│   │   ├── metagenomics/
│   │   └── chip-seq/
│   └── reference/                  # if not shared in shared/
│
├── notebooks/                      # free exploration, not part of the pipelines
│
└── results/
    ├── rna-seq/{data-run}/
    ├── wgs/{data-run}/
    ├── metagenomics/{data-run}/
    └── chip-seq/{data-run}/
```

## Structure principles

- **`shared/`** holds everything common to multiple pipelines: base environments, reusable functions, shared reference genomes. If the same function ends up in two pipelines, it belongs here.
- **`pipelines/{name}/envs/`** stays "thin": it extends `shared/envs/base.yml` and adds only pipeline-specific tools (e.g. STAR for RNA-seq, GATK for WGS), instead of reinstalling everything from scratch.
- **`data/`** and reference genomes are **not versioned in git** (too large): add them to `.gitignore`, with a download script to repopulate them.
- **`results/{pipeline}/{data-run}/`** separates outputs by pipeline and by individual run, avoiding overwrites.
- Adding a new pipeline (e.g. long-read/Nanopore) only means creating `pipelines/long-read/`, reusing the existing `shared/` code.

## General flow of each pipeline

```
Raw FASTQ
   → QC and trimming (FastQC / fastp / MultiQC)
   → Alignment (STAR / BWA / Bowtie2 / salmon)
   → Analysis-specific output (gene counts / variants / peaks / taxonomy)
   → Statistical analysis in R (DESeq2 / VariantAnnotation / phyloseq / ChIPseeker)
```

## Tool summary by category

### Core tools (common to all pipelines)

| Tool | Input | Output | Purpose |
|---|---|---|---|
| FastQC | Raw FASTQ | HTML report | QC of raw reads |
| MultiQC | Reports from other tools | Aggregated HTML report | Summarizes QC across all samples |
| fastp | Raw FASTQ | Cleaned FASTQ + report | Removes adapters/low-quality bases |
| samtools | SAM/BAM | Sorted/indexed BAM | Alignment file manipulation |
| bcftools | VCF/BCF | Filtered VCF | Variant file manipulation |
| seqtk | FASTQ | Downsampled FASTQ | Subsampling for quick testing |

### RNA-seq

| Tool | Input | Output | Purpose |
|---|---|---|---|
| STAR | FASTQ + genome | BAM | Accurate alignment (RAM-heavy) |
| salmon / kallisto | FASTQ + transcriptome | Abundance table | Lightweight pseudo-alignment |
| featureCounts | BAM + GTF | Count matrix | Counts reads per gene |
| DESeq2 / edgeR (R) | Count matrix | Differentially expressed genes | Statistical analysis |

### WGS/WES

| Tool | Input | Output | Purpose |
|---|---|---|---|
| BWA-MEM | FASTQ + genome | BAM | DNA alignment |
| GATK | BAM | VCF | Variant calling (heavy) |
| VEP / ANNOVAR | VCF | Annotated VCF/table | Variant annotation |
| VariantAnnotation (R) | VCF | Annotated R objects | Variant analysis in R |

### Metagenomics / 16S

| Tool | Input | Output | Purpose |
|---|---|---|---|
| QIIME2 | FASTQ (16S) | Abundance tables, taxonomy | Full microbiome pipeline |
| DADA2 (R) | FASTQ (16S) | ASVs | High-resolution bacterial sequence identification |
| Kraken2 | FASTQ (shotgun) | Per-read taxonomic classification | Identifies the origin of each read |
| phyloseq (R) | Abundance tables | Plots, statistics | Microbiome analysis |

### ChIP-seq / epigenetics

| Tool | Input | Output | Purpose |
|---|---|---|---|
| Bowtie2 | FASTQ + genome | BAM | Short-read alignment |
| MACS2 | BAM | Peaks | Finds protein-DNA binding regions |
| ChIPseeker (R) | Peaks | Annotation, plots | Biological interpretation of peaks |

### Environment and workflow management

| Tool | Input | Output | Purpose |
|---|---|---|---|
| conda / mamba | `.yml` file | Isolated environment | Python dependency management |
| renv (R) | R scripts | `renv.lock` | Reproducible R dependency management |
| Snakemake / Nextflow | Rules/pipeline | Automatic step execution | Pipeline orchestration |
| DVC / git-lfs | Large files | Data versioning | Manages files too large for plain git |

## Practical notes (local environment, laptop)

- Prefer lightweight tools where possible (salmon/kallisto instead of STAR, bcftools instead of GATK when sufficient)
- Develop/debug on subsamples with `seqtk` before running pipelines on full datasets
- If datasets grow, Snakemake/Nextflow allow moving execution to a cluster/cloud without rewriting the pipelines

## When to split code into multiple files (Python/R)

Rule of thumb: **split by functional responsibility, not by file size.**

**Split into separate files when:**
- The function is reused across multiple pipelines (belongs in `shared/lib/`)
- It represents a logically distinct pipeline step (e.g. `load_data.R`, `normalize.R`, `differential_expression.R`)
- It has heavy dependencies that other code shouldn't have to import
- It's meant to grow into an internal package/library
- You plan to write targeted tests mirroring the source files

**Keep in one file when:**
- It's a throwaway/exploratory script
- A handful of small, tightly related helper functions with no reuse elsewhere
- It's the pipeline's entry point (orchestration "glue" code), even if it calls functions defined elsewhere

**Heuristic**: ask "if I change this function, do I usually change the others in the same file too?" If yes, keep them together; if they evolve independently, split them.

**Practical guideline for this project**: everything going into `shared/lib/` deserves separate files by domain (e.g. `fastq_io.py`, `alignment_utils.py`, `plotting.py`), since it's meant to be imported selectively. Code inside `pipelines/{name}/scripts/` can stay more consolidated by step, since it's pipeline-specific and runs sequentially anyway. Don't over-split too early — refactor into separate files once functions start naturally clustering by theme (usually after the second pipeline reuses something).
