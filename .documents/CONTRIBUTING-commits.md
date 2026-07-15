# Commit Message Guidelines

This project follows the **Conventional Commits** standard, adapted for an RNA-seq / NGS bioinformatics pipeline, to keep the history clear, consistent, and easy to automate (changelogs, release notes, CI triggers, etc.).

## Structure

```
<type>(<optional scope>): <short description>

<optional body>

<optional footer>
```

## Available Types

| Type       | When to use it                                                        |
|------------|--------------------------------------------------------------------------|
| `feat`     | New feature (new analysis step, new module, new CLI option)             |
| `fix`      | Bug fix (wrong parameter, incorrect output, crash on edge case)         |
| `docs`     | Documentation, README, usage examples, workflow diagrams                |
| `style`    | Formatting, linting, whitespace (no logic change)                       |
| `refactor` | Code restructuring without changing scientific output                  |
| `test`     | Adding or updating unit/integration tests, test datasets                |
| `perf`     | Performance improvements (speed, memory, parallelization)               |
| `chore`    | Maintenance: dependencies, `.gitignore`, CI config, containers          |
| `data`     | Changes to reference data (genome, annotation, indexes)                 |

## Suggested Scopes

Since this is an NGS/RNA-seq tool, scopes can reflect pipeline stages:

- `qc` – quality control (FastQC, MultiQC, trimming)
- `align` – alignment/mapping (STAR, HISAT2, Salmon, Bowtie2)
- `counts` – read counting / quantification (featureCounts, HTSeq, tximport)
- `dea` – differential expression analysis (DESeq2, edgeR, limma)
- `annotation` – gene/transcript annotation, GTF/GFF handling
- `pipeline` – workflow orchestration (Nextflow, Snakemake, CWL)
- `report` – output reports, plots, MultiQC summaries
- `config` – config files, sample sheets, parameter files
- `docker` / `env` – containers, conda environments, dependencies

## Examples

```
feat(align): add support for STAR two-pass mode

fix(counts): correct strandedness parameter for featureCounts

docs(pipeline): add Nextflow usage example for paired-end reads

refactor(dea): simplify DESeq2 contrast generation logic

perf(align): enable multi-threaded alignment with STAR

chore(env): update conda environment.yml with samtools 1.20

data(annotation): update GTF to GENCODE v46

test(qc): add unit test for adapter trimming step

fix(pipeline): resolve crash on samples with zero mapped reads
```

## Breaking Changes

If the commit changes the pipeline's input/output format or CLI in an incompatible way:

```
feat(pipeline)!: switch count matrix output to CSV instead of TSV
```

or, in the footer:

```
BREAKING CHANGE: output count matrices are now comma-separated (.csv) instead of tab-separated (.tsv); downstream scripts must be updated
```

## Best Practices

1. **Short subject line** (max ~50-72 characters), imperative mood: `add`, not `added`
2. **Blank line** between subject and body
3. The **body explains why**, especially for changes affecting scientific results (e.g., changing a statistical threshold, aligner parameter, or normalization method) — reviewers and future users need to know the rationale
4. **Reference issues/tickets** in the footer: `Closes #42`, `Refs #17`
5. Keep commits **atomic**: don't mix pipeline logic changes with reference data updates in the same commit
6. If a commit changes results reproducibility (e.g., updates a reference genome, tool version, or random seed), mention it explicitly in the body

## Useful Commands

Short commit:
```bash
git commit -m "fix(align): correct read group tag in STAR output"
```

Commit with extended body (opens the editor) — recommended when a change affects scientific results:
```bash
git commit
```