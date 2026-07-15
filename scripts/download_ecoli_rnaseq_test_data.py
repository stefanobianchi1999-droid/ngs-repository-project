#!/usr/bin/env python3
"""
download_ecoli_rnaseq_test_data.py

Downloads test data for an E. coli RNA-seq pipeline (FastQC -> fastp -> STAR/salmon
-> featureCounts -> DESeq2):

  1. Reference genome + GTF annotation: E. coli K-12 MG1655 (NCBI assembly
     GCF_000005845.2 / ASM584v2) -> shared/reference/ecoli/
  2. RNA-seq FASTQ reads (paired or single-end, .fastq.gz) from ENA, given
     one or more run accessions -> data/raw/rna-seq/

Only uses the standard library (urllib, gzip) -- no extra installs needed.

--- How to find RNA-seq run accessions for E. coli ---
1. Go to https://www.ebi.ac.uk/ena/browser/text-search?query=Escherichia%20coli%20RNA-Seq
   (or search directly on https://www.ebi.ac.uk/ena/browser/)
2. Filter results by "Library strategy: RNA-Seq" and pick a study/run that matches
   your organism/condition of interest. Prefer smaller runs for a first test.
3. Note the Run accession (starts with SRR/ERR/DRR, e.g. "SRR1234567").
4. Pass it to this script:  python download_ecoli_rnaseq_test_data.py SRR1234567 SRR1234568

Usage:
    python scripts/download_ecoli_rnaseq_test_data.py [RUN_ACCESSION ...] [--force]

If no accession is given, only the genome + GTF are downloaded.
Already-downloaded files are skipped unless --force is passed.
Run this from the root of the project (see sequencing-project-structure.md).
"""

import argparse
import json
import shutil
import sys
import urllib.error
import urllib.request
from pathlib import Path

GENOME_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
    "GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
)
GTF_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
    "GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gtf.gz"
)
ENA_FILEREPORT_API = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={accession}&result=read_run&fields=fastq_ftp&format=json"
)

REF_DIR = Path("shared/reference/ecoli")
FASTQ_DIR = Path("data/raw/rna-seq")

# Network timeout for a single read/connect operation, in seconds.
TIMEOUT = 60


def _human(size: int) -> str:
    """Format a byte count as a short human-readable string."""
    value = float(size)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024 or unit == "TB":
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} TB"


def download(url: str, dest: Path, force: bool = False) -> None:
    """Download `url` to `dest`, atomically and idempotently.

    Skips the download if `dest` already exists (unless `force`). Writes to a
    temporary `.part` file first and renames on success, so an interrupted
    download never leaves a truncated file that looks complete.
    """
    if dest.exists() and not force:
        print(f"  [skip] {dest} already exists ({_human(dest.stat().st_size)})")
        return

    print(f"Downloading {url}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with urllib.request.urlopen(url, timeout=TIMEOUT) as response, open(tmp, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
    except Exception:
        tmp.unlink(missing_ok=True)  # don't leave a partial file behind
        raise
    tmp.replace(dest)
    print(f"  -> saved to {dest} ({_human(dest.stat().st_size)})")


def download_genome_and_annotation(force: bool = False) -> None:
    print("== Downloading E. coli K-12 MG1655 genome + GTF annotation ==")
    download(GENOME_URL, REF_DIR / "ecoli_k12_mg1655.fasta.gz", force)
    download(GTF_URL, REF_DIR / "ecoli_k12_mg1655.gtf.gz", force)
    print("Note: files are kept gzipped (.gz) — most aligners (STAR, salmon) accept")
    print("gzipped FASTA/annotation directly, or gunzip -k if a tool needs plain text.")


def fetch_ena_fastq_urls(accession: str) -> list[str]:
    url = ENA_FILEREPORT_API.format(accession=accession)
    with urllib.request.urlopen(url, timeout=TIMEOUT) as response:
        data = json.load(response)
    if not data:
        raise ValueError(f"No ENA record found for accession {accession}")
    fastq_ftp = data[0].get("fastq_ftp", "")
    if not fastq_ftp:
        raise ValueError(f"No FASTQ files listed for accession {accession}")
    return [f"https://{u}" if not u.startswith("http") else u for u in fastq_ftp.split(";")]


def download_rnaseq_fastq(accessions: list[str], force: bool = False) -> int:
    """Download FASTQ files for the given accessions. Returns the failure count."""
    print(f"== Downloading RNA-seq FASTQ for {len(accessions)} accession(s) from ENA ==")
    failures = 0
    for accession in accessions:
        try:
            urls = fetch_ena_fastq_urls(accession)
        except Exception as exc:
            print(f"  [!] Skipping {accession}: {exc}")
            failures += 1
            continue
        for url in urls:
            dest = FASTQ_DIR / Path(url).name  # already ends in .fastq.gz
            try:
                download(url, dest, force)
            except Exception as exc:
                print(f"  [!] Failed to download {url}: {exc}")
                failures += 1
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download E. coli genome/annotation and optional RNA-seq FASTQ runs.",
    )
    parser.add_argument(
        "accessions", nargs="*", metavar="RUN_ACCESSION",
        help="ENA/SRA run accessions (SRR/ERR/DRR...). If omitted, only the genome is fetched.",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-download files even if they already exist locally.",
    )
    args = parser.parse_args()

    failures = 0
    try:
        download_genome_and_annotation(args.force)
    except (urllib.error.URLError, OSError) as exc:
        print(f"  [!] Failed to download genome/annotation: {exc}")
        failures += 1

    if args.accessions:
        failures += download_rnaseq_fastq(args.accessions, args.force)
    else:
        print()
        print("No run accession given — skipping FASTQ download.")
        print("Find one on https://www.ebi.ac.uk/ena/browser/ (filter: Library strategy = RNA-Seq)")
        print("then re-run: python scripts/download_ecoli_rnaseq_test_data.py SRRxxxxxxx")

    print()
    if failures:
        print(f"Done with {failures} failure(s) — see [!] messages above.")
    else:
        print("Done. Try it out once you have FASTQ files:")
        print(f"  fastqc {FASTQ_DIR}/*.fastq.gz -o results/rna-seq/fastqc_raw/")
        print(
            "  fastp -i <sample>_1.fastq.gz -I <sample>_2.fastq.gz \\\n"
            "        -o <sample>_1.clean.fastq.gz -O <sample>_2.clean.fastq.gz \\\n"
            "        -h <sample>_fastp.html -j <sample>_fastp.json"
        )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
