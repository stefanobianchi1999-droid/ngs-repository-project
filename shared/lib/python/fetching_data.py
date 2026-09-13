"""Utility to download text datasets from a URL and save them to disk."""

import argparse
from pathlib import Path

import requests

# Maps a short dataset name to its (url, default output path).
DATASETS = {
    "expression": (
        "https://schatz-lab.org/teaching/exercises/rnaseq/rnaseq.1.expression/expression.txt",
        "exercises/expression.txt",
    ),
}

DEFAULT_DATASET = "expression"


def fetch_dataset(url: str, output_path: str) -> Path:
    """Download the file at the given URL and save it to output_path."""
    response = requests.get(url, timeout=30)
    response.raise_for_status()

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_bytes(response.content)

    return output_file


def main() -> None:
    parser = argparse.ArgumentParser(description="Download a known dataset by name.")
    parser.add_argument(
        "dataset",
        nargs="?",
        default=DEFAULT_DATASET,
        choices=DATASETS.keys(),
        help="Name of the dataset to download",
    )
    parser.add_argument("--output", help="Override the destination file path")
    args = parser.parse_args()

    url, default_output = DATASETS[args.dataset]
    output_path = args.output or default_output

    saved_path = fetch_dataset(url, output_path)
    print(f"Dataset saved to: {saved_path}")


if __name__ == "__main__":
    main()
