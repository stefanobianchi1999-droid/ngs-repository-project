"""Utility to download text datasets from a URL and save them to disk."""

import argparse
from pathlib import Path

import requests


def fetch_dataset(url: str, output_path: str) -> Path:
    """Download the file at the given URL and save it to output_path."""
    response = requests.get(url, timeout=30)
    response.raise_for_status()

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_bytes(response.content)

    return output_file


def main() -> None:
    parser = argparse.ArgumentParser(description="Download a dataset from a URL.")
    parser.add_argument("url", help="URL of the file to download")
    parser.add_argument("output", help="Destination file path")
    args = parser.parse_args()

    saved_path = fetch_dataset(args.url, args.output)
    print(f"Dataset saved to: {saved_path}")


if __name__ == "__main__":
    main()
