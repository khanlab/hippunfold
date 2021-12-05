"""
console-script for snakebids app
"""
from pathlib import Path
from snakebids.app import SnakeBidsApp


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir"""
    app = SnakeBidsApp("../", skip_parse_args=True)
    return app.parser


def main():
    """Runs snakemake"""
    app = SnakeBidsApp(Path(__file__).resolve().parents[1])  # to get repository root
    app.run_snakemake()


if __name__ == "__main__":
    main()
