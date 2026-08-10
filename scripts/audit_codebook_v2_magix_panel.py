from __future__ import annotations

import argparse
from pathlib import Path

from cisgrammar.capselex import write_json
from cisgrammar.capselex_codebook_v2 import audit_codebook_v2_focal_panel


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit official Codebook v2 focal MAGIX BED files")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--data-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write_json(
        audit_codebook_v2_focal_panel(args.config, args.data_directory),
        args.output,
    )


if __name__ == "__main__":
    main()
