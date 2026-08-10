from __future__ import annotations

import argparse
import json
from pathlib import Path

from cisgrammar.capselex import write_json
from cisgrammar.capselex_codebook_chip import build_gpzn_manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Resolve focal Codebook GPZN bigWigs from GEO metadata")
    parser.add_argument("--soft", type=Path, required=True)
    parser.add_argument("--filelist", type=Path, required=True)
    parser.add_argument("--panel-config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    panel = json.loads(args.panel_config.read_text(encoding="utf-8"))["panel"]
    focal_tfs = [record["tf"] for record in panel]
    write_json(build_gpzn_manifest(args.soft, args.filelist, focal_tfs), args.output)


if __name__ == "__main__":
    main()
