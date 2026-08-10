from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description="Download only whitelisted public assets with checksums")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    records = json.loads(args.manifest.read_text(encoding="utf-8"))
    args.output_directory.mkdir(parents=True, exist_ok=True)
    for record in records["assets"]:
        destination = args.output_directory / record["filename"]
        if not destination.exists() or (record.get("sha256") and sha256(destination) != record["sha256"]):
            temporary = destination.with_suffix(destination.suffix + ".part")
            urllib.request.urlretrieve(record["url"], temporary)
            temporary.replace(destination)
        observed = sha256(destination)
        if record.get("sha256") and observed != record["sha256"]:
            raise RuntimeError(f"checksum mismatch for {destination}")
        if record.get("size_bytes") and destination.stat().st_size != record["size_bytes"]:
            raise RuntimeError(f"size mismatch for {destination}")
        print(f"{destination}\t{destination.stat().st_size}\t{observed}")


if __name__ == "__main__":
    main()
