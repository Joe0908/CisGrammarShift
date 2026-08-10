from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def download_asset(record: dict[str, object], output_directory: Path) -> dict[str, object]:
    destination = output_directory / str(record["filename"])
    expected_size = int(record["size_bytes"]) if record.get("size_bytes") else None
    expected_sha256 = str(record["sha256"]) if record.get("sha256") else None
    complete = destination.exists() and (expected_size is None or destination.stat().st_size == expected_size)
    if complete and expected_sha256 is not None:
        complete = sha256(destination) == expected_sha256

    if not complete:
        temporary = destination.with_suffix(destination.suffix + ".part")
        if destination.exists():
            destination.replace(temporary)
        resume_from = temporary.stat().st_size if temporary.exists() else 0
        if expected_size is not None and resume_from > expected_size:
            temporary.unlink()
            resume_from = 0
        request = urllib.request.Request(str(record["url"]))
        if resume_from:
            request.add_header("Range", f"bytes={resume_from}-")
        with urllib.request.urlopen(request) as response:
            partial_response = getattr(response, "status", None) == 206
            mode = "ab" if resume_from and partial_response else "wb"
            with temporary.open(mode) as handle:
                shutil.copyfileobj(response, handle, length=1024 * 1024)
        if expected_size is not None and temporary.stat().st_size != expected_size:
            raise RuntimeError(f"size mismatch for partial download {temporary}")
        temporary.replace(destination)

    observed = sha256(destination)
    if expected_sha256 is not None and observed != expected_sha256:
        raise RuntimeError(f"checksum mismatch for {destination}")
    if expected_size is not None and destination.stat().st_size != expected_size:
        raise RuntimeError(f"size mismatch for {destination}")
    return {
        **record,
        "size_bytes": destination.stat().st_size,
        "sha256": observed,
        "sha256_source": "manifest" if expected_sha256 else "computed_after_download",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Download only whitelisted public assets with checksums")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--resolved-manifest", type=Path)
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    records = json.loads(args.manifest.read_text(encoding="utf-8"))
    args.output_directory.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [
            executor.submit(download_asset, record, args.output_directory)
            for record in records["assets"]
        ]
        resolved = []
        for future in futures:
            resolved_record = future.result()
            resolved.append(resolved_record)
            destination = args.output_directory / str(resolved_record["filename"])
            print(
                f"{destination}\t{resolved_record['size_bytes']}\t{resolved_record['sha256']}",
                flush=True,
            )
    if args.resolved_manifest is not None:
        args.resolved_manifest.parent.mkdir(parents=True, exist_ok=True)
        args.resolved_manifest.write_text(
            json.dumps({**records, "assets": resolved}, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
