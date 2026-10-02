#!/usr/bin/env python3
"""Compatibility adapter for SEPP pplacer jobs using EPA-ng.

PICRUSt2's bundled SEPP invokes pplacer with a reference alignment, a
RAxML-info file, a reference subtree, and a query-only aligned FASTA.  The
bundled pplacer binary segfaults on the current Linux runtime before reading
the inputs.  This adapter preserves SEPP's decomposition/alignment/merging
workflow but performs each subtree placement with the official EPA-ng binary
and the matching PICRUSt2 ``*.model`` file.

Only the pplacer arguments emitted by SEPP 4.5.5 are accepted.  Unknown
arguments cause a hard failure so that an interface change cannot silently
alter the analysis.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def read_fasta_ids(path: Path) -> list[str]:
    identifiers = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].strip().split()[0]
                if not identifier:
                    raise ValueError(f"Empty FASTA identifier in {path}")
                identifiers.append(identifier)
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"Duplicate query identifiers in {path}")
    return identifiers


def validate_jplace(path: Path, expected_ids: list[str]) -> None:
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    required_fields = {
        "edge_num",
        "likelihood",
        "like_weight_ratio",
        "distal_length",
        "pendant_length",
    }
    if not required_fields.issubset(set(payload.get("fields", []))):
        raise ValueError(f"Unexpected jplace fields in {path}: {payload.get('fields')}")
    observed_ids: list[str] = []
    for placement in payload.get("placements", []):
        if "n" in placement:
            observed_ids.extend(str(value) for value in placement["n"])
        elif "nm" in placement:
            observed_ids.extend(str(value[0]) for value in placement["nm"])
        else:
            raise ValueError(f"Placement lacks n/nm names in {path}")
    if len(observed_ids) != len(set(observed_ids)):
        raise ValueError(f"Duplicate placement identifiers in {path}")
    if set(observed_ids) != set(expected_ids):
        missing = sorted(set(expected_ids) - set(observed_ids))
        extra = sorted(set(observed_ids) - set(expected_ids))
        raise ValueError(
            f"Query/placement mismatch in {path}: "
            f"missing={missing[:8]}, extra={extra[:8]}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("-j", type=int, default=1)
    parser.add_argument("-r", dest="reference_alignment", required=True)
    parser.add_argument("-s", dest="raxml_info", required=True)
    parser.add_argument("-t", dest="reference_tree", required=True)
    parser.add_argument("--groups")
    parser.add_argument("query_alignment")
    args, unknown = parser.parse_known_args()
    if unknown:
        parser.error(f"unsupported pplacer arguments: {' '.join(unknown)}")
    return args


def main() -> int:
    args = parse_args()

    script_path = Path(__file__).resolve()
    project_root = script_path.parents[2]
    bundled_epa_ng = project_root / "tools" / "micromamba_root" / "bin" / "epa-ng"
    epa_ng_on_path = shutil.which("epa-ng")
    epa_ng = Path(epa_ng_on_path).resolve() if epa_ng_on_path else bundled_epa_ng

    raxml_info = Path(args.raxml_info).resolve()
    if not raxml_info.name.endswith(".raxml_info"):
        raise ValueError(
            f"Expected a PICRUSt2 *.raxml_info path, received: {raxml_info}"
        )
    model_file = raxml_info.with_name(
        raxml_info.name.removesuffix(".raxml_info") + ".model"
    )

    required = {
        "EPA-ng executable": epa_ng,
        "reference alignment": Path(args.reference_alignment),
        "query alignment": Path(args.query_alignment),
        "reference tree": Path(args.reference_tree),
        "EPA-ng model": model_file,
    }
    missing = [f"{label}: {path}" for label, path in required.items() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required input(s):\n" + "\n".join(missing))

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    query_alignment = Path(args.query_alignment).resolve()
    query_ids = read_fasta_ids(query_alignment)
    expected_jplace = Path(str(query_alignment).replace("fasta", "jplace"))

    temp_dir = Path(
        tempfile.mkdtemp(
            prefix=f".epa_ng_compat_{query_alignment.stem}_{os.getpid()}_",
            dir=out_dir,
        )
    )
    command = [
        str(epa_ng),
        "-s",
        str(Path(args.reference_alignment).resolve()),
        "-q",
        str(query_alignment),
        "-t",
        str(Path(args.reference_tree).resolve()),
        "-m",
        str(model_file),
        "--out-dir",
        str(temp_dir),
        "--threads",
        str(max(1, args.j)),
        "--filter-min-lwr",
        "0",
        "--filter-max",
        "10",
        "--redo",
    ]

    print(
        "SEPP placement compatibility adapter: "
        f"EPA-ng subtree placement for {query_alignment.name}",
        flush=True,
    )
    lock_path = Path(
        os.environ.get(
            "PICRUST2_EPA_COMPAT_LOCK",
            "/tmp/picrust2_sepp_epa_ng_compat.lock",
        )
    )
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        # SEPP can schedule one placement job per worker.  A whole-reference
        # EPA-ng run is memory intensive even after SEPP decomposition, so the
        # adapter serializes the placement executables while retaining
        # parallel HMM search/alignment upstream.
        with lock_path.open("w", encoding="utf-8") as lock_handle:
            print(f"Waiting for EPA-ng memory lock: {lock_path}", flush=True)
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
            print(f"Acquired EPA-ng memory lock: {lock_path}", flush=True)
            subprocess.run(command, check=True)
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
        generated = temp_dir / "epa_result.jplace"
        if not generated.exists() or generated.stat().st_size == 0:
            raise RuntimeError(f"EPA-ng did not create a valid jplace file: {generated}")
        validate_jplace(generated, query_ids)
        shutil.move(str(generated), str(expected_jplace))
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    print(
        f"Placement written to {expected_jplace}; "
        f"validated {len(query_ids)} query identifiers",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
