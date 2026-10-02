#!/usr/bin/env python3
"""Resume the split-reference PICRUSt2 pipeline from validated placement trees.

This utility is only a runtime-recovery helper.  It skips one or both
``place_seqs.py`` calls and copies trees produced previously from the exact
same FASTA, reference release, placement settings, and random seed.  Every
downstream PICRUSt2 step runs normally.

For a clean reproduction without a cache, use ``03_run_picrust2_primary.sh``.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import runpy
import shutil
import sys
import time

import picrust2.pipeline as pipeline


def option_value(command: list[str], option: str) -> str:
    try:
        return command[command.index(option) + 1]
    except (ValueError, IndexError) as error:
        raise ValueError(f"Expected {option} in PICRUSt2 command: {command}") from error


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--cached-bac-tree", required=True, type=Path)
    parser.add_argument("--cached-arc-tree", type=Path)
    parser.add_argument("--cache-record", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    cache_args, pipeline_args = parser.parse_known_args()

    cached_tree = cache_args.cached_bac_tree.resolve()
    if not cached_tree.is_file() or cached_tree.stat().st_size == 0:
        raise FileNotFoundError(f"Validated bacterial placement tree not found: {cached_tree}")
    cached_arc_tree = (
        cache_args.cached_arc_tree.resolve()
        if cache_args.cached_arc_tree is not None
        else None
    )
    if cached_arc_tree is not None and (
        not cached_arc_tree.is_file() or cached_arc_tree.stat().st_size == 0
    ):
        raise FileNotFoundError(
            f"Validated archaeal placement tree not found: {cached_arc_tree}"
        )

    original_system_call = pipeline.system_call_check
    used_domains: set[str] = set()

    def system_call_with_cache(command, *args, **kwargs):
        nonlocal used_domains
        if (
            isinstance(command, list)
            and command
            and Path(command[0]).name == "place_seqs.py"
            and Path(option_value(command, "--ref_dir")).name in {"bac_ref", "arc_ref"}
        ):
            reference_name = Path(option_value(command, "--ref_dir")).name
            domain = "bacterial" if reference_name == "bac_ref" else "archaeal"
            source_tree = cached_tree if domain == "bacterial" else cached_arc_tree
            if source_tree is None:
                return original_system_call(command, *args, **kwargs)
            if domain in used_domains:
                raise RuntimeError(f"The {domain} placement cache was requested more than once.")
            destination = Path(option_value(command, "--out_tree")).resolve()
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source_tree, destination)
            intermediate = Path(option_value(command, "--intermediate")).resolve()
            intermediate.mkdir(parents=True, exist_ok=True)
            (intermediate / f"CACHED_{domain.upper()}_PLACEMENT.txt").write_text(
                f"The {domain} tree was copied from a validated placement run made "
                "with the identical study FASTA, PICRUSt2-SC reference, and "
                "min_align 0.8. See the external cache record and domain-specific "
                "placement audit for checksums and commands.\n",
                encoding="utf-8",
            )
            used_domains.add(domain)
            print(
                f"Reused validated {domain} placement tree: {source_tree} -> {destination}",
                file=sys.stderr,
                flush=True,
            )
            return None
        return original_system_call(command, *args, **kwargs)

    pipeline.system_call_check = system_call_with_cache
    picrust_entry = shutil.which("picrust2_pipeline.py")
    if not picrust_entry:
        raise FileNotFoundError("picrust2_pipeline.py is not available on PATH.")
    sys.argv = [picrust_entry, *pipeline_args]
    started = time.monotonic()
    runpy.run_path(picrust_entry, run_name="__main__")
    elapsed_seconds = time.monotonic() - started
    expected_cached_domains = {"bacterial"}
    if cached_arc_tree is not None:
        expected_cached_domains.add("archaeal")
    if used_domains != expected_cached_domains:
        raise RuntimeError("PICRUSt2 completed without consuming the bacterial tree cache.")

    cache_args.cache_record.parent.mkdir(parents=True, exist_ok=True)
    study_fasta = Path(option_value(pipeline_args, "-s")).resolve()
    archaeal_record = (
        f"cached_archaeal_tree={cached_arc_tree}\n"
        f"cached_archaeal_tree_sha256={sha256_file(cached_arc_tree)}\n"
        if cached_arc_tree is not None
        else ""
    )
    cache_args.cache_record.write_text(
        f"cached_bacterial_tree={cached_tree}\n"
        f"cached_bacterial_tree_size={cached_tree.stat().st_size}\n"
        f"cached_bacterial_tree_sha256={sha256_file(cached_tree)}\n"
        f"{archaeal_record}"
        f"study_fasta={study_fasta}\n"
        f"study_fasta_sha256={sha256_file(study_fasta)}\n"
        "downstream_pipeline_args=" + " ".join(pipeline_args) + "\n",
        encoding="utf-8",
    )
    abundance_table = Path(option_value(pipeline_args, "-i")).resolve()
    cache_args.run_manifest.parent.mkdir(parents=True, exist_ok=True)
    cache_args.run_manifest.write_text(
        "PICRUSt2 primary run with validated placement caches\n"
        "====================================================\n"
        "The cached trees are completed placements from the exact same study "
        "FASTA, PICRUSt2-SC references, and min_align threshold. Bacterial "
        "placement used SEPP decomposition plus EPA-ng subtree placement; "
        "archaeal placement used PICRUSt2's whole-reference EPA-ng path. Only "
        "the redundant place_seqs.py calls were skipped; every downstream "
        "PICRUSt2 step ran normally.\n\n"
        f"elapsed_seconds={elapsed_seconds:.3f}\n"
        f"study_fasta={study_fasta}\n"
        f"study_fasta_sha256={sha256_file(study_fasta)}\n"
        f"abundance_table={abundance_table}\n"
        f"abundance_table_sha256={sha256_file(abundance_table)}\n"
        f"cached_bacterial_tree={cached_tree}\n"
        f"cached_bacterial_tree_sha256={sha256_file(cached_tree)}\n"
        f"{archaeal_record}"
        "pipeline_args=" + " ".join(pipeline_args) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
