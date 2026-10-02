#!/usr/bin/env python3
"""Finish SEPP after HMM stages when only subtree placement remains.

The script is deliberately strict.  It reconstructs SEPP's deterministic
reference-tree decomposition from the original reference files, verifies that
every reconstructed placement subset has the same leaf set and topology as
the retained SEPP subtree, maps each query alignment to its matching reference
alignment, runs the EPA-ng compatibility adapter at bounded concurrency,
validates every jplace file, invokes SEPP's standard JSON merger, validates the
merged query set, and finally grafts placements with gappa.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

from pplacer_epa_ng_compat import read_fasta_ids, validate_jplace

sys.setrecursionlimit(200_000)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def alignment_columns(path: Path) -> int:
    sequence_parts: list[str] = []
    started = False
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if started:
                    break
                started = True
            elif started:
                sequence_parts.append(line.strip())
    if not sequence_parts:
        raise ValueError(f"No sequence found in alignment: {path}")
    return len("".join(sequence_parts))


def chunk_number(path: Path) -> int:
    match = re.search(r"pplacer\.extended\.(\d+)\.", path.name)
    if not match:
        raise ValueError(f"Cannot parse SEPP fragment-chunk number: {path}")
    return int(match.group(1))


def map_query_backbones(problem_dir: Path) -> list[tuple[int, Path, Path]]:
    queries = sorted(
        problem_dir.glob("pplacer.extended.*.fasta"),
        key=chunk_number,
    )
    backbones = list(problem_dir.glob("pplacer.backbone.*.fasta"))
    if len(queries) != len(backbones):
        raise ValueError(
            f"{problem_dir.name}: {len(queries)} query alignments but "
            f"{len(backbones)} reference alignments"
        )
    backbone_columns = {path: alignment_columns(path) for path in backbones}
    unused = set(backbones)
    mapping = []
    for query in queries:
        columns = alignment_columns(query)
        candidates = [path for path in unused if backbone_columns[path] == columns]
        if not candidates:
            raise ValueError(
                f"{problem_dir.name}: no {columns}-column backbone for {query.name}"
            )
        # Files are written query-then-backbone by SEPP's JoinAlignJobs.
        # Inode proximity resolves the rare case of equal alignment lengths.
        backbone = min(
            candidates,
            key=lambda path: abs(path.stat().st_ino - query.stat().st_ino),
        )
        unused.remove(backbone)
        mapping.append((chunk_number(query), query, backbone))
    if unused:
        raise ValueError(f"{problem_dir.name}: unassigned backbones remain")
    return mapping


def parse_tree(path: Path):
    import dendropy

    return dendropy.Tree.get(
        path=str(path),
        schema="newick",
        preserve_underscores=True,
    )


def tree_leaf_names(tree) -> set[str]:
    return {node.taxon.label for node in tree.leaf_node_iter()}


def verify_subtree(expected_problem, observed_tree_path: Path) -> dict:
    import dendropy
    from dendropy.calculate import treecompare

    expected_newick = expected_problem.subtree.compose_newick(labels=False) + ";"
    namespace = dendropy.TaxonNamespace()
    expected_tree = dendropy.Tree.get(
        data=expected_newick,
        schema="newick",
        taxon_namespace=namespace,
        preserve_underscores=True,
    )
    observed_tree = dendropy.Tree.get(
        path=str(observed_tree_path),
        schema="newick",
        taxon_namespace=namespace,
        preserve_underscores=True,
    )
    expected_tree.encode_bipartitions()
    observed_tree.encode_bipartitions()
    symmetric_difference = int(
        treecompare.symmetric_difference(expected_tree, observed_tree)
    )
    expected_leaves = tree_leaf_names(expected_tree)
    observed_leaves = tree_leaf_names(observed_tree)
    if expected_leaves != observed_leaves or symmetric_difference != 0:
        raise ValueError(
            f"{expected_problem.label}: reconstructed/retained subtree mismatch; "
            f"leaf_difference={len(expected_leaves ^ observed_leaves)}, "
            f"symmetric_difference={symmetric_difference}"
        )
    return {
        "reference_leaves": len(expected_leaves),
        "symmetric_difference": symmetric_difference,
    }


def run_placement(
    adapter: Path,
    raxml_info: Path,
    tree: Path,
    chunk: int,
    query: Path,
    backbone: Path,
    log_dir: Path,
) -> dict:
    expected_jplace = Path(str(query).replace("fasta", "jplace"))
    query_ids = read_fasta_ids(query)
    if expected_jplace.exists() and expected_jplace.stat().st_size > 0:
        validate_jplace(expected_jplace, query_ids)
        status = "validated_existing"
    else:
        command = [
            str(adapter),
            "--out-dir",
            str(query.parent),
            "-j",
            "1",
            "-r",
            str(backbone),
            "-s",
            str(raxml_info),
            "-t",
            str(tree),
            "--groups",
            "10",
            str(query),
        ]
        log_path = log_dir / f"{query.parent.name}_chunk_{chunk}.log"
        with log_path.open("w", encoding="utf-8") as log_handle:
            subprocess.run(
                command,
                check=True,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
            )
        validate_jplace(expected_jplace, query_ids)
        status = "completed"
    return {
        "problem": query.parent.name,
        "chunk": chunk,
        "queries": len(query_ids),
        "alignment_columns": alignment_columns(query),
        "query_alignment": str(query),
        "reference_alignment": str(backbone),
        "jplace": str(expected_jplace),
        "jplace_sha256": sha256_file(expected_jplace),
        "status": status,
    }


def reconstruct_decomposition(args: argparse.Namespace):
    reconstruction_out = Path(
        tempfile.mkdtemp(prefix="sepp_decomposition_reconstruction_")
    )
    sys.argv = [
        "run_sepp.py",
        "--tree",
        str(args.reference_tree),
        "--raxml",
        str(args.raxml_info),
        "--cpu",
        str(args.sepp_cpu),
        "--molecule",
        "dna",
        "--outdir",
        str(reconstruction_out),
        "-seed",
        str(args.seed),
        "--alignment",
        str(args.reference_alignment),
        "--fragment",
        str(args.filtered_study_fasta),
    ]
    from sepp.exhaustive import ExhaustiveAlgorithm

    algorithm = ExhaustiveAlgorithm()
    algorithm.check_options()
    algorithm.root_problem = algorithm.build_subproblems()
    return algorithm.root_problem, reconstruction_out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sepp-temp-root", required=True, type=Path)
    parser.add_argument("--filtered-study-fasta", required=True, type=Path)
    parser.add_argument("--reference-alignment", required=True, type=Path)
    parser.add_argument("--reference-tree", required=True, type=Path)
    parser.add_argument("--raxml-info", required=True, type=Path)
    parser.add_argument("--adapter", required=True, type=Path)
    parser.add_argument("--json-merger", required=True, type=Path)
    parser.add_argument("--sepp-output-dir", required=True, type=Path)
    parser.add_argument("--out-tree", required=True, type=Path)
    parser.add_argument("--log-dir", required=True, type=Path)
    parser.add_argument("--audit-json", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--sepp-cpu", type=int, default=6)
    parser.add_argument("--seed", type=int, default=297834)
    args = parser.parse_args()

    for path in [
        args.sepp_temp_root,
        args.filtered_study_fasta,
        args.reference_alignment,
        args.reference_tree,
        args.raxml_info,
        args.adapter,
        args.json_merger,
    ]:
        if not path.exists():
            raise FileNotFoundError(path)
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")

    args.log_dir.mkdir(parents=True, exist_ok=True)
    args.sepp_output_dir.mkdir(parents=True, exist_ok=True)
    args.audit_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_tree.parent.mkdir(parents=True, exist_ok=True)

    root_problem, reconstruction_out = reconstruct_decomposition(args)
    reconstructed = {problem.label: problem for problem in root_problem.get_children()}
    retained_dirs = {
        path.name: path
        for path in args.sepp_temp_root.glob("P_*")
        if path.is_dir()
    }
    if set(reconstructed) != set(retained_dirs):
        raise ValueError(
            "Reconstructed and retained SEPP placement subsets differ: "
            f"reconstructed={sorted(reconstructed)}, retained={sorted(retained_dirs)}"
        )

    all_jobs = []
    subset_audit = {}
    expected_query_ids = read_fasta_ids(args.filtered_study_fasta)
    observed_query_ids: list[str] = []
    for label in sorted(reconstructed, key=lambda value: int(value.split("_")[1])):
        problem_dir = retained_dirs[label]
        tree_files = list(problem_dir.glob("pplacer.tree.*.tre"))
        if not tree_files:
            raise ValueError(f"No retained subtree for {label}")
        tree_hashes = {sha256_file(path) for path in tree_files}
        if len(tree_hashes) != 1:
            raise ValueError(f"Retained subtree copies differ within {label}")
        tree_path = tree_files[0]
        subset_audit[label] = verify_subtree(reconstructed[label], tree_path)
        mapping = map_query_backbones(problem_dir)
        subset_audit[label]["chunks"] = len(mapping)
        subset_audit[label]["tree_sha256"] = next(iter(tree_hashes))
        for chunk, query, backbone in mapping:
            query_ids = read_fasta_ids(query)
            observed_query_ids.extend(query_ids)
            all_jobs.append(
                (
                    args.adapter.resolve(),
                    args.raxml_info.resolve(),
                    tree_path.resolve(),
                    chunk,
                    query.resolve(),
                    backbone.resolve(),
                    args.log_dir.resolve(),
                )
            )
    if len(observed_query_ids) != len(set(observed_query_ids)):
        raise ValueError("A query ASV occurs in more than one SEPP placement job.")
    if set(observed_query_ids) != set(expected_query_ids):
        raise ValueError(
            "Retained SEPP placement jobs do not cover the filtered study FASTA: "
            f"missing={len(set(expected_query_ids) - set(observed_query_ids))}, "
            f"extra={len(set(observed_query_ids) - set(expected_query_ids))}"
        )

    placement_rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(run_placement, *job) for job in all_jobs]
        for future in as_completed(futures):
            placement_rows.append(future.result())
            completed = len(placement_rows)
            print(
                f"Validated EPA-ng subtree placements: {completed}/{len(all_jobs)}",
                flush=True,
            )
    placement_rows.sort(key=lambda row: (int(row["problem"].split("_")[1]), row["chunk"]))

    merge_lines = [root_problem.subtree.compose_newick(labels=True) + ";"]
    by_key = {(row["problem"], row["chunk"]): row for row in placement_rows}
    for label in sorted(reconstructed, key=lambda value: int(value.split("_")[1])):
        problem = reconstructed[label]
        for chunk in range(root_problem.fragment_chunks):
            row = by_key[(label, chunk)]
            merge_lines.append(
                problem.subtree.compose_newick(labels=True)
                + ";\n"
                + row["jplace"]
            )
    merge_lines.extend(["", ""])
    merge_input = "\n".join(merge_lines)
    merged_jplace = args.sepp_output_dir / "output_placement.json"
    merger_log = args.log_dir / "sepp_json_merger.log"
    with merger_log.open("w", encoding="utf-8") as log_handle:
        subprocess.run(
            ["java", "-jar", str(args.json_merger), "-", "-", str(merged_jplace)],
            input=merge_input,
            text=True,
            check=True,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
        )
    validate_jplace(merged_jplace, expected_query_ids)

    gappa_log = args.log_dir / "gappa_graft.log"
    with gappa_log.open("w", encoding="utf-8") as log_handle:
        subprocess.run(
            [
                "gappa",
                "examine",
                "graft",
                "--jplace-path",
                str(merged_jplace),
                "--fully-resolve",
                "--out-dir",
                str(args.sepp_output_dir),
            ],
            check=True,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    generated_tree = merged_jplace.with_suffix(".newick")
    if not generated_tree.exists() or generated_tree.stat().st_size == 0:
        raise RuntimeError(f"gappa did not produce the expected tree: {generated_tree}")
    shutil.copy2(generated_tree, args.out_tree)
    output_leaf_names = tree_leaf_names(parse_tree(args.out_tree))
    reference_leaf_names = tree_leaf_names(parse_tree(args.reference_tree))
    expected_output_leaves = reference_leaf_names | set(expected_query_ids)
    if output_leaf_names != expected_output_leaves:
        raise ValueError(
            "Grafted-tree leaf audit failed: "
            f"missing={len(expected_output_leaves - output_leaf_names)}, "
            f"extra={len(output_leaf_names - expected_output_leaves)}"
        )

    audit = {
        "status": "passed",
        "reason_for_resume": (
            "The six-process SEPP run completed HMM search and alignment but "
            "six simultaneous EPA-ng subtree placements exceeded the 14-GiB "
            "memory limit. Retained HMM outputs were completed at bounded "
            "EPA-ng concurrency without recomputing or changing alignments."
        ),
        "seed": args.seed,
        "sepp_cpu_used_for_retained_hmm_stage": args.sepp_cpu,
        "epa_ng_max_concurrent_jobs": args.workers,
        "placement_subsets": len(reconstructed),
        "placement_jobs": len(placement_rows),
        "validated_query_ASVs": len(expected_query_ids),
        "validated_output_tree_leaves": len(output_leaf_names),
        "subsets": subset_audit,
        "placements": placement_rows,
        "merged_jplace": str(merged_jplace),
        "merged_jplace_sha256": sha256_file(merged_jplace),
        "output_tree": str(args.out_tree),
        "output_tree_sha256": sha256_file(args.out_tree),
        "reconstruction_temporary_directory": str(reconstruction_out),
    }
    args.audit_json.write_text(
        json.dumps(audit, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps({key: audit[key] for key in [
        "status",
        "placement_subsets",
        "placement_jobs",
        "validated_query_ASVs",
        "merged_jplace_sha256",
        "output_tree_sha256",
    ]}, indent=2))


if __name__ == "__main__":
    main()
