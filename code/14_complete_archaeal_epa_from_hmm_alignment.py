#!/usr/bin/env python3
"""Complete archaeal PICRUSt2 placement from a validated HMM alignment."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.setrecursionlimit(200_000)

from picrust2.place_seqs import (  # noqa: E402
    check_alignments,
    gappa_jplace_to_newick,
    run_epa_ng,
)
from picrust2.util import read_fasta, read_stockholm, write_fasta  # noqa: E402

from pplacer_epa_ng_compat import read_fasta_ids, validate_jplace  # noqa: E402


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def newick_leaves(path: Path) -> list[str]:
    import dendropy

    tree = dendropy.Tree.get(
        path=str(path),
        schema="newick",
        preserve_underscores=True,
    )
    return [node.taxon.label for node in tree.leaf_node_iter()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-fasta", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--reference-tree", required=True, type=Path)
    parser.add_argument("--reference-model", required=True, type=Path)
    parser.add_argument("--query-stockholm", required=True, type=Path)
    parser.add_argument("--sepp-filtered-fasta", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--out-tree", required=True, type=Path)
    parser.add_argument("--audit-json", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--min-align", type=float, default=0.8)
    args = parser.parse_args()

    for path in [
        args.study_fasta,
        args.reference_fasta,
        args.reference_tree,
        args.reference_model,
        args.query_stockholm,
        args.sepp_filtered_fasta,
    ]:
        if not path.is_file():
            raise FileNotFoundError(path)

    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.out_tree.parent.mkdir(parents=True, exist_ok=True)
    args.audit_json.parent.mkdir(parents=True, exist_ok=True)

    aligned = read_stockholm(str(args.query_stockholm), clean_char=True)
    reference = read_fasta(str(args.reference_fasta))
    study = read_fasta(str(args.study_fasta))
    reference_ids = set(reference)
    reference_aligned = {sequence: aligned[sequence] for sequence in reference_ids}
    study_aligned = check_alignments(
        raw_seqs=study,
        aligned_seqs=aligned,
        min_align=args.min_align,
        verbose=True,
    )
    retained_ids = set(study_aligned)
    retained_from_sepp = set(read_fasta_ids(args.sepp_filtered_fasta))
    if retained_ids != retained_from_sepp:
        raise ValueError(
            "Recomputed and SEPP-stage archaeal alignment filters differ: "
            f"missing={len(retained_from_sepp - retained_ids)}, "
            f"extra={len(retained_ids - retained_from_sepp)}"
        )

    reference_msa = args.work_dir / "ref_seqs_hmmalign.fasta"
    study_msa = args.work_dir / "study_seqs_hmmalign.fasta"
    write_fasta(reference_aligned, str(reference_msa))
    write_fasta(study_aligned, str(study_msa))

    epa_dir = args.work_dir / "epa_out"
    run_epa_ng(
        tree=str(args.reference_tree),
        ref_msa_fastafile=str(reference_msa),
        study_msa_fastafile=str(study_msa),
        model=str(args.reference_model),
        out_dir=str(epa_dir),
        chunk_size=250,
        threads=args.threads,
        print_cmds=True,
    )
    parsed_jplace = epa_dir / "epa_result_parsed.jplace"
    validate_jplace(parsed_jplace, sorted(retained_ids))
    gappa_jplace_to_newick(
        jplace_file=str(parsed_jplace),
        outfile=str(args.out_tree),
        print_cmds=True,
    )

    output_leaves = newick_leaves(args.out_tree)
    expected_leaves = reference_ids | retained_ids
    if set(output_leaves) != expected_leaves or len(output_leaves) != len(
        expected_leaves
    ):
        raise ValueError(
            "Archaeal grafted-tree leaf audit failed: "
            f"missing={len(expected_leaves - set(output_leaves))}, "
            f"extra={len(set(output_leaves) - expected_leaves)}, "
            f"duplicates={len(output_leaves) - len(set(output_leaves))}"
        )

    failed_ids = sorted(set(study) - retained_ids)
    audit = {
        "status": "passed",
        "placement_implementation": (
            "PICRUSt2 whole-reference EPA-ng path using the existing "
            "domain-specific HMM alignment, standard 0.99 accumulated-LWR "
            "filter, maximum 100 placements, query chunk size 250, and gappa grafting."
        ),
        "input_ASVs": len(study),
        "ASVs_passing_archaeal_alignment_threshold": len(retained_ids),
        "ASVs_failing_archaeal_alignment_threshold": len(failed_ids),
        "failed_ASV_IDs": failed_ids,
        "reference_leaves": len(reference_ids),
        "output_tree_leaves": len(output_leaves),
        "unique_output_tree_leaves": len(set(output_leaves)),
        "missing_expected_leaves": 0,
        "unexpected_leaves": 0,
        "duplicate_leaf_names": 0,
        "min_align": args.min_align,
        "threads": args.threads,
        "study_fasta_sha256": sha256_file(args.study_fasta),
        "query_stockholm_sha256": sha256_file(args.query_stockholm),
        "parsed_jplace_sha256": sha256_file(parsed_jplace),
        "output_tree_sha256": sha256_file(args.out_tree),
    }
    args.audit_json.write_text(
        json.dumps(audit, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
