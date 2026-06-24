#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import statistics
from pathlib import Path
from multiprocessing import Pool
from ete3 import Tree
from itertools import combinations

# --- USER SETTINGS ---
MAX_WORKERS = 8  # number of parallel processes
MIN_GENES = 3    # only summarize orthogroups with at least this many genes/sequences


# ============================================================
# Tree processing
# ============================================================

def process_tree(fname, tree_dir):
    """
    Process one .tre file.

    Only trees with at least MIN_GENES leaves are included.

    Calculates:
        - median root-to-tip distance
        - Wiener index: sum of all pairwise leaf distances
        - treeness: internal branch length / total branch length
    """
    if not fname.endswith(".tre"):
        return None

    path = os.path.join(tree_dir, fname)

    try:
        tree = Tree(path, format=1)

        root = tree.get_tree_root()
        leaves = list(tree.iter_leaves())

        # only keep trees with at least 4 genes/leaves
        if len(leaves) < MIN_GENES:
            return None

        rtt = [root.get_distance(leaf) for leaf in leaves]
        median_rtt = statistics.median(rtt) if rtt else None

        wiener = sum(
            tree.get_distance(a, b)
            for a, b in combinations(leaves, 2)
        )

        branch_lengths = [
            n.dist
            for n in tree.traverse()
            if n.up is not None and n.dist is not None and n.dist > 0
        ]

        total_bl = sum(branch_lengths)

        internal_bl = sum(
            n.dist
            for n in tree.traverse()
            if (
                n.up is not None
                and not n.is_leaf()
                and n.dist is not None
                and n.dist > 0
            )
        )

        treeness = internal_bl / total_bl if total_bl > 0 else None
        if treeness >0.999:
            print(f"high treeness tree: {fname}")
        return median_rtt, wiener, treeness

    except Exception as e:
        print(f"Tree error: {fname} ({e})")

    return None


# ============================================================
# Gap processing
# ============================================================

def process_msa_gaps(msa_file):
    """
    Read one FASTA alignment file.

    Only alignments with at least MIN_GENES sequences are included.

    For each sequence, calculate:
        - number of gap regions
        - median gap-region length

    Only sequences with at least one gap are included in the gap outputs.
    """
    sequences = []

    with open(msa_file) as f:
        name = None
        seq = []

        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    sequences.append("".join(seq))

                name = line[1:]
                seq = []

            else:
                seq.append(line)

        if name is not None:
            sequences.append("".join(seq))

    # only keep alignments with at least 4 genes/sequences
    if len(sequences) < MIN_GENES:
        return [], []

    num_gaps_list = []
    median_gap_list = []

    for seq in sequences:
        gap_lengths = []
        current_gap = 0

        for char in seq:
            if char == "-":
                current_gap += 1
            else:
                if current_gap > 0:
                    gap_lengths.append(current_gap)
                    current_gap = 0

        if current_gap > 0:
            gap_lengths.append(current_gap)

        if gap_lengths:
            num_gaps_list.append(len(gap_lengths))
            median_gap_list.append(statistics.median(gap_lengths))

    return num_gaps_list, median_gap_list


def process_alignment_gap_stats(alignment_dir):
    """
    Process all alignments in alignment_files.

    Only alignments with at least MIN_GENES sequences are included.

    Returns:
        all_num_gaps, all_median_gaps
    """
    all_num_gaps = []
    all_median_gaps = []

    alignment_dir = Path(alignment_dir)

    if not alignment_dir.is_dir():
        print(f"Warning: alignment folder not found: {alignment_dir}")
        return all_num_gaps, all_median_gaps

    for file in alignment_dir.rglob("*.alignment.fa"):
        n, m = process_msa_gaps(file)
        all_num_gaps.extend(n)
        all_median_gaps.extend(m)

    return all_num_gaps, all_median_gaps


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", required=True, help="Simulation folder")
    args = parser.parse_args()

    sim_folder = args.folder

    stats_file = os.path.join(
        sim_folder,
        "Simulated_Orthogroup_statistics.txt"
    )

    tree_dir = os.path.join(
        sim_folder,
        "tree_files"
    )

    alignment_dir = os.path.join(
        sim_folder,
        "alignment_files"
    )

    output_dir = os.path.join(
        sim_folder,
        "simulation_summaries"
    )

    os.makedirs(output_dir, exist_ok=True)

    # ========================================================
    # 1. Parse table
    #    Only keep orthogroups with at least MIN_GENES genes
    # ========================================================

    num_duplications = []
    num_species = []
    num_genes = []

    with open(stats_file) as f:
        header = f.readline().strip().split("\t")

        required_cols = [
            "num_duplications",
            "num_species",
            "num_genes"
        ]

        for col in required_cols:
            if col not in header:
                raise ValueError(f"Missing column: {col}")

        dup_idx = header.index("num_duplications")
        species_idx = header.index("num_species")
        genes_idx = header.index("num_genes")

        total_rows = 0
        kept_rows = 0

        for line in f:
            parts = line.strip().split("\t")

            if len(parts) < len(header):
                continue

            total_rows += 1

            try:
                this_num_genes = int(parts[genes_idx])

                # only keep orthogroups with at least 3 genes
                if this_num_genes < MIN_GENES:
                    continue

                num_duplications.append(int(parts[dup_idx]))
                num_species.append(int(parts[species_idx]))
                num_genes.append(this_num_genes)

                kept_rows += 1

            except ValueError:
                continue

    print(
        f"Kept {kept_rows} / {total_rows} orthogroups "
        f"with at least {MIN_GENES} genes from statistics table."
    )

    # ========================================================
    # 2. Multiprocess tree processing
    #    Only trees with at least MIN_GENES leaves are retained
    # ========================================================

    tree_files = [
        f for f in os.listdir(tree_dir)
        if f.endswith(".tre")
    ]

    with Pool(processes=MAX_WORKERS) as pool:
        results = pool.starmap(
            process_tree,
            [(fname, tree_dir) for fname in tree_files]
        )

    median_rtt_values = []
    wiener_values = []
    treeness_values = []

    for r in results:
        if r is None:
            continue

        median_rtt, wiener, treeness = r

        if median_rtt is not None:
            median_rtt_values.append(median_rtt)

        if wiener is not None:
            wiener_values.append(wiener)

        if treeness is not None:
            treeness_values.append(treeness)

    print(
        f"Processed {len(median_rtt_values)} trees "
        f"with at least {MIN_GENES} leaves."
    )

    if not median_rtt_values:
        print("Warning: no valid trees processed.")

    # ========================================================
    # 3. Process alignment gaps
    #    Only alignments with at least MIN_GENES sequences kept
    # ========================================================

    num_gaps, median_gap_size = process_alignment_gap_stats(
        alignment_dir
    )

    if not num_gaps:
        print("Warning: no gap statistics collected.")

    # ========================================================
    # 4. Write output files
    # ========================================================

    def write_list(filename, data):
        with open(os.path.join(output_dir, filename), "w") as f:
            for x in data:
                f.write(f"{x}\n")

    write_list("duplication_counts.txt", num_duplications)
    write_list("num_species.txt", num_species)
    write_list("num_genes.txt", num_genes)

    write_list("median_rtt.txt", median_rtt_values)
    write_list("wiener_index.txt", wiener_values)
    write_list("treeness.txt", treeness_values)

    write_list("num_gaps.txt", num_gaps)
    write_list("median_gap_size.txt", median_gap_size)

    print(f"Wrote summaries to: {output_dir}")