#!/usr/bin/env python3
"""
benchmark_orthologs.py — Ortholog benchmarking against a simulated reference.

Usage:
    python benchmark_orthologs.py --output-folder <sim_dir> --tool-output <tool_dir> --tool <tool>
                                  [--csv <results.csv>]

Arguments:
    --output-folder / -o    Simulation output folder (contains Simulated_orthologs.txt).
    --tool-output   / -to   The tool's output folder.
    --tool          / -t    Tool name: orthofinder, fastoma, broccoli, sonicparanoid2
    --csv                   Optional path to write results as a CSV file.

Tool-specific file locations (relative to --tool-output):
    orthofinder   : **/Orthologues/  (searched recursively)
    fastoma       : orthologs.tsv.gz
    broccoli      : dir_step4/orthologous_pairs.txt
    sonicparanoid2: sp_default/runs/*/species_to_species_orthologs/
"""

import argparse
import csv
import glob
import gzip
import os
import sys
from itertools import product
from pathlib import Path


SUPPORTED_TOOLS = ("orthofinder", "fastoma", "broccoli", "sonicparanoid2")


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def pair(g1, g2):
    return tuple(sorted((g1.strip(), g2.strip())))


def load_true_pairs(path):
    pairs = set()
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sp1 = row["Species_1"].strip()
            sp2 = row["Species_2"].strip()
            if sp1 == sp2:
                continue
            g1 = row["Gene_1"].strip()
            g2 = row["Gene_2"].strip()
            if g1 and g2 and g1 != g2:
                pairs.add(pair(g1, g2))
    return pairs


# ---------------------------------------------------------------------------
# OrthoFinder
# ---------------------------------------------------------------------------

def find_orthofinder_orthologues_dir(tool_output):
    matches = glob.glob(str(tool_output / "**" / "Orthologues"), recursive=True)
    if not matches:
        raise FileNotFoundError(
            f"Could not find an Orthologues folder under {tool_output}.\n"
            "Expected: <tool_output>/<any>/Results_*/Orthologues"
        )
    if len(matches) > 1:
        matches = sorted(matches, key=os.path.getmtime, reverse=True)
        print(f"Multiple Orthologues folders found; using most recent: {matches[0]}")
    return Path(matches[0])


def split_genes(x):
    if x is None or x.strip() == "":
        return []
    return [g.strip() for g in x.split(",") if g.strip()]


def load_orthofinder_pairs(folder):
    pairs = set()
    for tsv in Path(folder).rglob("*.tsv"):
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            cols = reader.fieldnames
            if not cols:
                continue
            focal_cols = [c for c in cols if c not in {"Orthogroup", "Species", "Orthologs"}]
            if len(focal_cols) != 1:
                continue
            focal_col = focal_cols[0]
            for row in reader:
                focal_genes    = split_genes(row[focal_col])
                ortholog_genes = split_genes(row["Orthologs"])
                for g1, g2 in product(focal_genes, ortholog_genes):
                    if g1 != g2:
                        pairs.add(pair(g1, g2))
    return pairs


# ---------------------------------------------------------------------------
# FastOMA
# ---------------------------------------------------------------------------

def load_fastoma_pairs(path):
    pairs = set()
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            g1, g2 = parts[0], parts[1]
            if g1 and g2 and g1 != g2:
                pairs.add(pair(g1, g2))
    return pairs


# ---------------------------------------------------------------------------
# Broccoli
# ---------------------------------------------------------------------------

def load_broccoli_pairs(path):
    pairs = set()
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            g1, g2 = parts[0], parts[1]
            if g1 and g2 and g1 != g2:
                pairs.add(pair(g1, g2))
    return pairs


# ---------------------------------------------------------------------------
# SonicParanoid2
# ---------------------------------------------------------------------------

def genes_from_sonic_field(x):
    """OrthoA/OrthoB fields: gene1 score gene2 score ... — genes at even indices."""
    tokens = x.strip().split()
    return tokens[0::2]


def find_sonic_ortholog_dir(tool_output):
    for base in [tool_output / "sp_default", tool_output]:
        runs_dir = base / "runs"
        if not runs_dir.is_dir():
            continue
        for run_dir in runs_dir.iterdir():
            if not run_dir.is_dir():
                continue
            ortholog_dir = run_dir / "species_to_species_orthologs"
            if ortholog_dir.is_dir():
                return ortholog_dir
    raise FileNotFoundError(
        f"No species_to_species_orthologs folder found under {tool_output}.\n"
        "Expected: <tool_output>/sp_default/runs/*/species_to_species_orthologs"
    )


def load_sonicparanoid_pairs(tool_output):
    pairs = set()
    ortholog_dir = find_sonic_ortholog_dir(tool_output)
    print(f"Ortholog dir:    {ortholog_dir}")
    for file in ortholog_dir.rglob("*"):
        if not file.is_file():
            continue
        with open(file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            if reader.fieldnames is None:
                continue
            if "OrthoA" not in reader.fieldnames or "OrthoB" not in reader.fieldnames:
                continue
            for row in reader:
                genes_a = genes_from_sonic_field(row["OrthoA"])
                genes_b = genes_from_sonic_field(row["OrthoB"])
                for g1, g2 in product(genes_a, genes_b):
                    if g1 and g2 and g1 != g2:
                        pairs.add(pair(g1, g2))
    return pairs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark ortholog predictions against a simulated reference.\n\n"
            "Requires two folders:\n"
            "  --output-folder  : simulation folder containing Simulated_orthologs.txt\n"
            "  --tool-output    : the tool's output folder"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--output-folder", "-o", required=True, dest="output_folder",
                        help="Simulation output folder (contains Simulated_orthologs.txt).")
    parser.add_argument("--tool-output", "-to", required=True, dest="tool_output",
                        help="Tool output folder.")
    parser.add_argument("--tool", "-t", required=True, dest="tool", choices=SUPPORTED_TOOLS,
                        help=f"Orthology tool. Supported: {', '.join(SUPPORTED_TOOLS)}")
    parser.add_argument("--csv", default=None, dest="csv",
                        help="Optional path to write results as a CSV file.")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    output_folder = Path(args.output_folder).resolve()
    tool_output   = Path(args.tool_output).resolve()
    tool          = args.tool.lower()

    for label, folder in [("output-folder", output_folder), ("tool-output", tool_output)]:
        if not folder.is_dir():
            print(f"ERROR: --{label} not found: {folder}")
            sys.exit(1)

    ref_file = output_folder / "Simulated_orthologs.txt"
    if not ref_file.exists():
        print(f"ERROR: reference file not found: {ref_file}")
        sys.exit(1)

    print(f"Tool:            {tool}")
    print(f"Simulation dir:  {output_folder}")
    print(f"Tool output dir: {tool_output}")

    true_pairs = load_true_pairs(ref_file)

    # --- Load inferred pairs ---
    if tool == "orthofinder":
        orthologues_dir = find_orthofinder_orthologues_dir(tool_output)
        print(f"Orthologues:     {orthologues_dir}")
        inferred_pairs = load_orthofinder_pairs(orthologues_dir)
    elif tool == "fastoma":
        ortholog_file = tool_output / "orthologs.tsv.gz"
        if not ortholog_file.exists():
            print(f"ERROR: FastOMA ortholog file not found: {ortholog_file}"); sys.exit(1)
        print(f"Ortholog file:   {ortholog_file}")
        inferred_pairs = load_fastoma_pairs(ortholog_file)
    elif tool == "broccoli":
        ortholog_file = tool_output / "dir_step4" / "orthologous_pairs.txt"
        if not ortholog_file.exists():
            print(f"ERROR: Broccoli ortholog file not found: {ortholog_file}"); sys.exit(1)
        print(f"Ortholog file:   {ortholog_file}")
        inferred_pairs = load_broccoli_pairs(ortholog_file)
    elif tool == "sonicparanoid2":
        inferred_pairs = load_sonicparanoid_pairs(tool_output)

    tp = len(true_pairs & inferred_pairs)
    fp = len(inferred_pairs - true_pairs)
    fn = len(true_pairs - inferred_pairs)

    precision = tp / (tp + fp) if tp + fp > 0 else 0
    recall    = tp / (tp + fn) if tp + fn > 0 else 0
    f1        = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0

    results = {
        "tool":              tool,
        "simulation_folder": str(output_folder),
        "tool_output":       str(tool_output),
        "TP":                tp,
        "FP":                fp,
        "FN":                fn,
        "Precision":         round(precision, 6),
        "Recall":            round(recall, 6),
        "F1":                round(f1, 6),
    }

    # --- Print results ---
    print("\nRESULTS\n")
    print(f"TP:        {tp}")
    print(f"FP:        {fp}")
    print(f"FN:        {fn}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"F1:        {f1:.4f}")

    # --- Write CSV ---
    if args.csv:
        csv_path = Path(args.csv)
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with open(csv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(results.keys()))
            writer.writeheader()
            writer.writerow(results)
        print(f"\nResults written to: {csv_path}")


if __name__ == "__main__":
    main()
