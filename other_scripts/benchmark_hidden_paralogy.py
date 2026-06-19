#!/usr/bin/env python3
"""
benchmark_hidden_paralogy.py — Hidden paralogy error rate benchmarking.

What is hidden paralogy?
    A gene duplication creates copies 1 and 2 in ancestor of species A and B.
    Species A loses copy 2; species B loses copy 1. A_copy1 and B_copy2 survive
    in the same orthogroup but their LCA in the gene tree was a Duplication, not
    a Speciation. A tool that calls them orthologs is making a hidden paralogy error.

Definitions:
    False-friend pair:
        species(gene1) != species(gene2)
        AND both genes are in the same simulated orthogroup
        AND LCA(gene1, gene2) in the gene tree has VERTEXTYPE=Duplication

    Hidden paralogy error rate:
        numerator   = false-friend pairs called as orthologs by the tool
        denominator = false-friend pairs where both genes are in the same
                      inferred orthogroup (the tool had the opportunity to classify them)

Usage:
    python benchmark_hidden_paralogy.py
        --sim-folder  <simulation_output_dir>
        --tool-output <tool_output_dir>
        --tool        orthofinder|fastoma|broccoli|sonicparanoid2
        [--csv        results.csv]

Tool-specific orthogroup file locations (relative to --tool-output):
    orthofinder   : **/Orthogroups/Orthogroups.txt  (searched recursively)
    fastoma       : RootHOGs.tsv
    broccoli      : dir_step3/orthologous_groups.txt
    sonicparanoid2: sp_default/runs/*/ortholog_groups/flat.ortholog_groups.tsv

Tool-specific ortholog file locations (relative to --tool-output):
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
import re
import sys
from collections import defaultdict
from itertools import product
from pathlib import Path


SUPPORTED_TOOLS = ("orthofinder", "fastoma", "broccoli", "sonicparanoid2")


# =============================================================================
# Gene name helper
# =============================================================================

def species_from_gene(gene):
    """
    Extract species name from a simulated gene name.
    e.g. Homo_sapiens_5_1 → Homo_sapiens  (strips last two tokens)
    """
    parts = gene.split("_")
    return "_".join(parts[:-2])


# =============================================================================
# Gene-tree parsing — find inter-species false-friend pairs
# =============================================================================

def parse_vtype_map(tree_str):
    """
    Extract {XTAGX_name: VERTEXTYPE} from the annotated Newick string produced
    by PrunedTree output.  Annotations have the form:
        XTAGX6_15:0.04[ID=6 HOST=15 VERTEXTYPE=Speciation DISCPT=(1,0)]
    """
    vtype = {}
    for m in re.finditer(r"(XTAGX\w+):[0-9e.+\-]+\[([^\]]+)\]", tree_str):
        vt = re.search(r"VERTEXTYPE=(\w+)", m.group(2))
        if vt:
            vtype[m.group(1)] = vt.group(1)
    return vtype


def strip_annotations(tree_str):
    """Remove all [...] annotation blocks so the result is plain Newick."""
    return re.sub(r"\[[^\]]*\]", "", tree_str).strip().rstrip(";")


def _split_children(s):
    """Split a Newick child string on commas at depth 0."""
    parts, depth, start = [], 0, 0
    for i, c in enumerate(s):
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        elif c == "," and depth == 0:
            parts.append(s[start:i])
            start = i + 1
    parts.append(s[start:])
    return parts


def parse_newick(s):
    """
    Parse a plain Newick string into a lightweight tuple tree:
        (node_name, [child_tuples], frozenset_of_leaf_names)
    Leaf names are the XTAGX identifiers.
    """
    s = s.strip()
    if s.startswith("("):
        depth, end = 0, -1
        for i, c in enumerate(s):
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    end = i
                    break
        children = [parse_newick(k) for k in _split_children(s[1:end])]
        m = re.match(r"(XTAGX\w+)", s[end + 1:])
        name = m.group(1) if m else None
        leaves = frozenset(lf for _, _, ls in children for lf in ls)
        return (name, children, leaves)
    else:
        m = re.match(r"(XTAGX\w+)", s)
        name = m.group(1) if m else None
        return (name, [], frozenset([name]) if name else frozenset())


def _collect_false_friends(node, vtype, rmap, pairs):
    """
    Recursively walk the tree.  At each Duplication node, add all inter-species
    pairs drawn from the cross-product of each pair of child subtrees.
    """
    name, children, _ = node
    if not children:
        return

    if vtype.get(name) == "Duplication":
        for i in range(len(children)):
            for j in range(i + 1, len(children)):
                _, _, li = children[i]
                _, _, lj = children[j]
                for xn1 in li:
                    for xn2 in lj:
                        g1 = rmap.get(xn1)
                        g2 = rmap.get(xn2)
                        if g1 and g2 and species_from_gene(g1) != species_from_gene(g2):
                            pairs.add(frozenset([g1, g2]))

    for child in children:
        _collect_false_friends(child, vtype, rmap, pairs)


def _load_relabelmap(path):
    """Load XTAGX_name → final_gene_name from relabelmap TSV."""
    rmap = {}
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                rmap[parts[0]] = parts[1]
    return rmap


def _ogs_with_duplications(sim_dir):
    """Return set of integer OG indices that have num_duplications > 0."""
    stats_file = Path(sim_dir) / "Simulated_Orthogroup_statistics.txt"
    with_dups = set()
    with open(stats_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if int(row.get("num_duplications", 0)) > 0:
                with_dups.add(int(row["orthogroup"]))
    return with_dups


def get_all_false_friends(sim_dir):
    """
    Return a set of frozensets, each containing two gene names (from different
    species) whose LCA in the simulated gene tree is a Duplication node.
    """
    sim_dir = Path(sim_dir)
    temp_dir = sim_dir / "temporary_files"
    all_pairs = set()

    for og_idx in _ogs_with_duplications(sim_dir):
        tree_file = temp_dir / f"{og_idx}.pruned.tree"
        rmap_file = temp_dir / f"{og_idx}.relabelmap.txt"
        if not tree_file.exists() or not rmap_file.exists():
            continue

        tree_str = tree_file.read_text()
        vtype    = parse_vtype_map(tree_str)
        tree     = parse_newick(strip_annotations(tree_str))
        rmap     = _load_relabelmap(rmap_file)

        _collect_false_friends(tree, vtype, rmap, all_pairs)

    return all_pairs


# =============================================================================
# Tool orthogroup membership: gene → og_id
# =============================================================================

def _gene_to_og_orthofinder(tool_output):
    matches = glob.glob(
        str(tool_output / "**" / "Orthogroups" / "Orthogroups.txt"), recursive=True
    )
    if not matches:
        raise FileNotFoundError(
            f"Orthogroups.txt not found under {tool_output}.\n"
            "Expected: <tool_output>/**/Orthogroups/Orthogroups.txt"
        )
    og_file = sorted(matches, key=os.path.getmtime, reverse=True)[0]
    print(f"OG file:         {og_file}")
    g2og = {}
    with open(og_file) as f:
        for line in f:
            line = line.strip()
            if not line or ":" not in line:
                continue
            og_id, rest = line.split(":", 1)
            og_id = og_id.strip()
            for gene in rest.split():
                gene = gene.strip().rstrip(",")
                if gene:
                    g2og[gene] = og_id
    return g2og


def _gene_to_og_fastoma(tool_output):
    og_file = tool_output / "RootHOGs.tsv"
    if not og_file.exists():
        raise FileNotFoundError(f"FastOMA RootHOGs.tsv not found: {og_file}")
    print(f"OG file:         {og_file}")
    g2og = {}
    with open(og_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            og_id = parts[0].rstrip(":").strip()
            gene  = parts[1].strip()
            if gene and og_id.lower() not in {"hog", "roothog", "orthogroup"}:
                g2og[gene] = og_id
    return g2og


def _gene_to_og_broccoli(tool_output):
    og_file = tool_output / "dir_step3" / "orthologous_groups.txt"
    if not og_file.exists():
        raise FileNotFoundError(f"Broccoli OG file not found: {og_file}")
    print(f"OG file:         {og_file}")
    g2og = {}
    with open(og_file) as f:
        next(f, None)  # skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            og_id = parts[0]
            for gene in parts[1:]:
                if gene and gene != "*":
                    g2og[gene] = og_id
    return g2og


def _gene_to_og_sonicparanoid2(tool_output):
    for base in [tool_output / "sp_default", tool_output]:
        pattern = str(base / "runs" / "*" / "ortholog_groups" / "flat.ortholog_groups.tsv")
        matches = glob.glob(pattern)
        if matches:
            break
    if not matches:
        raise FileNotFoundError(
            f"flat.ortholog_groups.tsv not found under {tool_output}."
        )
    og_file = sorted(matches, key=os.path.getmtime, reverse=True)[0]
    print(f"OG file:         {og_file}")
    g2og = {}
    with open(og_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        if "group_id" not in header:
            raise ValueError(f"'group_id' column not found in {og_file}")
        group_idx = header.index("group_id")
        sp_indices = [i for i in range(len(header)) if i != group_idx]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= group_idx:
                continue
            og_id = parts[group_idx].strip()
            if not og_id:
                continue
            for idx in sp_indices:
                if idx >= len(parts):
                    continue
                for gene in parts[idx].split(","):
                    gene = gene.strip()
                    if gene and gene not in {"", "*", "NA"}:
                        g2og[gene] = og_id
    return g2og


_GENE_TO_OG = {
    "orthofinder":   _gene_to_og_orthofinder,
    "fastoma":       _gene_to_og_fastoma,
    "broccoli":      _gene_to_og_broccoli,
    "sonicparanoid2": _gene_to_og_sonicparanoid2,
}


def load_gene_to_og(tool, tool_output):
    return _GENE_TO_OG[tool](tool_output)


# =============================================================================
# Tool ortholog pairs (gene, gene) called explicitly as orthologs
# =============================================================================

def _pair(g1, g2):
    return frozenset([g1.strip(), g2.strip()])


def _split_genes(x):
    if not x or not x.strip():
        return []
    return [g.strip() for g in x.split(",") if g.strip()]


def _load_ortholog_pairs_orthofinder(tool_output):
    matches = glob.glob(str(tool_output / "**" / "Orthologues"), recursive=True)
    if not matches:
        raise FileNotFoundError(f"No Orthologues folder found under {tool_output}.")
    orthologues_dir = sorted(matches, key=os.path.getmtime, reverse=True)[0]
    print(f"Orthologues:     {orthologues_dir}")
    pairs = set()
    for tsv in Path(orthologues_dir).rglob("*.tsv"):
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            cols = reader.fieldnames or []
            focal = [c for c in cols if c not in {"Orthogroup", "Species", "Orthologs"}]
            if len(focal) != 1:
                continue
            focal_col = focal[0]
            for row in reader:
                for g1, g2 in product(_split_genes(row[focal_col]), _split_genes(row["Orthologs"])):
                    if g1 != g2:
                        pairs.add(_pair(g1, g2))
    return pairs


def _load_ortholog_pairs_fastoma(tool_output):
    path = tool_output / "orthologs.tsv.gz"
    if not path.exists():
        raise FileNotFoundError(f"FastOMA ortholog file not found: {path}")
    print(f"Ortholog file:   {path}")
    pairs = set()
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] != parts[1]:
                pairs.add(_pair(parts[0], parts[1]))
    return pairs


def _load_ortholog_pairs_broccoli(tool_output):
    path = tool_output / "dir_step4" / "orthologous_pairs.txt"
    if not path.exists():
        raise FileNotFoundError(f"Broccoli ortholog file not found: {path}")
    print(f"Ortholog file:   {path}")
    pairs = set()
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] != parts[1]:
                pairs.add(_pair(parts[0], parts[1]))
    return pairs


def _load_ortholog_pairs_sonicparanoid2(tool_output):
    for base in [tool_output / "sp_default", tool_output]:
        runs_dir = base / "runs"
        if not runs_dir.is_dir():
            continue
        for run_dir in runs_dir.iterdir():
            ortholog_dir = run_dir / "species_to_species_orthologs"
            if ortholog_dir.is_dir():
                print(f"Ortholog dir:    {ortholog_dir}")
                pairs = set()
                for f_path in ortholog_dir.rglob("*"):
                    if not f_path.is_file():
                        continue
                    with open(f_path) as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        if not reader.fieldnames:
                            continue
                        if "OrthoA" not in reader.fieldnames:
                            continue
                        for row in reader:
                            ga = row["OrthoA"].split()[0::2]
                            gb = row["OrthoB"].split()[0::2]
                            for g1, g2 in product(ga, gb):
                                if g1 and g2 and g1 != g2:
                                    pairs.add(_pair(g1, g2))
                return pairs
    raise FileNotFoundError(
        f"No species_to_species_orthologs folder found under {tool_output}."
    )


_LOAD_ORTHOLOG_PAIRS = {
    "orthofinder":   _load_ortholog_pairs_orthofinder,
    "fastoma":       _load_ortholog_pairs_fastoma,
    "broccoli":      _load_ortholog_pairs_broccoli,
    "sonicparanoid2": _load_ortholog_pairs_sonicparanoid2,
}


def load_ortholog_pairs(tool, tool_output):
    return _LOAD_ORTHOLOG_PAIRS[tool](tool_output)


# =============================================================================
# Metric computation
# =============================================================================

def compute_hidden_paralogy_rate(false_friends, gene_to_og, ortholog_pairs):
    """
    Returns:
        ff_total        : total false-friend pairs in the simulation
        ff_same_og      : false-friend pairs where both genes are in the same tool OG
        ff_called_ortho : false-friend pairs explicitly called as orthologs
        error_rate      : ff_called_ortho / ff_same_og  (None if ff_same_og == 0)
    """
    ff_total        = len(false_friends)
    ff_same_og      = 0
    ff_called_ortho = 0

    for pair in false_friends:
        g1, g2 = tuple(pair)
        og1 = gene_to_og.get(g1)
        og2 = gene_to_og.get(g2)
        in_same_og = (og1 is not None) and (og1 == og2)

        if in_same_og:
            ff_same_og += 1
            if pair in ortholog_pairs:
                ff_called_ortho += 1

    error_rate = ff_called_ortho / ff_same_og if ff_same_og > 0 else None
    return ff_total, ff_same_og, ff_called_ortho, error_rate


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark hidden paralogy error rate against a simulated reference.\n\n"
            "A false-friend pair is two genes from different species in the same\n"
            "simulated orthogroup whose LCA in the gene tree is a Duplication node.\n\n"
            "Error rate = (false friends called as orthologs) / (false friends in same tool OG)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--sim-folder", "-s", required=True, dest="sim_folder",
        help="Simulation output folder (contains temporary_files/ and Simulated_Orthogroup_statistics.txt).",
    )
    parser.add_argument(
        "--tool-output", "-to", required=True, dest="tool_output",
        help="Tool output folder.",
    )
    parser.add_argument(
        "--tool", "-t", required=True, dest="tool", choices=SUPPORTED_TOOLS,
        help=f"Orthology tool. Supported: {', '.join(SUPPORTED_TOOLS)}",
    )
    parser.add_argument(
        "--csv", default=None, dest="csv",
        help="Optional path to write results as a CSV file.",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    sim_folder  = Path(args.sim_folder).resolve()
    tool_output = Path(args.tool_output).resolve()
    tool        = args.tool.lower()

    for label, path in [("sim-folder", sim_folder), ("tool-output", tool_output)]:
        if not path.is_dir():
            print(f"ERROR: --{label} not found: {path}")
            sys.exit(1)

    stats_file = sim_folder / "Simulated_Orthogroup_statistics.txt"
    if not stats_file.exists():
        print(f"ERROR: {stats_file} not found — is this a valid simulation folder?")
        sys.exit(1)

    print(f"Tool:            {tool}")
    print(f"Simulation dir:  {sim_folder}")
    print(f"Tool output dir: {tool_output}")

    # --- 1. Extract false-friend pairs ---
    print("\nScanning gene trees for false-friend pairs...")
    false_friends = get_all_false_friends(sim_folder)
    print(f"False-friend pairs found: {len(false_friends)}")

    if not false_friends:
        print(
            "\nNo cross-species false-friend pairs found in this simulation.\n"
            "Hidden paralogy error rate is undefined (no false friends to evaluate)."
        )
        if args.csv:
            csv_path = Path(args.csv)
            csv_path.parent.mkdir(parents=True, exist_ok=True)
            with open(csv_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["tool", "sim_folder", "tool_output",
                                                    "ff_total", "ff_same_og",
                                                    "ff_called_orthologs", "error_rate"])
                w.writeheader()
                w.writerow({"tool": tool, "sim_folder": str(sim_folder),
                             "tool_output": str(tool_output),
                             "ff_total": 0, "ff_same_og": 0,
                             "ff_called_orthologs": 0, "error_rate": "NA"})
            print(f"Results written to: {csv_path}")
        return

    # --- 2. Load tool orthogroup membership ---
    print("\nLoading tool orthogroup membership...")
    gene_to_og = load_gene_to_og(tool, tool_output)
    print(f"Genes in tool OGs: {len(gene_to_og)}")

    # --- 3. Load tool ortholog pairs ---
    print("\nLoading tool ortholog pairs...")
    ortholog_pairs = load_ortholog_pairs(tool, tool_output)
    print(f"Ortholog pairs:  {len(ortholog_pairs)}")

    # --- 4. Compute metric ---
    ff_total, ff_same_og, ff_called_ortho, error_rate = compute_hidden_paralogy_rate(
        false_friends, gene_to_og, ortholog_pairs
    )

    error_rate_str = f"{error_rate:.6f}" if error_rate is not None else "NA (no false friends in same tool OG)"

    print("\nRESULTS\n")
    print(f"False-friend pairs (simulation):             {ff_total}")
    print(f"False-friend pairs in same tool OG:          {ff_same_og}   [denominator]")
    print(f"False-friend pairs called as orthologs:      {ff_called_ortho}   [numerator]")
    print(f"Hidden paralogy error rate:                  {error_rate_str}")

    results = {
        "tool":                  tool,
        "sim_folder":            str(sim_folder),
        "tool_output":           str(tool_output),
        "ff_total":              ff_total,
        "ff_same_og":            ff_same_og,
        "ff_called_orthologs":   ff_called_ortho,
        "error_rate":            error_rate if error_rate is not None else "NA",
    }

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
