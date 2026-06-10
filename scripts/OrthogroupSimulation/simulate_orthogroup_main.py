#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Orthogroup simulation entrypoint for Orthosim.py.

Important design:
  - Orthosim.py is the master wrapper.
  - Orthosim.py requires a config for --Orthogroup-simulation.
  - Orthosim.py may also receive extra CLI flags; those override the input config.
  - Orthosim.py writes the final merged config to <output>/config.txt.
  - This script reads that final merged config and runs the simulation.

This script does not try to rediscover a generic input folder.
It adapts the Orthosim.py config values into the legacy args expected by
orthogroup_simulation_utils.py.
"""

import argparse
import os
import sys
import time
import multiprocessing as mp
from pathlib import Path

import orthogroup_simulation_utils as og


def read_config(config_path):
    """
    Read key=value config file.

    Blank lines and comment lines are ignored.
    """
    config = {}

    with open(config_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            if "=" not in line:
                raise ValueError(
                    f"Bad config line {line_number} in {config_path}: {line}"
                )

            key, value = line.split("=", 1)
            config[key.strip()] = value.strip()

    return config


def require_config(config, key):
    """
    Return a required config value or fail clearly.
    """
    value = config.get(key)

    if value is None or str(value).strip() == "":
        raise ValueError(f"Required value missing from Orthosim config: {key}")

    return value


def resolve_path(path_value, project_root):
    """
    Resolve paths from config.

    Absolute paths are used as-is.
    Relative paths are interpreted relative to the Orthosim.py project root.
    """
    path = Path(path_value)

    if path.is_absolute():
        return str(path)

    return str((Path(project_root) / path).resolve())


def first_existing_file(candidates, description):
    """
    Return first existing file from a list of candidate paths.
    """
    for candidate in candidates:
        if candidate and os.path.isfile(candidate):
            return candidate

    tried = "\n".join(f"  {candidate}" for candidate in candidates if candidate)
    raise ValueError(f"Could not find {description}. Tried:\n{tried}")


def infer_project_root(config):
    """
    Orthosim.py writes SimPath=<directory containing Orthosim.py>.
    That is treated as the project root.

    Fallback is two directories above this script:
      scripts/OrthogroupSimulation/../..
    """
    if config.get("SimPath"):
        return os.path.abspath(config["SimPath"])

    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(script_dir, "..", ".."))


def write_runtime_parameter_file(config, output_dir, project_root, threads):
    """
    orthogroup_simulation_utils.LoadParameters() expects a parameter file
    containing key=value lines.

    The Orthosim output config is already the final merged config, so we copy
    it into a runtime parameter file and add derived internal defaults.

    These derived values do not need to be user-facing config entries:
      iqtree_path  = <project_root>/bin/iqtree3
      sagephy_path = <project_root>/bin/sagephy-1.0.0.jar
    """
    runtime_config = dict(config)

    runtime_config["threads"] = str(threads)

    runtime_config.setdefault(
        "iqtree_path",
        os.path.join(project_root, "bin", "iqtree3"),
    )

    runtime_config.setdefault(
        "sagephy_path",
        os.path.join(project_root, "bin", "sagephy-1.0.0.jar"),
    )

    # The current utils use sigma_log_mean and sigma_log_sd.
    # Keep compatibility with configs that only have max_sigma2.
    runtime_config.setdefault("sigma_log_mean", "0")

    if "sigma_log_sd" not in runtime_config:
        if "max_sigma2" in runtime_config:
            runtime_config["sigma_log_sd"] = runtime_config["max_sigma2"]
        else:
            runtime_config["sigma_log_sd"] = "1"

    runtime_parameter_file = os.path.join(
        output_dir,
        "orthogroup_simulation_parameters.txt",
    )

    with open(runtime_parameter_file, "w") as out:
        for key in sorted(runtime_config):
            out.write(f"{key}={runtime_config[key]}\n")

    return runtime_parameter_file


def configure_legacy_args(args, config):
    """
    Convert Orthosim.py config names into the args names expected by
    orthogroup_simulation_utils.py.

    The utils currently expect:
      args.o
      args.s
      args.p
      args.f
      args.PFAM
      args.d
      args.gap_profile
    """
    project_root = infer_project_root(config)

    args.o = os.path.abspath(args.output)
    args.n = int(require_config(config, "Orthogroups"))
    args.s = resolve_path(require_config(config, "ultrametric_tree"), project_root)
    args.gap_profile = resolve_path(require_config(config, "gap_file"), project_root)

    # pfam_species allows domain species to differ from the species-tree label.
    # If absent, use species.
    domain_species = require_config(config, "species")
    args.domain_species = domain_species

    args.p = write_runtime_parameter_file(
        config=config,
        output_dir=args.o,
        project_root=project_root,
        threads=args.threads,
    )

    pfam_root = config.get("pfam_dir")
    if pfam_root:
        pfam_root = resolve_path(pfam_root, project_root)
    else:
        pfam_root = os.path.join(project_root, "PFAM")

    # Your stated PFAM layout:
    #   PFAM/pfam_genomes
    #   PFAM/pfam_results/csv_files/species_models_fixed
    #
    # Optional explicit overrides are supported, but not required.
    genome_candidates = []

    if config.get("pfam_genome_file"):
        genome_candidates.append(resolve_path(config["pfam_genome_file"], project_root))

    genome_candidates.extend([
        os.path.join(pfam_root, "pfam_genomes", f"{domain_species}.fa"),
        os.path.join(pfam_root, "pfam_genomes", f"{domain_species}.fasta"),
        os.path.join(pfam_root, "pfam_genomes", f"{domain_species}.faa"),
    ])

    args.f = first_existing_file(
        genome_candidates,
        f"PFAM genome FASTA for species '{domain_species}'",
    )

    domain_csv_candidates = []

    if config.get("pfam_domain_csv"):
        domain_csv_candidates.append(resolve_path(config["pfam_domain_csv"], project_root))

    domain_csv_candidates.extend([
        os.path.join(
            pfam_root,
            "pfam_results",
            "csv_files",
            "species_models_fixed",
            f"{domain_species}_pfam_models.csv",
        ),
        os.path.join(
            pfam_root,
            "pfam_results",
            "csv_files",
            "species_models_fixed",
            f"{domain_species}.csv",
        ),
    ])

    args.d = first_existing_file(
        domain_csv_candidates,
        f"PFAM domain-model CSV for species '{domain_species}'",
    )

    # orthogroup_simulation_utils.ExtractPfamDomains() needs args.PFAM.
    # This is the PFAM scan result text file.
    # We support explicit config override plus common locations.
    pfam_scan_candidates = []

    if config.get("pfam_scan_file"):
        pfam_scan_candidates.append(resolve_path(config["pfam_scan_file"], project_root))

    pfam_scan_candidates.extend([
        os.path.join(pfam_root, "pfam_results", f"{domain_species}_pfam.txt"),
        os.path.join(pfam_root, "pfam_results", f"{domain_species}.pfam.txt"),
        os.path.join(pfam_root, "pfam_results", f"{domain_species}.txt"),
        os.path.join(pfam_root, "pfam_results", "txt_files", f"{domain_species}_pfam.txt"),
        os.path.join(pfam_root, "pfamscan_results", f"{domain_species}_pfam.txt"),
    ])

    args.PFAM = first_existing_file(
        pfam_scan_candidates,
        f"PFAM scan results text file for species '{domain_species}'",
    )

    args.project_root = project_root
    args.pfam_root = pfam_root

    return args


def create_subfolders(output_dir):
    """
    Orthosim.py creates the top-level output directory.
    The simulation script creates only its own subdirectories.
    """
    os.makedirs(output_dir, exist_ok=True)

    for folder in [
        "temporary_files",
        "proteome_files",
        "alignment_files",
        "tree_files",
    ]:
        os.makedirs(os.path.join(output_dir, folder), exist_ok=True)


def check_inputs_and_load_parameters(args):
    """
    Fail early before multiprocessing starts.
    """
    required_files = {
        "Orthosim config": args.config,
        "species tree": args.s,
        "gap profile": args.gap_profile,
        "runtime parameter file": args.p,
        "PFAM genome FASTA": args.f,
        "PFAM scan results": args.PFAM,
        "PFAM domain-model CSV": args.d,
    }

    for description, path in required_files.items():
        if not os.path.isfile(path):
            raise ValueError(f"{description} not found: {path}")

    # Load runtime parameters into orthogroup_simulation_utils globals.
    og.args = args
    og.LoadParameters(args.p)

    if not os.path.isfile(str(og.sagephy_path)):
        raise ValueError(f"SagePhy jar not found: {og.sagephy_path}")

    if not os.path.isfile(str(og.iqtree_path)):
        raise ValueError(f"IQ-TREE executable not found: {og.iqtree_path}")

    if not os.access(str(og.iqtree_path), os.X_OK):
        raise ValueError(f"IQ-TREE exists but is not executable: {og.iqtree_path}")


def print_detected_inputs(args):
    print("Detected Orthosim-controlled inputs:")
    print(f"  Config:             {args.config}")
    print(f"  Project root:       {args.project_root}")
    print(f"  Output directory:   {args.o}")
    print(f"  Runtime parameters: {args.p}")
    print(f"  Orthogroups:        {args.n}")
    print(f"  Threads:            {args.threads}")
    print(f"  Species tree:       {args.s}")
    print(f"  Gap profile:        {args.gap_profile}")
    print(f"  PFAM root:          {args.pfam_root}")
    print(f"  PFAM species:       {args.domain_species}")
    print(f"  Genome FASTA:       {args.f}")
    print(f"  PFAM scan results:  {args.PFAM}")
    print(f"  Domain CSV:         {args.d}")
    print(f"  IQ-TREE:            {og.iqtree_path}")
    print(f"  SagePhy:            {og.sagephy_path}")


def main():
    print("---------------------------------------------------")
    print("Welcome to the Orthogroup Simulation\n")

    start_time = time.time()
    current_hour = time.localtime().tm_hour

    if 6 <= current_hour < 12:
        print("Good morning!")
    elif 12 <= current_hour < 18:
        print("Good afternoon!")
    else:
        print("Good evening!")

    parser = argparse.ArgumentParser(
        description=(
            "Orthogroup simulation worker. "
            "This script is called by Orthosim.py and reads the final "
            "Orthosim-generated output/config.txt."
        )
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Path to the final Orthosim-generated config.txt",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output directory created by Orthosim.py",
    )

    parser.add_argument(
        "--threads",
        required=True,
        type=int,
        help="Thread/process count passed by Orthosim.py",
    )

    args = parser.parse_args()

    args.config = os.path.abspath(args.config)
    args.output = os.path.abspath(args.output)

    if args.threads < 1:
        raise ValueError("--threads must be >= 1")

    if not os.path.isfile(args.config):
        raise ValueError(f"Orthosim config not found: {args.config}")

    if not os.path.isdir(args.output):
        raise ValueError(
            f"Output directory not found: {args.output}\n"
            "Orthosim.py should create this directory before this script is called."
        )

    config = read_config(args.config)
    args = configure_legacy_args(args, config)

    create_subfolders(args.o)
    check_inputs_and_load_parameters(args)

    # Override legacy hardcoded THREADS = 64 in the utils module.
    og.THREADS = args.threads

    print_detected_inputs(args)

    mp.set_start_method("spawn", force=True)

    og.RunOrthogroupMultiProc(
        n=args.n,
        outdir=args.o,
        threads=args.threads,
    )

    print("\n")

    og.CopyAlignments()
    print("Written alignments\n")

    og.CopyTrees()
    print("Written trees\n")

    og.BuildProteomes()
    print("Written proteomes\n")

    og.ConcatOrthologs()
    print("Written orthologs\n")

    og.SaveOrthogroups()
    print("Written orthogroups\n")

    og.BuildOrthogroupStats()
    print("Written orthogroup stats\n")

    og.WriteRunLog()
    print("Written log file\n")

    elapsed = time.time() - start_time
    print(f"Output written to: {args.o}")
    print(f"Elapsed time: {elapsed:.2f} seconds\n")
    print("Have a smashing day!\n")


if __name__ == "__main__":
    main()
