#!/usr/bin/env python3
"""
run_tools.py — Run OrthoFinder, FastOMA, SonicParanoid2, and Broccoli on a
simulated proteome dataset.

Usage:
    python run_tools.py --sim-folder /path/to/sim_folder [--threads 64]

What it does:
    1. Copies proteome_files/*.fa  →  proteome/*.fasta
    2. Runs OrthoFinder on proteome_files/
    3. Runs FastOMA       on proteome/  (via nextflow, conda env fastoma-nf)
    4. Runs SonicParanoid2 on proteome/ (conda env sonic)
    5. Runs Broccoli       on proteome/ (conda env env-broccoli)
       Broccoli's output dirs (dir_step1 … dir_step4) land in sim_folder.

Output locations (relative to --sim-folder):
    OrthoFinder   : proteome_files/OrthoFinder/Results_I1.5/
    FastOMA       : fastoma_output/
    SonicParanoid2: sp_default/
    Broccoli      : dir_step1/ … dir_step4/
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Machine-specific fixed paths — edit these once for your system
# ---------------------------------------------------------------------------

NEXTFLOW         = Path("/local/home/zool2506/SIM_TRAINING/OrthoTrain_analysis/nextflow")
OMAMER_DB        = Path("/local/home/zool2506/SIM_TRAINING/OrthoTrain_analysis/Metazoa.h5")
BROCCOLI_SCRIPT  = Path("/local/home/zool2506/SIM_TRAINING/OrthoTrain_analysis/Broccoli-master/broccoli.py")
NEXTFLOW_WORKDIR = Path("/local/home/zool2506/w1")
SPECIES_TREE = ("/local/home/zool2506/SIM_TRAINING/OrthoTrain_analysis/empiricaldata_qfoOrthobench/OrthoBench/fastoma_input/species_tree.nwk")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd, cwd=None):
    cmd = [str(x) for x in cmd]
    print(f"\n>>> {' '.join(cmd)}\n")
    subprocess.run(cmd, check=True, cwd=str(cwd) if cwd else None)


def run_in_conda(env_name, cmd, cwd=None):
    run(["conda", "run", "-n", env_name, "--no-capture-output"] + [str(x) for x in cmd], cwd=cwd)


def section(title):
    bar = "=" * 60
    print(f"\n{bar}\n{title}\n{bar}")


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

def prepare_proteomes(sim_dir):
    """Copy proteome_files/*.fa → proteome/*.fasta.

    OrthoFinder reads from proteome_files/ directly (.fa extension is fine).
    FastOMA, Broccoli, and SonicParanoid2 expect .fasta files in proteome/.
    """
    src = sim_dir / "proteome_files"
    dst = sim_dir / "proteome"

    if not src.is_dir():
        print(f"ERROR: proteome_files/ not found in sim folder: {src}")
        sys.exit(1)

    dst.mkdir(parents=True, exist_ok=True)
    fa_files = sorted(src.glob("*.fa"))

    if not fa_files:
        print(f"ERROR: no .fa files found in {src}")
        sys.exit(1)

    for f in fa_files:
        shutil.copy(f, dst / f"{f.stem}.fasta")
    ## copy species tree
    shutil.copy(SPECIES_TREE,sim_dir)
    
    print(f"Prepared {len(fa_files)} proteomes:")
    print(f"  Source : {src}")
    print(f"  Dest   : {dst}")
    return src, dst


# ---------------------------------------------------------------------------
# Tools
# ---------------------------------------------------------------------------

def run_orthofinder(src, threads):
    section("OrthoFinder")
    run([
        "orthofinder",
        "-f", src,
        "-t", threads,
        "-a", threads,
        "-I", "1.5",
        "-n", "I1.5",
    ])


def run_fastoma(sim_dir, threads):
    section("FastOMA")
    NEXTFLOW_WORKDIR.mkdir(parents=True, exist_ok=True)
    output_dir = sim_dir / "fastoma_output"
    run_in_conda("fastoma-nf", [
        NEXTFLOW, "run", "dessimozlab/FastOMA",
        "-r", "v0.4.0",
        "-profile", "conda",
        "-w", NEXTFLOW_WORKDIR,
        "--input_folder", sim_dir,
        "--output_folder", output_dir,
        "--omamer_db", OMAMER_DB,
        "--force_pairwise_ortholog_generation", "true",
        f"-process.cpus={threads}",
    ])


def run_sonicparanoid(sim_dir, dst, threads):
    section("SonicParanoid2")
    run_in_conda("sonic", [
        "sonicparanoid",
        "-i", dst,
        "-o", sim_dir / "sp_default",
        "-m", "default",
        "-t", threads,
    ])


def run_broccoli(sim_dir, dst, threads):
    section("Broccoli")
    # cwd=sim_dir so Broccoli's output dirs (dir_step1 … dir_step4) land there
    run_in_conda("env-broccoli", [
        "python", BROCCOLI_SCRIPT,
        "-dir", dst,
        "-threads", threads,
    ], cwd=sim_dir)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Run OrthoFinder, FastOMA, SonicParanoid2, and Broccoli on a simulated dataset.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--sim-folder", "-s",
        required=True,
        dest="sim_folder",
        help="Simulation output folder (must contain proteome_files/).",
    )
    parser.add_argument(
        "--threads", "-t",
        type=int,
        default=64,
        dest="threads",
        help="Number of threads to pass to each tool (default: 64).",
    )

    args = parser.parse_args()

    sim_dir = Path(args.sim_folder).resolve()
    if not sim_dir.is_dir():
        print(f"ERROR: --sim-folder not found: {sim_dir}")
        sys.exit(1)

    print(f"Sim folder : {sim_dir}")
    print(f"Threads    : {args.threads}")

    # Warn early if any fixed paths are missing
    for label, path in [
        ("NEXTFLOW",        NEXTFLOW),
        ("OMAMER_DB",       OMAMER_DB),
        ("BROCCOLI_SCRIPT", BROCCOLI_SCRIPT),
    ]:
        if not path.exists():
            print(f"WARNING: {label} not found: {path}")

    # Step 1: prepare proteomes
    src, dst = prepare_proteomes(sim_dir)

    # Steps 2–5: run tools
    run_orthofinder(src, args.threads)
    run_fastoma(sim_dir, args.threads)
    run_sonicparanoid(sim_dir, dst, args.threads)
    run_broccoli(sim_dir, dst, args.threads)

    section("Done")
    print(f"All tools finished. Outputs in: {sim_dir}")
    print()
    print("Expected output locations:")
    print(f"  OrthoFinder   : {src}/OrthoFinder/Results_I1.5/")
    print(f"  FastOMA       : {sim_dir}/fastoma_output/")
    print(f"  SonicParanoid2: {sim_dir}/sp_default/")
    print(f"  Broccoli      : {sim_dir}/dir_step1/ … dir_step4/")


if __name__ == "__main__":
    main()
