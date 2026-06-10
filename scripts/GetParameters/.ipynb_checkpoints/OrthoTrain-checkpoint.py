#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import time
import sys
import subprocess
import shutil

import MakeOrthogroupSequenceFiles
import CallDuplications
import EstimateParameters


# Ensure the script's directory is in PYTHONPATH
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

### python scripts/OrthoTrain.py -f Results_Apr15 -o output -pf pfam.txt -t 8 -DO --assign model_proteome/

def main():
    TRAIN_ART = r"""
            (  ) (@@) ( )  (@)  ()    @@    O     @     O     @
       (@@@)
   (    )
(@@@@)

     ====        ________                ___________
 _D _|  |_______/        \__I_I_____===__|_________|
  |(_)---  |   H\________/ |   |        =|___ ___|     
  /     |  |   H  |  |     |   |         ||_| |_||     
 |      |  |   H  |__--------------------| [___] |  
 | ________|___H__/__|_____/[][]~\_______|       |  
 |/ |   |-----------I_____I [][] []  D   |=======|
__/ =| o |=-~~\  /~~\  /~~\  /~~\ ____Y___________|
 |/-=|___|=O=====O=====O=====O   |_____/~\___/
  \_/      \__/  \__/  \__/  \__/

        All aboard the OrthoTrain!
"""
    print(TRAIN_ART)
    parser = argparse.ArgumentParser(
    prog='OrthoTrain.py',
    description=(
        "Welcome to OrthoTrain\n"
        "Train realistic orthogroup simulations from OrthoFinder outputs."
    ),
    epilog=(
        "Example usage:\n"
        "  OrthoTrain.py -f Results_Apr15 -o output-t 8 \n"
        "  OrthoTrain.py -f Results_Apr15 -o output -t 8 \n"
    ),
    formatter_class=argparse.RawDescriptionHelpFormatter
)

    # REQUIRED
    parser.add_argument(
        '-f',
        dest='folder',
        type=str,
        required=True,
        help="Path to the base OrthoFinder results folder"
    )

    parser.add_argument(
        '-o',
        dest='output',
        type=str,
        required=True,
        help="Output folder"
    )

    parser.add_argument(
        '-t',
        dest='threads',
        type=int,
        default=8,
        help="Number of threads"
    )
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    ortho_folder_path = os.path.abspath(args.folder)
    output_folder_path = os.path.abspath(args.output)
    n_threads = args.threads

    # --- VALIDATION ---
    if not os.path.isdir(ortho_folder_path):
        print(f"ERROR: OrthoFinder folder not found: {ortho_folder_path}")
        sys.exit(1)

    # ensure output exists
    os.makedirs(output_folder_path, exist_ok=True)

    def get_single_fasta(folder):
        fasta_files = [
            f for f in os.listdir(folder)
            if f.endswith(".fa") or f.endswith(".fasta")
        ]
        if len(fasta_files) == 0:
            print(f"ERROR: Model proteome folder must contain exactly one FASTA file: {folder}")
            sys.exit(1)
        if len(fasta_files) > 1:
            print(f"ERROR: Model proteome folder must contain exactly one FASTA file: {folder}")
            sys.exit(1)
        return os.path.join(folder, fasta_files[0])

    # --- WELCOME ---
    print("---------------------------------------------------")
    print("Welcome to OrthoTrain\n")
    print(f"Base OrthoFinder: {ortho_folder_path}")
    #print(f"Domain OrthoFinder: {domain_results_path}")
    print(f"Output: {output_folder_path}")

    start_time = time.time()

    current_hour = time.localtime().tm_hour
    if 6 <= current_hour < 12:
        print("Good morning!")
    elif 12 <= current_hour < 18:
        print("Good afternoon!")
    else:
        print("Good evening!")

    print("---------------------------------------------------")


    # OrthoFinder logic
    print("Making Orthogroup sequence files ...")
    MakeOrthogroupSequenceFiles.main(
        ortho_folder_path,
        output_folder_path,
        n_threads
    )
    print("Orthogroup sequence files created.")
    print("---------------------------------------------------")

    print("Aligning and Treeing to Estimate parameters.")
    done_file = os.path.join(output_folder_path,"TrainingResults","Orthogroup_analysis","DONE")
    if os.path.exists(done_file):
        print("Parameter estimation already completed — skipping.")
    else:
        EstimateParameters.main(
            ortho_folder_path,
            output_folder_path,
            n_threads
        )
    print("Alignment and Tree parameters estimated.")
    print("---------------------------------------------------")

    print("Rooting gene trees and counting duplications.")
    CallDuplications.main(
        ortho_folder_path,
        output_folder_path,
        n_threads,
    )
    print("Gene tree rooting and duplication counting complete.")
    print("---------------------------------------------------")

    print("Preparing SimulationInputs...")
    sim_dir = os.path.join(output_folder_path, "TrainingResults", "SimulationInputs")
    os.makedirs(sim_dir, exist_ok=True)
    def safe_copy(src, dst_name):
        if not os.path.exists(src):
            print(f"WARNING: missing file: {src}")
            return
        dst = os.path.join(sim_dir, dst_name)
        shutil.copy2(src, dst)
    # --- key inputs from pipeline ---
    species_tree = os.path.join(
        output_folder_path,
        "TrainingResults",
        "Orthogroup_analysis",
        "rescaled_ultrametric_species_tree.nwk"
    )
    sim_params = os.path.join(
        output_folder_path,
        "TrainingResults",
        "Orthogroup_analysis",
        "simulation_parameters.txt"
    )
    gap_profile = os.path.join(
        output_folder_path,
        "TrainingResults",
        "Orthogroup_analysis",
        "gap_position_profile_counts.tsv"
    )
    # --- copy everything ---
    safe_copy(species_tree, "species_tree.nwk")
    safe_copy(sim_params, "simulation_parameters.txt")
    safe_copy(gap_profile, "gap_position_profile_counts.tsv")
    # keep original filenames for these
    print(f"SimulationInputs ready at: {sim_dir}")
    print("---------------------------------------------------")

    print("Done.")

    elapsed_time = time.time() - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)

    print(f"---Finished! This run took {minutes} minutes {seconds} seconds.")
    print(f"---Files have landed in {output_folder_path}")
    print("---Thank you for choosing OrthoTrain.")
    print("---OrthoTrain: Kelly Lab (2026), bioRxiv ")


if __name__ == "__main__":
    main()
