#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import time
import sys
import subprocess
import shutil

import SonicParanoid2_MakeOrthogroupSequenceFiles
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
    prog='SonicParanoid2Train.py',
    description=(
        "Welcome to OrthoTrain\n"
        "Train realistic orthogroup simulations from SonicParanoid2 outputs."
    ),
    formatter_class=argparse.RawDescriptionHelpFormatter
)

    # REQUIRED
    parser.add_argument(
        '-so',
        dest='sonic_output',
        type=str,
        required=True,
        help="Path to the base SonicParanoid2 results folder"
    )
    
    parser.add_argument(
        '-si',
        dest='sonic_input',
        type=str,
        required=True,
        help="Path to the base SonicParanoid2 input folder"
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

    sonic_out = os.path.abspath(args.sonic_output)
    sonic_in = os.path.abspath(args.sonic_input)
    output_folder_path = os.path.abspath(args.output)
    n_threads = args.threads


    # --- VALIDATION ---
    if not os.path.isdir(sonic_out):
        print(f"ERROR: SonicParanoid2 output folder not found: {sonic_out}")
        sys.exit(1)

    if not os.path.isdir(sonic_in):
        print(f"ERROR: SonicParanoid2 input folder not found: {sonic_in}")
        sys.exit(1)

    # ensure output exists
    os.makedirs(output_folder_path, exist_ok=True)

    # --- WELCOME ---
    print("---------------------------------------------------")
    print("Welcome to OrthoTrain\n")
    print(f"Base SonicParanoid2: {sonic_out}")
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
    print("Making SonicParanpoid2 sequence and alignment files ...")
    SonicParanoid2_MakeOrthogroupSequenceFiles.main(
        sonic_out,
        sonic_in,
        output_folder_path,
        n_threads,
        )
    print("Orthogroup sequence files created.")

    print("Aligning and Treeing to Estimate parameters.")
    EstimateParameters.main(    
            sonic_out,  
            output_folder_path,    
            n_threads)
    print("Parameter estimation complete.")
    print("Alignment and Tree parameters estimated.")
    print("---------------------------------------------------")

    print("Rooting gene trees and counting duplications.")
    CallDuplications.main(
        sonic_out,
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
