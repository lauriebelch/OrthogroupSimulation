 #!/usr/bin/env python3

import argparse
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import subprocess
import uuid
import shutil
import argparse
# --- HARDCODED PATHS ---
#SIM_INPUT_FOLDER = "/local/home/biol0216/Simulation/Tuning/GA/GA/sim_training_v1/scripts/SimulationInputs"
#TRAINING_ANALYSIS = "/local/home/biol0216/Simulation/Tuning/GA/GA/sim_training_v1/scripts/Orthogroup_analysis"

#SIM_INPUT_FOLDER = "/local/home/zool2506/SIM_TRAINING/OrthoTrain/OrthoBench/fastoma_train/TrainingResults/SimulationInputs"
#TRAINING_ANALYSIS = "/local/home/zool2506/SIM_TRAINING/OrthoTrain/OrthoBench/fastoma_train/TrainingResults/Orthogroup_analysis"

#SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
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


def modify_param_file(complete_parameters,param_file,current_file_path,default_species_pick,run_dir, dup_rate, loss_rate, gbc_val):
    new_param_file = os.path.join(run_dir,"parameter_file.txt")
    with open(param_file) as f:
        lines = f.readlines()
    ######### needs...
    """
    ultrametric_tree=Training_FastOMA/TrainingResults/SimulationInputs/species_tree.nwk
    species=Mnemiopsis_leidyi
    Orthogroups=100
    threads=4
    """

    with open(new_param_file, "w") as f:
        tree_path = "ultrametric_tree=" + os.path.join(complete_parameters["Orthogroup_train_results"],"TrainingResults","SimulationInputs","species_tree.nwk")
        species = "species="+ default_species_pick
        Orthogroups  = "Orthogroups=" + complete_parameters['Orthogroups']
        gap_file=os.path.join(complete_parameters["Orthogroup_train_results"],"TrainingResults","SimulationInputs","gap_position_profile_counts.tsv")
        threads = "threads=1"
        output = "SimPath=" + run_dir
        f.write(f"{tree_path}\n")
        f.write(f"{species}\n")
        f.write(f"{Orthogroups}\n")
        f.write(f"{threads}\n")
        f.write(f"{output }\n")

        for line in lines:
            if line.startswith("max_duplication_rate="):
                f.write(f"max_duplication_rate={dup_rate}\n")
            elif line.startswith("max_loss_rate="):
                f.write(f"max_loss_rate={loss_rate}\n")
            elif line.startswith("gbc="):
                f.write(f"gbc={gbc_val}\n")
            elif line.startswith("sagephy_path"):
                new_sagephy_path = os.path.join(run_dir,"bin","sagephy-1.0.0.jar")
                new_sagephy_path_command = "sagephy_path=" + new_sagephy_path
                f.write(new_sagephy_path_command+"\n")

            elif line.startswith("iqtree"):
                new_iqtree_path = os.path.join(run_dir,"bin","iqtree3")
                new_iqtree_path_command = "iqtree_path=" + new_iqtree_path
                f.write(new_iqtree_path_command+"\n")

            else:
                f.write(line)



def main(dup_rate,loss_rate,gbc_val,complete_parameters,current_file_path,default_species_pick):
    #parser = argparse.ArgumentParser()
    #parser.add_argument("-d", type=float, required=True)
    #parser.add_argument("-l", type=float, required=True)
    #parser.add_argument("-gbc", type=float, required=True)

    #args = parser.parse_args()

    #dup_rate = args.d
    #loss_rate = args.l
    #gbc_val = args.gbc


    # --- create run folder ---
    base_out = os.path.abspath(os.path.join(complete_parameters["output"],"Simulation_Temp_Files"))
    #ORTHOTRAIN_DIR = os.path.dirname(SCRIPT_DIR)
    #base_out = os.path.join(ORTHOTRAIN_DIR, "simulation_runs")    
    os.makedirs(base_out, exist_ok=True)

    run_dir = os.path.join(base_out, uuid.uuid4().hex)
    os.makedirs(run_dir)


    # --- modify parameter file INSIDE copied folder ---
    param_file = os.path.join(complete_parameters["Orthogroup_train_results"],"TrainingResults","SimulationInputs", "simulation_parameters.txt")
    modify_param_file(complete_parameters,param_file,current_file_path,default_species_pick,run_dir, dup_rate, loss_rate, gbc_val)


    # --- paths to scripts ---
    #simulate_script = os.path.join(SCRIPT_DIR, "simulate_orthogroup.py")
    #stats_script = os.path.join(SCRIPT_DIR, "GetSimulationStats.py")
    #compare_script = os.path.join(SCRIPT_DIR, "CompareDistributions.py")

    # --- copy binary files to new location ---
    binary_bin = os.path.join(current_file_path,"bin")
    run_bin = os.path.join(run_dir,"bin")
    shutil.copytree(binary_bin,run_bin)


    #sim_out = os.path.join(run_dir, "sim")
	
    # --- run simulation ---
    #cmd_sim = [
    #    "python", simulate_script,
    #    "--folder", sim_inputs_copy,
    #    "-n", "500",
    #    "-o", sim_out,
    #    "-d", "Mnemiopsis_leidyi"
    #]
    #subprocess.run(cmd_sim, check=True, cwd=SCRIPT_DIR)

    # --- run simulation ---
    import scripts.OrthogroupSimulation.orthogroup_simulation_utils as og

    config_file = read_config(os.path.join(run_dir, "parameter_file.txt"))
    print(config_file)
		

    quit()






    # --- stats ---
    cmd_stats = [
        "python", stats_script,
        "-f", sim_out
    ]

    subprocess.run(cmd_stats, check=True, cwd=SCRIPT_DIR)

    # --- compare ---
    cmd_compare = [
        "python", compare_script,
        "-s", sim_out,
        "-e", TRAINING_ANALYSIS
    ]

    result = subprocess.run(
        cmd_compare,
        capture_output=True,
        text=True,
        check=True,
        cwd=SCRIPT_DIR
    )

    # --- extract score ---
    score = result.stdout.strip().splitlines()[-1]

    return score,run_dir


if __name__ == "__main__":
    main()