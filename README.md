# OrthoSim

OrthoSim is a Python package for training and running realistic orthogroup simulations.

It is designed to generate simulated orthogroups that resemble empirical orthogroups while providing a complete ground truth for orthogroups and orthologs.

## Contents

- [Main workflows](#main-workflows)
- [Pipeline summary](#pipeline-summary)
- [Basic usage](#basic-usage)
- [Command-line options](#command-line-options)
  - [General options](#general-options)
  - [Orthogroup simulation workflow](#orthogroup-simulation-workflow)
  - [Get-parameters workflow](#get-parameters-workflow)
  - [Genetic algorithm workflow](#genetic-algorithm-workflow)
  - [PFAM domain extraction workflow](#pfam-domain-extraction-workflow)
- [Config files](#config-files)
- [Simulation config parameters](#simulation-config-parameters)
- [How the simulations are generated](#how-the-simulations-are-generated)
- [Output](#output)

<p align="left">
  <img src="https://github.com/lauriebelch/OrthogroupSimulation/blob/main/OrthoSim.png" alt="OrthoSim logo" width="350">
</p>

> **Note:** This README is currently a work in progress.

## Main workflows

OrthoSim currently provides four main workflows.

| Workflow | Command | Purpose |
|---|---|---|
| Parameter extraction | `--Get-Parameters` | Extract parameters from empircal orthogroups (sequence-evolution, tree properties etc.) |
| Genetic algorithm tuning | `--GA` | Tune parameters that cannot be directly extracted from empirical orthogroups (duplication and loss rate, orthogroup birth bias). |
| Orthogroup simulation | `--Orthogroup-simulation` | Simulate orthogroups using a known set of parameters, provided in a config file. |
| PFAM domain extraction | `--PFAM` | Extract PFAM domain information for a proteome, allowing simulations to work differently on domain and non-domain regions. (note that we provide PFAM domain info for some model organisms) |

## Pipeline summary

The intended workflow is:

```text
Empirical proteomes
        v
Orthogroup inference
        v
OrthoSim --Get-Parameters
        v
OrthoSim --GA
        v
[OPTIONAL] OrthoSim --PFAM (to extract domain info for non-model species
        v
OrthoSim --Orthogroup-simulation
        v
OUTPUT = Simulated orthogroups with complete ground truth
```

## Basic usage

Run the OrthoSim.py script with one workflow flag and the required options for that workflow.

```bash
python OrthoSim.py --Orthogroup-simulation [options]
```

## Command-line options

### General options

These options are shared by multiple workflows.

| Option | Meaning |
|---|---|---|
| `--output OUTPUT` | Path to the output directory |
| `--threads THREADS` | Number of threads. Default is `1`.|

## Orthogroup simulation workflow

This workflow generates simulated orthogroups using a supplied species tree, starting species, and simulation parameters.

### Example command

```bash
python OrthoSim.py \
  --Orthogroup-simulation \
  --output ExampleSimulation \
  --Config ExampleData/simulation_config.txt \
  --threads 4
```

Parameters can also be supplied directly on the command line. If the same parameter is present both in the config file and as a command-line flag, the command-line value is used.

For example:

```bash
python OrthoSim.py \
  --Orthogroup-simulation \
  --output ExampleSimulation \
  --Config ExampleData/simulation_config.txt \
  --ultrametric-tree ExampleData/species_tree.nwk \
  --species Mnemiopsis_leidyi \
  --max-duplication-rate 0.5 \
  # remaining options ...
  --threads 4
```

### Orthogroup simulation command-line options

| Option | Required? | Meaning |
|---|---:|---|
| `--output OUTPUT` | Yes | Output directory |
| `--Config PATH_TO_CONFIG_FILE` | Optional | Path to a config file containing simulation parameters. |
| `--threads THREADS` | Optional | Number of threads to use. Default is `1`. |

These options are required unless provided in the config file
| Option | Meaning |
|---|---|
| `--ultrametric-tree TREE_PATH` | Path to the input ultrametric species tree in Newick format. |
| `--species SPECIES_NAME` | Starting/focal species used to seed simulated orthogroups with real protein sequences. |
| `--Orthogroups NUMBER_OF_ORTHOGROUPS` | Number of orthogroups to simulate. |
| `--prop-invar-mean FLOAT` | Mean parameter for the distribution of the proportion of invariant sites. |
| `--prop-invar-sd FLOAT` | Standard deviation parameter for the distribution of the proportion of invariant sites. |
| `--gamma-shape-mean FLOAT`  | Mean parameter for the gamma shape distribution controlling among-site rate variation. |
| `--gamma-shape-sd FLOAT` | Standard deviation parameter for the gamma shape distribution. |
| `--max-indel-insert FLOAT` | Maximum insertion rate. |
| `--max-indel-delete FLOAT` | Maximum deletion rate. |
| `--indel-size FLOAT` | Indel-size parameter |
| `--max-duplication-rate FLOAT` | Maximum gene duplication rate. |
| `--max-loss-rate FLOAT` | Maximum gene loss rate. |
| `--max-transfer-rate FLOAT`| Maximum horizontal transfer rate. |
| `--replacement-prob FLOAT` | Transfer replacement probability. |
| `--leaf-sampling-probability FLOAT`| Probability of retaining/sampling leaves in the simulated gene tree. |
| `--relax-model STR` | Branch-rate relaxation model, for example `ACRY07`. |
| `--max-start-rate FLOAT` | Starting branch-rate value for branch relaxation. Usually `1`. |
| `--sigma-log-mean FLOAT` | Mean of the log-space distribution used to sample branch-relaxation sigma. |
| `--sigma-log-sd FLOAT` | Standard deviation of the log-space distribution used to sample branch-relaxation sigma. |
| `--gbc FLOAT` | Orthogroup birth-bias / phylostratigraphy parameter. |
| `--gap-file PATH_TO_GAP_PROFILE` | Path to the empirical gap-position profile file. |


## Get-parameters workflow

The `--Get-Parameters` workflow extracts empirical parameters from orthogroups inferred by an orthogroup inference tool.

### Example command

```bash
python OrthoSim.py \
  --Get-Parameters \
  --output ParameterRun \
  --tool OF3 \
  --tool-output ExampleData/OrthoFinderResults \
  --tool-tree ExampleData/species_tree.nwk \
  --tools-proteomes ExampleData/Proteomes \
  --threads 4
```
### Get-Parameters command-line options

| Option | Required? | Meaning |
|---|---:|---|
| `--output OUTPUT` | Yes | Output directory for extracted parameters. |
| `--tool OF3,FASTOMA,BROCOLI,SONICPARANOID2` | Yes | Orthogroup inference tool used to generate the empirical orthogroups. |
| `--tool-output PATH_TO_TOOL_OUTPUT` | Yes | Path to the output directory or results file from the orthogroup inference tool. |
| `--tool-tree PATH_TO_TOOL_TREE` | Required for Broccoli and SonicParanoid2 | Path to the species tree associated with the orthogroup inference run. |
| `--tools-proteomes PATH_TO_TOOLS_PROTEOMES` | Required for FastOMA, Broccoli, and SonicParanoid2 | Path to the proteomes used as input to the orthogroup inference tool. |
| `--threads THREADS` | Optional | Number of threads to use. Default is `1`. |

## Genetic algorithm workflow

The `--GA` workflow tunes parameters that cannot be directly estimated from empirical orthogroup assignments.

### Example command

```bash
python OrthoSim.py \
  --GA \
  --threads 4
```

### Genetic algorithm command-line options

| Option | Required? | Meaning |
|---|---:|---|
| `--output OUTPUT` | Yes | Output directory for the genetic algorithm run. |
| `--threads THREADS` | Optional | Number of threads to use. Default is `1`. |

## PFAM domain extraction workflow

The `--PFAM` workflow extracts PFAM domain information for a given proteome.

PFAM annotations are used by the simulation workflow so that domain and non-domain regions can be simulated separately. Domain regions can be assigned PFAM-specific substitution models and lower evolutionary rates, while non-domain regions use the background sequence-evolution model.

### Example command

```bash
python OrthoSim.py \
  --PFAM \
  --output PFAMRun \
  --Proteome ExampleData/Mnemiopsis_leidyi.fa \
  --threads 4
```

### PFAM command-line options

| Option | Required? | Meaning |
|---|---:|---|
| `--output OUTPUT` | Yes | Output directory for the PFAM run. |
| `--Proteome PATH_TO_PROTEOME` | Yes | Path to the proteome FASTA file. |
| `--threads THREADS` | Optional | Number of threads to use. Default is `1`. |

## Config files

The orthogroup simulation workflow can be run using a config file.

Config files use simple `key=value` formatting.

Example:

```text
ultrametric_tree=ExampleData/species_tree.nwk
species=Mnemiopsis_leidyi
Orthogroups=100
```

## Simulation config parameters

| Parameter | Meaning |
|---|---|
| `ultrametric_tree` | Path to the ultrametric species tree in Newick format. |
| `species` | Starting species used to seed simulated orthogroups with real protein sequences. |
| `Orthogroups` | Number of orthogroups to simulate. |
| `prop_invar_mean` | Mean parameter for the distribution of the proportion of invariant sites. |
| `prop_invar_sd` | Standard deviation parameter for the distribution of the proportion of invariant sites. |
| `gamma_shape_mean` | Mean parameter for the gamma shape distribution controlling among-site rate variation. |
| `gamma_shape_sd` | Standard deviation parameter for the gamma shape distribution controlling among-site rate variation. |
| `max_indel_insert` | Maximum insertion rate. |
| `max_indel_delete` | Maximum deletion rate. |
| `indel_size` | Indel-size parameter used during sequence simulation. |
| `max_duplication_rate` | Maximum duplication rate used during gene-tree simulation. |
| `max_loss_rate` | Maximum loss rate used during gene-tree simulation. |
| `max_transfer_rate` | Maximum horizontal transfer rate. |
| `replacement_prob` | Probability of transfer replacement. |
| `leaf_sampling_probability` | Probability of retaining or sampling leaves in the simulated gene tree. |
| `relax_model` | Branch-rate relaxation model used during simulation. For example, `ACRY07`. |
| `max_start_rate` | Starting branch-rate value for branch relaxation. |
| `sigma_log_mean` | Mean of the log-space distribution used to sample the branch-relaxation sigma parameter.  |
| `sigma_log_sd` | Standard deviation of the log-space distribution used to sample the branch-relaxation sigma parameter. |
| `gbc` | Orthogroup birth-bias / phylostratigraphy parameter controlling where orthogroups tend to originate on the species tree. |
| `gap_file` | Path to the empirical gap-position profile file. This is used to reposition simulated gap blocks so that simulated alignments better match empirical gap-position patterns. |

## How the simulations are generated

Explain how they are generated

## Output

## Benchmarking orthology inference methods

## Current notes and TODOs

## Citation

## Contact
