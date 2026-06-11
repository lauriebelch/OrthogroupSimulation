# OrthoSim

OrthoSim is a Python package for training and running realistic orthogroup simulations.

It is designed to generate simulated orthogroups that resemble empirical orthogroups while providing a complete ground truth for orthogroups and orthologs.

<p align="left">
  <img src="https://github.com/lauriebelch/OrthogroupSimulation/blob/main/OrthoSim.png" alt="OrthoSim logo" width="350">
</p>

> **Note:** This README is currently a work in progress.

## Main workflows

OrthoSim currently provides four main workflows.

| Workflow | Command | Purpose |
|---|---|---|
| Orthogroup simulation | `--Orthogroup-simulation` | Simulate orthogroups using a known set of parameters, provided in a config file. |
| Parameter extraction | `--Get-Parameters` | Extract parameters from empircal orthogroups (sequence-evolution, tree properties etc.) |
| Genetic algorithm tuning | `--GA` | Tune parameters that cannot be directly extracted from empirical orthogroups (duplication and loss rate, orthogroup birth bias). |
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

| Option | Meaning | Used by |
|---|---|---|
| `--output OUTPUT` | Path to the output directory | All workflows |
| `--threads THREADS` | Number of threads. Default is `1`. | All workflows |

## Workflow 1: Orthogroup simulation

The orthogroup simulation workflow generates simulated orthogroups using a supplied species tree, starting species, and simulation parameters.

### Example command

```bash
python OrthoSim.py \
  --Orthogroup-simulation \
  --output ExampleSimulation \
  --Config ExampleData/simulation_config.txt \
  --threads 4
```

Parameters can also be supplied directly on the command line:

```bash
python OrthoSim.py \
  --Orthogroup-simulation \
  --output ExampleSimulation \
  --ultrametric-tree ExampleData/species_tree.nwk \
  --species Mnemiopsis_leidyi \
  --Orthogroups 100 \
  --prop-invar-mean -1.2226 \
  --prop-invar-sd 0.4214 \
  --gamma-shape-mean 0.6443 \
  --gamma-shape-sd 0.3884 \
  --max-indel-insert 0.04 \
  --max-indel-delete 0.04 \
  --indel-size 1 \
  --max-duplication-rate 0.5 \
  --max-loss-rate 0.74 \
  --max-transfer-rate 0 \
  --replacement-prob 0.0 \
  --leaf-sampling-probability 1 \
  --relax-model ACRY07 \
  --max-start-rate 1 \
  --sigma-log-mean 0.972280 \
  --sigma-log-sd 0.524373 \
  --gbc 6.75 \
  --gap-file ExampleData/gap_position_profile_counts.tsv \
  --threads 4
```
**note that adding parameter options will overwrite the config file**

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
| `--max-indel-insert FLOAT` | Maximum insertion rate. Per-orthogroup insertion rates are sampled up to this value. |
| `--max-indel-delete FLOAT` | Maximum deletion rate. Per-orthogroup deletion rates are sampled up to this value. |
| `--indel-size FLOAT` | Indel-size parameter used by the sequence simulation. |
| `--max-duplication-rate FLOAT` | Maximum gene duplication rate. Usually tuned by the genetic algorithm. |
| `--max-loss-rate FLOAT` | Maximum gene loss rate. Usually tuned by the genetic algorithm. |
| `--max-transfer-rate FLOAT`| Maximum horizontal transfer rate. Default use is usually `0`. |
| `--replacement-prob FLOAT` | Transfer replacement probability. Default use is usually `0`. |
| `--leaf-sampling-probability FLOAT`| Probability of retaining/sampling leaves in the simulated gene tree. Usually `1`. |
| `--relax-model STR` | Branch-rate relaxation model, for example `ACRY07`. |
| `--max-start-rate FLOAT` | Starting branch-rate value for branch relaxation. Usually `1`. |
| `--sigma-log-mean FLOAT` | Mean of the log-space distribution used to sample branch-relaxation sigma. |
| `--sigma-log-sd FLOAT` | Standard deviation of the log-space distribution used to sample branch-relaxation sigma. |
| `--gbc FLOAT` | Orthogroup birth-bias / phylostratigraphy parameter. Usually tuned by the genetic algorithm. |
| `--gap-file PATH_TO_GAP_PROFILE` | Path to the empirical gap-position profile file. |


## Workflow 2: Get parameters

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

## Workflow 3: Genetic algorithm tuning

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

## Workflow 4: PFAM domain extraction

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
prop_invar_mean=-1.2226
prop_invar_sd=0.4214
gamma_shape_mean=0.6443
gamma_shape_sd=0.3884
max_indel_insert=0.04
max_indel_delete=0.04
indel_size=1
max_duplication_rate=0.5
max_loss_rate=0.74
max_transfer_rate=0
replacement_prob=0.0
leaf_sampling_probability=1
relax_model=ACRY07
max_start_rate=1
sigma_log_mean=0.972280
sigma_log_sd=0.524373
gbc=6.75
gap_file=ExampleData/gap_position_profile_counts.tsv
```

## Simulation config parameters

| Parameter | Meaning |
|---|---|
| `ultrametric_tree` | Path to the ultrametric species tree in Newick format. This tree provides the species framework for gene-tree simulation. |
| `species` | Starting or focal species used to seed simulated orthogroups with real protein sequences. |
| `Orthogroups` | Number of orthogroups to simulate. |
| `prop_invar_mean` | Mean parameter for the distribution of the proportion of invariant sites. This is estimated from empirical orthogroups using IQ-TREE results from trimmed alignments. |
| `prop_invar_sd` | Standard deviation parameter for the distribution of the proportion of invariant sites. |
| `gamma_shape_mean` | Mean parameter for the gamma shape distribution controlling among-site rate variation. This is estimated from empirical orthogroups. |
| `gamma_shape_sd` | Standard deviation parameter for the gamma shape distribution controlling among-site rate variation. |
| `max_indel_insert` | Maximum insertion rate. For each simulated orthogroup, the insertion rate is sampled between zero and this maximum. |
| `max_indel_delete` | Maximum deletion rate. For each simulated orthogroup, the deletion rate is sampled between zero and this maximum. |
| `indel_size` | Indel-size parameter used during sequence simulation. Current simulations use a geometric indel-size model. |
| `max_duplication_rate` | Maximum duplication rate used during gene-tree simulation. This is usually tuned by the genetic algorithm. |
| `max_loss_rate` | Maximum loss rate used during gene-tree simulation. This is usually tuned by the genetic algorithm. |
| `max_transfer_rate` | Maximum horizontal transfer rate. Current default use is usually `0`, but the parameter is adjustable. |
| `replacement_prob` | Probability of transfer replacement. Current default use is usually `0`, but the parameter is adjustable. |
| `leaf_sampling_probability` | Probability of retaining or sampling leaves in the simulated gene tree. Usually set to `1`. |
| `relax_model` | Branch-rate relaxation model used during simulation. For example, `ACRY07`. |
| `max_start_rate` | Starting branch-rate value for branch relaxation. Usually set to `1`. |
| `sigma_log_mean` | Mean of the log-space distribution used to sample the branch-relaxation sigma parameter. Estimated from empirical gene trees. |
| `sigma_log_sd` | Standard deviation of the log-space distribution used to sample the branch-relaxation sigma parameter. Estimated from empirical gene trees. |
| `gbc` | Orthogroup birth-bias / phylostratigraphy parameter controlling where orthogroups tend to originate on the species tree. Usually tuned by the genetic algorithm. |
| `gap_file` | Path to the empirical gap-position profile file. This is used to reposition simulated gap blocks so that simulated alignments better match empirical gap-position patterns. |

## How the simulations are generated

Explain how they are generated

## Output

## Benchmarking orthology inference methods

## Current notes and TODOs

## Citation

## Contact
