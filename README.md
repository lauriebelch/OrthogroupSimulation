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
|---|---|
| `--output` | Path to the output directory |
| `--threads` | Number of threads. Default is `1`.|

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
| `--output` | Yes | Output directory |
| `--Config` | Optional | Path to a config file containing simulation parameters. |
| `--threads` | Optional | Number of threads to use. Default is `1`. |

The options below are required unless provided in the [config file](#config-files). See [Simulation config parameters](#simulation-config-parameters) for an explanation of each parameter.

| Option | Meaning |
|---|---|
| `--ultrametric-tree` | Path to ultrametric species tree in Newick format. |
| `--species` | Starting species used to seed simulated orthogroups |
| `--Orthogroups` | Number of orthogroups to simulate. |
| `--prop-invar-mean` | Mean proportion of invariant sites. |
| `--prop-invar-sd ` | Standard deviation of the proportion of invariant sites. |
| `--gamma-shape-mean `  | Mean gamma shape distribution (controlling among-site rate variation.) |
| `--gamma-shape-sd ` | Standard deviation of gamma shape distribution. |
| `--max-indel-insert ` | Maximum insertion rate. |
| `--max-indel-delete ` | Maximum deletion rate. |
| `--indel-size ` | Indel-size parameter |
| `--max-duplication-rate ` | Maximum gene duplication rate. |
| `--max-loss-rate ` | Maximum gene loss rate. |
| `--max-transfer-rate `| Maximum horizontal transfer rate. |
| `--replacement-prob ` | Transfer replacement probability. |
| `--leaf-sampling-probability `| Probability of retaining leaves in the simulated gene tree. |
| `--relax-model ` | Branch-rate relaxation model, for example `ACRY07`. |
| `--max-start-rate ` | Starting branch-rate value for branch relaxation. |
| `--sigma-log-mean ` | Mean of the log-space distribution for branch-relaxation sigma. |
| `--sigma-log-sd ` | Standard deviation of the log-space distribution for branch-relaxation sigma. |
| `--gbc ` | Orthogroup birth-bias / phylostratigraphy parameter. |
| `--gap-file ` | Path to the empirical gap-position profile file. |


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

#### Required options

| Option | Allowed / default | Meaning |
|---|---|---|
| `--output` | NA | Output directory for extracted parameters. |
| `--tool` | `OF3`, `FASTOMA`, `BROCOLI`, `SONICPARANOID2` | Orthogroup inference tool used to generate the empirical orthogroups. |
| `--tool-output` | NA | Path to the output directory from the orthogroup inference run. |

#### Additional requirements by tool

| Tool | Additional required options | Meaning |
|---|---|---|
| `OF3` | None | OrthoFinder output contains the information needed for this workflow. |
| `FASTOMA` | `--tools-proteomes` | Input proteomes are needed to map genes back to species. |
| `BROCOLI` | `--tool-tree`, `--tools-proteomes` | The species tree and input proteomes are needed to process Broccoli output. |
| `SONICPARANOID2` | `--tool-tree`, `--tools-proteomes` | The species tree and input proteomes are needed to process SonicParanoid2 output. |

#### Optional options

| Option | Allowed / default | Meaning |
|---|---|---|
| `--threads` | Default: `1` | Number of threads to use. |

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
| `--output ` | Yes | Output directory for the genetic algorithm run. |
| `--threads` | Optional | Number of threads to use. Default is `1`. |

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
| `--output ` | Yes | Output directory for the PFAM run. |
| `--Proteome ` | Yes | Path to the proteome FASTA file. |
| `--threads ` | Optional | Number of threads to use. Default is `1`. |

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

OrthoSim simulates orthogroups by combining empirical parameter sampling with gene-tree and sequence simulation.

Each simulated orthogroup is seeded with a real protein sequence from the chosen starting species. A gene tree is then generated using SaGePhy, with duplication, loss, and orthogroup birth parameters taken from the simulation config file. Branch lengths are relaxed using SaGePhy branch-rate relaxation.

Protein sequence evolution is simulated using IQ-TREE AliSim. The background model uses empirical sequence-evolution parameters, including the proportion of invariant sites and gamma-distributed among-site rate variation.

OrthoSim divides each seed protein into PFAM-domain and non-domain regions. Domain regions are simulated using PFAM-specific substitution models and lower evolutionary rates, while non-domain regions use the background model.

## Output

## Benchmarking orthology inference methods

OrthoSim simulations can be used to benchmark orthology inference methods because the complete ground truth is known.

For each simulated dataset, OrthoSim records:
- the true orthogroup membership of every simulated gene
- the true ortholog pairs among all simulated genes
- the simulated gene trees used to generate the data

This means that the output of an orthology inference method can be compared directly against the true simulated relationships. 

## Citation

**All Aboard the OrthoTrain: empirically trained simulations for benchmarking orthology inference**
L. Belcher, J. H. Homes, and S. Kelly 

## Contact
