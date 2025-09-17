# OrthoSim

OrthoSim is a realistic orthogroup simulator, providing proteomes where we know exactly which genes are in which orthogroup, and whether they are orthologs or paralogs

## How do I use it?

usage: simulate_orthogroup.py [-h] -s -f -PFAM -d -n -p -o

You must provide the following inputs

```-s``` is the species tree - which must be ultrametric and have species names in the format genus_species

```-f``` is the starting genome (protein fasta)

```-n``` is the number of orthogroups you wish to simulate

```-p``` is the customisable parameter file

```-o``` is the output directory

```-PFAM``` is a pfam_scan output file containing information on domains and motifs in the genes in your starting genome 

```-d``` is a list of substitution models to be used in PFAM domain regions. It must be a csv file with a column 'alisim_model'


You can test Orthogroup simulations by first downloading the script and the example data in the 'inputs' folder

```python simulate_orthogroup.py -s inputs/ultrametric_fungi_tree.nwk -f inputs/Saccharomyces_cerevisiae.fasta -p inputs/simulation_parameters.txt -o fun1 -n 100 -PFAM inputs/fungi.pfam.txt -d inputs/domain_models.csv```

## What outputs do I get?

- alignment_files = a multiple sequence alignment for each simulated orthogroup
- tree_files = a gene tree for each simulated orthogroup
- proteome_files = a proteome for each species in the species tree
- Simulated_orthogroups = a file listing genes in each simulated orthogroup
- Simulated_orthologs = ortholog pairs
- Simulated_orthogroup_statistics = summary file showing the parameters selected for each orthogroup
- temporary_files = intermediate files generated during the simulations

Genes in proteomes are named in the format genus_species_og_index


## How are the simulations generated?

For each orthogroup, the basic approach is as follows;
1) select a random gene from the starting genome
2) simulate a gene tree with duplications/losses etc (SagePhy)
3) relax the branch lengths on that gene tree (SagePhy)
4) simulate an alignment with different substitution models for domain and non-domain regions
5) add each gene to the relevant proteome

For full details, please read our paper

## How do I get a PFAM domain file?

Run pfam scan [https://anaconda.org/bioconda/pfam_scan] on your starting genome


## How do I get substitution models for domains?

The reccommended way is to estimate real parameters from an Orthofinder run

First run OrthoFinder3.1 on some proteomes, then use the scripts provided in the 'tools' folder to generate the file you need

```python domain_param_estimation.py```
```python model_reader.py```

## What do all of the variables in the simulation_parameters file mean?

Lets go through the example:

*Proportion of invariants sites for non-domain regions are drawn from a log-normal distribution. Note that these values are defined on the log-scale.*

- prop_invar_mean=-2.65120136307523

- prop_invar_sd=1.06026766870849

*Gamma shape alphas for the non-domain regions are drawn from a log-normal distribution. Note that these values are defined on the log-scale.*

- gamma_shape_mean=0.2252

- gamma_shape_sd=0.7189

*Insertion and deletion rates are drawn from a uniform distribution with a maximum of the value you provide. The units are 'number of in/dels per 100 substitutions'*

- max_indel_insert=0.04

- max_indel_delete=0.04

*Indel sizes are drawn from a mean of the value you provide, and a variance determined by iqtree's alisim*

- indel_size=1

*Duplication, loss, and transfer rates are passed to SagePhy. they have no units*

- max_duplication_rate=0.5

- max_loss_rate=0.74

- max_transfer_rate=0

*When a transfer occurs, it can either add to the genome, or replace another gene*

- replacement_prob=0.0

*You might want the simulations to miss some genes*

- leaf_sampling_probability=1

*Sagephy branch length relaxation allows several models and rates*

- relax_model=ACRY07

- max_start_rate=0.5

- max_sigma2=0.1

*gbc determines how biased orthogroup birth is to the root of the species tree. if you set this as 50000 almost all orthogroups will be born at the root*

- gbc=6.75

*you must provide paths to sagephy and iqtree. Note that iqtree must be the version that allows indels with partition files [https://github.com/iqtree/iqtree3/issues/87]*

- sagephy_path=/local/home/zool2506/Simulations/python_sim/SaGePhy-master/sagephy-1.0.0.jar

- iqtree_path=/local/home/zool2506/Simulations/python_sim/iqtree-3.0.1-Linux/bin/iqtree3


## Can I use these simulations to benchmark my orthology inference method?

Yes. Ask Yi for details

## All Orthogroup simulation options

options:
  -h, --help  show this help message and exit
  -s S        Species tree (newick
  -f F        Starting genome (fasta)
  -PFAM PFAM  PFAM_scan results
  -d D        Domain substitution models
  -n N        Number of Orthogroups
  -p P        Parameter file
  -o O        Output directory
