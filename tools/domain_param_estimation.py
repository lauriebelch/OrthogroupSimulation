#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 07:42:08 2025

@author: user
"""
import os
import subprocess
import re

alignments_folder = "/local/home/zool2506/Simulations/OPTIMIZING_SIMS/FUNGI/orthofinder_input/OrthoFinder/Results_Jul30/MultipleSequenceAlignments"
output_folder = "/local/home/zool2506/Simulations/python_sim/domain_param_estimation/Subset_MSA"
genome_file = "/local/home/zool2506/Simulations/python_sim/inputs/Saccharomyces_cerevisiae.fasta"

def ReadGenome(filename):
    genome = {}
    current_id = None
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:]
                current_id = current_id.split()[0]
                genome[current_id] = ""
            else:
                genome[current_id] += line
    return genome

def ReadMSA(filename):
    msa = {}
    current_id = None
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:]
                msa[current_id] = ""
            else:
                msa[current_id] += line
    if len(msa.keys()) < 4:
        return 0
    return msa

def ExtractPfamDomains(msa):
    pfam_file = "/local/home/zool2506/Simulations/python_sim/fungi.pfam.txt"
    # set gene_id 
    # filter alignment to only s. cerevisiae genes    
    filtered_ids = {k.split()[0] for k in msa.keys() if "Saccharomyces_cerevisiae" in k}
    filtered_ids2 = {gene_id.rsplit("_", 1)[-1] for gene_id in filtered_ids}
    if not filtered_ids:
        return 0    
    # results
    results = []
    types = {"Domain"}
    with open(pfam_file) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split()
            if cols[0] not in filtered_ids2:
                continue
            if types and cols[7] not in types:
                continue
            hit = {
                "seq_id": "Saccharomyces_cerevisiae_"+cols[0],
                "start": int(cols[1]),
                "end": int(cols[2]),
            }
            results.append(hit)
            return(results)
    return 0
    
def GetFocalProtein(res):
    a = res[0]["seq_id"]
    a = a.rsplit("_", 1)[-1]
    return genome[a]

def GetGeneName(res):
    a = res[0]["seq_id"]
    a = a.rsplit("_", 1)[-1]
    return a

def MapAlignmentBLAST(raw_seq, alignment, domain, output_folder):
    seq_id = domain[0]["seq_id"]
    aligned_seq = alignment[seq_id]
    start = domain[0]["start"]
    end   = domain[0]["end"]
    
    ### use blast to map start and end to the alignment 
    fasta_path = os.path.join(output_folder, f"{seq_id}_ungapped.fasta")
    db_prefix  = os.path.join(output_folder, f"{seq_id}_db")
    query_path = os.path.join(output_folder, f"{seq_id}_query.fasta")
    out_path   = os.path.join(output_folder, f"{seq_id}_blast.txt")
    # write ungapped sequence from alignment
    record_postions_of_gaps = []
    for i in range(0,len(aligned_seq)):
        if aligned_seq[i] == "-":
            record_postions_of_gaps.append(i)
    ungapped_seq = "".join(c for c in aligned_seq if c != "-")
    with open(fasta_path, "w") as f:
        f.write(f">{seq_id}\n{ungapped_seq}\n")
    os.system("diamond makedb --quiet --in %s -d %s" % (fasta_path, db_prefix))
    # write query (domain subsequence)
    query_seq = raw_seq[start-1:end]
    with open(query_path, "w") as f:
        f.write(">query\n%s\n" % query_seq)
    # run diamond blastp (tabular outfmt)
    command = "diamond blastp -d %s -q %s -o %s --outfmt 5 --quiet --evalue 0.05 -k 0 --max-hsps 0 --threads 1" % (db_prefix, query_path, out_path)
    os.system(command)
    # parse first hit
    domain_hits = []
    alignment_locations = []
    with open(out_path) as results_file:
            for line in results_file:
                if "<Hsp_hseq>" in line:
                    domain_hits.append(line.split(">")[1].split("<")[0].replace("-",""))
                    break
    START = 0
    END = 0
    ungapped = ungapped_seq
    for sub_seq in domain_hits:
        if sub_seq in ungapped:
            START = ungapped.index(sub_seq)
            END = ungapped.index(sub_seq) + len(sub_seq)
    var = 0
    for i in record_postions_of_gaps:
        if i < START:
            START = START + 1
            var = var + 1
        if START < i and i < END:
            END = END + 1
    alignment_locations.append([START,(END+var)])
    g_start = START
    g_end = END+var
    if not g_start:
        return 0,0
    if not g_end:
        return 0,0
    return g_start, g_end

def SubsetAlignment(alignment, domain, raw_protein, output_folder):
    region_start, region_end = MapAlignmentBLAST(raw_protein, alignment, domain, output_folder)
    if region_end == 0:
        return 0
    #print(region_start)
    #print(region_end)
    return {k: v[region_start:region_end+1] for k, v in alignment.items()}

def SaveSubsetAlignment(alignment, filename):
    with open(filename, "w") as f:
        for header, seq in alignment.items():
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
                
def RunIQtree(fasta):
    #cmd = ["iqtree3", "--quiet", "-s", fasta, "-m", "LG+G4+I",
     #      "-T", "128", "--fast"]
    trimmed_file = outname+".fasta"
    trimal_cmd_trim = ["trimal", "-in", fasta, "-out", trimmed_file, "-automated1"]
    subprocess.run(trimal_cmd_trim, check=True)
    if not os.path.exists(trimmed_file):
        return 0
    result_file  = trimmed_file + ".iqtree"   # IQ-TREE always writes this
    if os.path.exists(result_file):
        print(f"Skipping {fasta} (already processed)")
        return
    # remove sequences that are all gaps
    cleaned = {}
    with open(trimmed_file) as f:
        current = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line
                cleaned[current] = ""
            else:
                cleaned[current] += line
    cleaned = {h:s for h,s in cleaned.items() if s.replace("-", "") != ""}
    if len(cleaned)<3:
        return
    with open(trimmed_file, "w") as f:
        for h, s in cleaned.items():
            f.write(f"{h}\n")
            f.write(f"{s}\n")
    cmd = ["iqtree3", "--quiet", "-s", trimmed_file, "-mset", "LG+I,JTT+G",
           "-T", "18","--cmax", "4", "--fast"]
    subprocess.run(cmd, check=True)

# this will be the main code

genome = ReadGenome(genome_file)

for fname in os.listdir(alignments_folder):
    if not fname.endswith(".fa"):
        continue
    alignment_file = os.path.join(alignments_folder, fname)
    print(f"\nProcessing {alignment_file}")
    # load the MSA
    msa = ReadMSA(alignment_file)
    if msa == 0:
        print("  skipped: too few sequences")
        continue
    # look for pfam domains
    res = ExtractPfamDomains(msa)
    if res == 0:
        print("  skipped: no matching Pfam domains")
        continue
    # subset alignment
    my_protein = GetFocalProtein(res)
    gene_name = GetGeneName(res)
    subset_align = SubsetAlignment(msa, res, my_protein, output_folder)
    if subset_align == 0:
        print("  skipped: cant map alignment to domain")
        continue
    # make output name based on input file
    base = os.path.splitext(fname)[0]
    outname = os.path.join(output_folder, f"{gene_name}_subset.fasta")
    # save subset alignment
    SaveSubsetAlignment(subset_align, outname)
    # run iqtree modelfinder
    RunIQtree(outname)
    print(f"  finished: results in {outname}")
