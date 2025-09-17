import os, csv, re

def ParseIQtree(path):
    text = open(path).read()
    model  = re.search(r'Model of substitution:\s*(\S+)', text)
    pinv   = re.search(r'Proportion of invariable sites:\s*([\d.]+)', text)
    alpha  = re.search(r'Gamma shape alpha:\s*([\d.]+)', text)
    rates  = re.search(r'Site proportion and rates:\s*(.+)', text)

    model_str = model.group(1) if model else ""
    I_str     = pinv.group(1) if pinv else ""
    alpha_str = alpha.group(1) if alpha else ""
    rates_str = rates.group(1) if rates else ""

    # base = model string minus trailing +I/+Gx/+Rx
    toks = model_str.split("+")
    keep = [t for t in toks if not re.match(r'^(I|G\d+|R\d+)$', t)]
    base = "+".join(keep)

    # --- build AliSim-ready string ---
    parts = [base]
    if I_str:
        parts.append(f"I{{{I_str}}}")
    if alpha_str:
        m2 = re.search(r'G(\d+)', model_str)
        if m2:
            k = m2.group(1)
            parts.append(f"G{k}{{{alpha_str}}}")
    if rates_str:
        pairs = re.findall(r'\(([\d.]+),([\d.]+)\)', rates_str)
        k = len(pairs)
        flat = ",".join([",".join(p) for p in pairs])
        parts.append(f"R{k}" + "{" + flat + "}")

    alisim_model = "+".join([p for p in parts if p])

    return {
        "file": os.path.basename(path).rsplit("_", 1)[0],
        "model": model_str,
        "I": I_str,
        "alpha": alpha_str,
        "rates": rates_str,
        "alisim_model": alisim_model
    }


folder = "/local/home/zool2506/Simulations/python_sim/domain_param_estimation/Subset_MSA"
rows = [ParseIQtree(os.path.join(folder,f)) for f in os.listdir(folder) if f.endswith(".iqtree")]
out_file = "/local/home/zool2506/Simulations/python_sim/domain_param_estimation/domain_models.csv"

with open(out_file,"w",newline="") as out:
    writer = csv.DictWriter(out, fieldnames=["file","model","I","alpha","rates", "alisim_model"])
    writer.writeheader()
    writer.writerows(rows)