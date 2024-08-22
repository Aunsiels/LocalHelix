import pandas as pd
from tqdm import tqdm


def read_my_heritage(path):
    df = pd.read_csv(path, comment='#')
    res = dict()
    for row in tqdm(df.itertuples(), total=len(df), desc="Loading MyHeritage data"):
        res[row.RSID] = {"chromosome": row.CHROMOSOME,
                         "position": row.POSITION,
                         "genotype": list(sorted(row.RESULT))}
    return res


def read_23andme(path):
    df = pd.read_csv(path, comment='#', header=0, names=["rsid", "chromosome", "position", "genotype"], sep="\t")
    res = dict()
    for row in tqdm(df.itertuples(), total=len(df), desc="Loading 23andme data"):
        res[row.rsid] = {"chromosome": row.chromosome,
                         "position": row.position,
                         "genotype": list(sorted(row.genotype))}
    return res


def read_ancestry(path):
    df = pd.read_csv(path, comment='#', sep="\t")
    res = dict()
    for row in tqdm(df.itertuples(), total=len(df), desc="Loading Ancestry data"):
        res[row.rsid] = {"chromosome": row.chromosome,
                         "position": row.position,
                         "genotype": list(sorted([row.allele1, row.allele2]))}
    return res


def auto_load_dna(path):
    with open(path, "r") as f:
        data = f.readline().lower()
    if "myheritage" in data:
        return read_my_heritage(path)
    if "23andme" in data:
        return read_23andme(path)
    if "ancestrydna" in data:
        return read_ancestry(path)
    raise ValueError("DNA parser can only handle MyHeritage, 23andme, and AncestryDNA")
