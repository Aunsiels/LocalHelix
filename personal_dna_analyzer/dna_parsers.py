import pandas as pd
from tqdm import tqdm

from personal_dna_analyzer.utils import get_complement


def read_my_heritage(path):
    df = pd.read_csv(path, comment='#')
    res = dict()
    for row in tqdm(df.itertuples(), total=len(df), desc="Loading MyHeritage data"):
        genotype = list(sorted(row.RESULT))
        res[row.RSID] = {"forward": genotype,
                         "backward": get_complement(genotype)}
    return res


def read_23andme(path):
    df = pd.read_csv(path, comment='#', header=0, names=["rsid", "chromosome", "position", "genotype"], sep="\t")
    res = dict()
    for row in tqdm(df.itertuples(), total=len(df), desc="Loading 23andme data"):
        genotype = list(sorted(row.genotype))
        res[row.rsid] = {"forward": genotype,
                         "backward": get_complement(genotype)}
    return res


def read_ancestry(path):
    df = pd.read_csv(path, comment='#', sep="\t")
    res = dict()
    for row in tqdm(df.itertuples(), total=len(df), desc="Loading Ancestry data"):
        genotype = list(sorted([row.allele1, row.allele2]))
        res[row.rsid] = {"forward": genotype,
                         "backward": get_complement(genotype)}
    return res


def get_full_genotype(rs, alleles):
    return "Rs" + rs[2:] + "(" + alleles[0].upper() + ";" + \
        alleles[1].upper() + ")"


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
