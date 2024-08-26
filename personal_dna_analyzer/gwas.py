import os.path
import re

import pandas as pd
from tqdm import tqdm
import urllib.request

from personal_dna_analyzer.utils import my_hook

GWAS_ASSOCIATIONS_FILENAME = "gwas_associations.tsv"


def get_gwas_html(rs):
    gwas = rs["GWAS"]
    if not gwas:
        return "<b>GWAS Traits</b>: None"
    res = "<b>GWAS Traits</b>: <br><ul>"
    for al, value in gwas.items():
        res += "<li>" + al + ": "
        temp = []
        for key, count_or_beta in value.items():
            count, or_beta = count_or_beta
            mapped_trait, trait_uri, text = key
            text = (text + ", ") if text else text
            if text:
                or_beta_text = "Reported Betas:"
            else:
                or_beta_text = "Reported <span text=\">1 accentuate the trait, <1 reduces the trait\">Odd Ratio</span>:"
            or_beta_text += " [" + ", ".join(str(x) for x in or_beta) + "]"
            temp.append('<a  class=\"link-dark\" href="' + trait_uri + '">' + mapped_trait + "</a> (" + text +
                        str(count) + " study(s), " + or_beta_text + ")")
        res += ", ".join(temp)
        res += "<br>"
    res += "</ul><br><b>GWAS page</b>: <a class=\"link-dark\" href=\"https://www.ebi.ac.uk/gwas/variants/" + \
           rs["rs"].split("(")[0] + \
           "\">" + rs["rs"].split("(")[0] + "</a>"
    return res


def get_gwas_traits():
    df = pd.read_csv("gwas_associations.tsv", sep="\t")
    df.rename(columns={x: x.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "").replace("-", "_")
              .replace("95%", "p").replace("[", "").replace("]", "")
                       for x in df.columns}, inplace=True)
    df.fillna("", inplace=True)
    res = dict()
    for row in tqdm(df.itertuples(), desc="Loading GWAS traits", total=len(df)):
        # disease = getattr(row, "DISEASE_TRAIT")
        snp = getattr(row, "SNPS")
        if "-" in getattr(row, "STRONGEST_SNP_RISK_ALLELE"):
            abnormal = getattr(row, "STRONGEST_SNP_RISK_ALLELE").split("-")[1]
        else:
            abnormal = ""
        text = getattr(row, "p_CI_TEXT")
        mapped_trait = getattr(row, "MAPPED_TRAIT")
        trait_uri = getattr(row, "MAPPED_TRAIT_URI")
        or_or_beta = getattr(row, "OR_or_BETA")
        key = (snp, abnormal)
        text = re.sub(r"\[.*]", "", text).strip()
        if key not in res:
            res[key] = dict()
        value = (mapped_trait, trait_uri, text)
        if value not in res[key]:
            res[key][value] = [0, []]
        res[key][value][0] += 1
        res[key][value][1].append(or_or_beta)
    return res


def initialize_gwas(force=False):
    if not os.path.exists(GWAS_ASSOCIATIONS_FILENAME) or force:
        with tqdm(unit='B', unit_scale=True, leave=True, miniters=1,
                  desc="Downloading GWAS") as t:
            urllib.request.urlretrieve("https://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
                                       GWAS_ASSOCIATIONS_FILENAME, my_hook(t))
