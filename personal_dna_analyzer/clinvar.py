import gzip
import json
import os
from collections import namedtuple
from xml.etree.ElementTree import iterparse
import urllib.request

import pandas as pd
from tqdm import tqdm

from personal_dna_analyzer.utils import my_hook, BASES

HAPLOTYPES_TSV = "clinvar_haplotypes.tsv"

CLINVAR_PATHOLOGIES_TSV = "clinvar_pathologies.tsv"

CLINVAR_VARIANTS_TSV = "clinvar_variants.tsv"

CLINVAR_XML = "ClinVarVCVRelease_00-latest.xml.gz"


def get_clinvar_rs_pathologies(pathology_mapping):
    df = pd.read_csv(CLINVAR_PATHOLOGIES_TSV, sep="\t")
    pathologies = dict()
    for row in tqdm(df.itertuples(), total=len(df), desc="Loading ClinVar pathologies"):
        if pd.isna(row.rs_id):
            continue
        id = "rs" + str(int(row.rs_id))
        values = pathology_mapping.get(row.variation_id, set())
        if id not in pathologies:
            pathologies[id] = set()
        pathologies[id].update(values)
    return pathologies


def get_associated_clinvar_pathologies():
    df = pd.read_csv("variant_summary.txt.gz", sep="\t",
                     usecols=["RS# (dbSNP)", "PhenotypeIDS", "PhenotypeList"])
    df["id"] = "rs" + df["RS# (dbSNP)"].astype(str)
    pathologies = dict()
    for row in tqdm(df.itertuples(), total=len(df), desc="Loading ClinVar variant data"):
        id = row.id
        values = get_phenotype_medgen_ids(row)
        if id not in pathologies:
            pathologies[id] = set()
        pathologies[id].update(values)
    return pathologies


def get_variants_mapping():
    df = pd.read_csv("variant_summary.txt.gz", sep="\t",
                     usecols=["RS# (dbSNP)", "Name", "VariationID", "PhenotypeIDS", "PhenotypeList"])
    df["id"] = "rs" + df["RS# (dbSNP)"].astype(str)
    variants = dict()
    for row in tqdm(df.itertuples(), total=len(df), desc="Loading ClinVar variant data"):
        id = row.id
        name = row.Name
        if ">" not in name:
            continue
        idx = name.find(">")
        abnormal = name[idx + 1]
        if id not in variants:
            variants[id] = []
        values = get_phenotype_medgen_ids(row)
        variants[id].append((row.VariationID, abnormal, values))
    return variants


def get_clinvar_variant_pathologies():
    df = pd.read_csv(CLINVAR_PATHOLOGIES_TSV, sep="\t")
    res = dict()
    for row in df.itertuples():
        id = row.variation_id
        if id not in res:
            res[id] = set()
        if not pd.isna(row.name):
            res[id].add(Pathology(row.condition_id, row.name, row.is_pathogenic, row.n_submissions, row.status))
    return res


def get_clinvar_variants(pathology_mapping):
    df = pd.read_csv(CLINVAR_VARIANTS_TSV, sep="\t")
    df = df[df['variation_type'].isin(["single nucleotide variant", "Haplotype"])]
    # All possible types: Complex, CompoundHeterozygote, copy number gain, copy number loss, Deletion, Diplotype,
    # Distinct chromosomes, Duplication, fusion, Haplotype, 'Haplotype, single variant', Indel, Insertion,
    # Inversion, Microsatellite, Phase unknown, protein only, single nucleotide variant, Tandem duplication,
    # Translocation, Variation
    variants = dict()
    for row in tqdm(df.itertuples(), total=len(df), desc="Loading ClinVar variant data"):
        if pd.isna(row.rs_id):
            continue
        id = "rs" + str(int(row.rs_id))
        if not row.allele or pd.isna(row.allele):
            continue
        abnormal = row.allele.strip()
        if id not in variants:
            variants[id] = []
        variants[id].append((row.variation_id, abnormal,
                             pathology_mapping.get(row.variation_id, set())))
    return variants


def get_phenotype_medgen_ids(row):
    ids = row.PhenotypeIDS.split("|")
    ids_filtered = []
    for i in ids:
        final_id = ""
        for pid in i.split(","):
            if pid.startswith("MedGen"):
                final_id = pid
                break
        ids_filtered.append(final_id)
    values = set(Pathology(str(x), y, "", 0, "")
                 for x, y in zip(ids_filtered, row.PhenotypeList.split("|")))
    return values


def get_haplotypes():
    df = pd.read_csv("clinvar_haplotypes.tsv", sep="\t")
    h_to_v = dict()
    v_to_h = dict()
    for row in tqdm(df.itertuples(), total=len(df), desc="Loading haplotypes"):
        haplotype = row.variation_id
        h_to_v[haplotype] = []
        for v in [int(x) for x in row.variants.split(";") if x]:
            if v not in v_to_h:
                v_to_h[v] = []
            v_to_h[v].append(haplotype)
            h_to_v[haplotype].append(v)
    return h_to_v, v_to_h


def find_all_haplotypes(variants, h_to_v, v_to_h):
    possible_haplotypes = set()
    for variant in variants:
        for haplotype in v_to_h.get(variant, []):
            possible_haplotypes.add(haplotype)
    variants = set(variants)
    res = []
    for haplotype in possible_haplotypes:
        if all(x in variants for x in h_to_v.get(haplotype, [])):
            res.append(haplotype)
    return res


def get_clinvar_var_link(text):
    return "Unknown" if text is None \
        else ('<a class=\"link-dark\" href="https://www.ncbi.nlm.nih.gov/clinvar/variation/' +
              str(text) + '">' +
              str(text) + '</a>')


def get_clinvar_variant_from_rs(rs_full, mapping):
    rs = rs_full.split("(")[0].lower()
    alleles = rs_full.split("(")[1][:-1].split(";")
    variants = mapping.get(rs, [])
    res_var = []
    res_pathos = []
    for variant, abnormal, pathos in variants:
        if abnormal and abnormal in alleles:
            res_var.append(variant)
            res_pathos.append(pathos)
    return res_var, res_pathos


def get_tree_xml():
    xml_file = gzip.open('ClinVarVCVRelease_00-latest.xml.gz', 'r')
    root = dict()
    current = root
    previous = []
    for event, row in iterparse(xml_file, events=("start", "end",)):
        if event == "start":
            if row.tag not in current:
                current[row.tag] = dict()
                current[row.tag]["count"] = 0
                current[row.tag]["attrs"] = dict()
                current[row.tag]["text_example"] = ""
            current[row.tag]["count"] += 1
            for k, v in row.attrib.items():
                if k not in current[row.tag]["attrs"]:
                    current[row.tag]["attrs"][k] = v
            if row.text and len(current[row.tag]["text_example"]) < len(row.text.strip()):
                current[row.tag]["text_example"] = row.text.strip()
            previous.append(current)
            current = current[row.tag]
        else:
            current = previous.pop()
        row.clear()
    xml_file.close()
    return root


def print_tree(tree, depth=0):
    for k in tree:
        if k not in ("count", "attrs", "text_example"):
            print(" " * depth * 4 + str(depth) + ":" + k)
            print(" " * depth * 4 + "count = " + str(tree[k]["count"]))
            print(" " * depth * 4 + "attrs = " + ", ".join(str(k) + ": " + str(v) for k, v in tree[k]["attrs"].items()))
            print(" " * depth * 4 + "text example = " + str(tree[k]["text_example"]))
            print_tree(tree[k], depth + 1)


def get_variation(variation):
    variation = str(variation)
    xml_file = gzip.open('ClinVarVCVRelease_00-latest.xml.gz', 'r')
    save = False
    rows = dict()
    current = rows
    previous = []
    for event, row in iterparse(xml_file, events=("start", "end",)):
        if event == "start":
            if row.tag == "VariationArchive" and row.attrib["VariationID"] == variation:
                save = True
                current = rows
        else:
            if row.tag == "VariationArchive" and save:
                break
        if save and event == "start":
            if row.tag not in current:
                current[row.tag] = []
            current[row.tag].append(dict())
            previous.append(current)
            current = current[row.tag][-1]
            current["attrs"] = row.attrib.copy()
            current["text"] = row.text.strip() if row.text else ""
        elif save and event == "end":
            current = previous.pop()
        row.clear()
    xml_file.close()
    with open("example3.json", "w") as f:
        json.dump(rows, f)


def process_clinvar_release():
    xml_file = gzip.open('ClinVarVCVRelease_00-latest.xml.gz', 'r')
    rows = dict()
    current = rows
    previous = []
    variation_file = open(CLINVAR_VARIANTS_TSV, "w")
    variation_file.write("\t".join(("variation_id", "variation_name", "variation_type", "rs_id", "allele")) + "\n")
    pathology_file = open(CLINVAR_PATHOLOGIES_TSV, "w")
    pathology_file.write("\t".join(("variation_id", "rs_id", "name", "condition_id",
                                    "status", "is_pathogenic", "n_submissions")) + "\n")
    haplotype_file = open(HAPLOTYPES_TSV, "w")
    haplotype_file.write("variation_id" + "\t" + "variants" + "\n")
    counter = 0
    for event, row in tqdm(iterparse(xml_file, events=("start", "end",)), total=1077308286):
        counter += 1
        if event == "start" and row.tag == "VariationArchive":
            rows = dict()
            current = rows
        elif event == "end" and row.tag == "VariationArchive":
            try:
                process_variation(rows, variation_file, pathology_file, haplotype_file)
            except KeyError:
                json.dump(rows, open("error.json", "w"))
                raise
        if event == "start":
            if row.tag not in current:
                current[row.tag] = []
            current[row.tag].append(dict())
            previous.append(current)
            current = current[row.tag][-1]
            current["attrs"] = row.attrib.copy()
            current["text"] = row.text.strip() if row.text else ""
        else:
            current = previous.pop()
        row.clear()
    xml_file.close()
    print("Number of entries:", counter)


def process_variation(var_dict, variant_file, pathology_file, haplotype_file):
    base = var_dict["VariationArchive"][0]
    if "ClassifiedRecord" not in base:
        return
    classified_record = base["ClassifiedRecord"][0]
    variation_id = base["attrs"]["VariationID"]
    variation_name = base["attrs"]["VariationName"]
    variation_type = base["attrs"]["VariationType"]
    rs_id = ""
    modification_allele = ""
    if "SimpleAllele" in classified_record:
        simple_allele = classified_record["SimpleAllele"][0]
        if "XRefList" in simple_allele:
            for xref in simple_allele["XRefList"][0]["XRef"]:
                if xref["attrs"]["DB"] == "dbSNP":
                    rs_id = xref["attrs"]["ID"]
        if "Location" in simple_allele:
            location = simple_allele["Location"][0]
            if "SequenceLocation" in location:
                for seq in location["SequenceLocation"]:
                    attrs_location = seq["attrs"]
                    if attrs_location["Assembly"] == "GRCh37":
                        if "alternateAlleleVCF" in attrs_location:
                            modification_allele = attrs_location["alternateAlleleVCF"]
                        elif "Strand" in attrs_location:
                            modification_allele = attrs_location["Strand"]
    variant_file.write("\t".join((variation_id, variation_name, variation_type, rs_id, modification_allele)) + "\n")
    if "Haplotype" in classified_record:
        haplotype = classified_record["Haplotype"][0]
        h_variants = []
        for allele in haplotype["SimpleAllele"]:
            h_variants.append(allele["attrs"]["VariationID"])
        haplotype_file.write(variation_id + "\t" + ";".join(h_variants) + "\n")
    for rcv in classified_record["RCVList"][0]["RCVAccession"]:
        conditions = rcv["ClassifiedConditionList"][0]["ClassifiedCondition"]
        for condition in conditions:
            name = condition["text"]
            if condition["attrs"]:
                condition_id = condition["attrs"]["DB"] + ":" + condition["attrs"]["ID"]
            else:
                condition_id = ""
            if "GermlineClassification" in rcv["RCVClassifications"][0]:
                # We ignore somatic
                classification = rcv["RCVClassifications"][0]["GermlineClassification"][0]
                if len(rcv["RCVClassifications"][0]["GermlineClassification"]) != 1:
                    print("Germline", variation_id)
                    raise KeyError
                status = classification["ReviewStatus"][0]["text"]
                is_pathogenic = classification["Description"][0]["text"]
                n_submissions = classification["Description"][0]["attrs"]["SubmissionCount"]
                pathology_file.write("\t".join((variation_id, rs_id, name, condition_id,
                                                status, is_pathogenic, n_submissions)) + "\n")


def get_clinvar_pathology_link(pathology):
    if pd.isna(pathology.condition_id):
        split = []
    else:
        split = pathology.condition_id.split(":")
    if len(split) == 2:
        link = '<a  class=\"link-dark\" href="https://www.ncbi.nlm.nih.gov/medgen/' + split[1] + '">' + \
               pathology.name + '</a>'
    else:
        link = '<a  class=\"link-dark\" href="https://www.ncbi.nlm.nih.gov/medgen/?term=' + \
               pathology.name.replace(" ", "+") + '">' + pathology.name + '</a>'
    if pathology.is_pathogenic:
        if not pd.isna(pathology.status):
            link += (" [<span title=\"" + pathology.status + "\">" +
                     str(pathology.n_submissions) + "</span>]")
        else:
            link += (" [<span>" +
                     str(pathology.n_submissions) + "</span>]")
    return link


def print_pathologies_html(pathos):
    groups = dict()
    for pathology in pathos:
        is_pathogenic = pathology.is_pathogenic if not pd.isnull(pathology.is_pathogenic) else "Unknown"
        if is_pathogenic not in groups:
            groups[is_pathogenic] = []
        groups[is_pathogenic].append(pathology)
    if not groups:
        return "None"
    res = ["<ul>"]
    for group, pathologies in groups.items():
        res.append("<li>")
        res.append(group)
        res.append(": ")
        pathologies_links = []
        for pathology in pathologies:
            link = get_clinvar_pathology_link(pathology)
            pathologies_links.append(link)
        res.append(", ".join(pathologies_links))
        res.append("</li>")
    res.append("</ul>")
    return "".join(res)


def initialize_clinvar(force=False):
    url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarVCVRelease_00-latest.xml.gz"
    if not os.path.exists(CLINVAR_XML) or force:
        with tqdm(unit='B', unit_scale=True, leave=True, miniters=1,
                  desc="Downloading ClinVar") as t:
            urllib.request.urlretrieve(url, CLINVAR_XML, my_hook(t))
    if not os.path.exists(CLINVAR_VARIANTS_TSV) or force:
        process_clinvar_release()


Pathology = namedtuple("Pathology", ["condition_id", "name", "is_pathogenic",
                                     "n_submissions", "status"])

if __name__ == '__main__':
    get_variation(17525)
    # process_clinvar_release()
