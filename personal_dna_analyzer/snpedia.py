import json
import os.path

import requests
from bs4 import BeautifulSoup
import pickledb
from tqdm import tqdm

from personal_dna_analyzer.dna_parsers import get_full_genotype

FILE_SNPS = 'snps.json'
FILE_GENOTYPES = "genotypes.json"
FILE_MEDICAL_CONDITIONS = "medical_conditions.json"
DATA_GENOTYPES = pickledb.load('data_genotypes.db', False)

URL_ENDPOINT = "https://bots.snpedia.com/api.php"


def get_all_category(category_name):
    params = {
        'action': 'query',
        'format': 'json',
        'list': "categorymembers",
        "cmtitle": "Category:" + category_name,
        "cmlimit": 500
    }
    response = requests.get(URL_ENDPOINT, params=params)
    data = response.json()
    members = [x["title"] for x in data['query']['categorymembers']]
    pbar = tqdm(desc="Downloading from SNPedia " + category_name)
    while "continue" in data:
        next_page = data["continue"]["cmcontinue"]
        params = {
            'action': 'query',
            'format': 'json',
            'list': "categorymembers",
            "cmtitle": "Category:" + category_name,
            "cmlimit": 500,
            "cmcontinue": next_page
        }
        response = requests.get(URL_ENDPOINT, params=params)
        data = response.json()
        members += [x["title"] for x in data['query']['categorymembers']]
        pbar.update(1)
    pbar.close()
    return members


def get_snp_names():
    return get_all_category("Is_a_snp")


def get_genotypes_names():
    return get_all_category("Is_a_genotype")


def get_medical_conditions_names():
    return get_all_category("Is_a_medical_condition")


def save_snps(force=False):
    if not os.path.exists(FILE_SNPS) or force:
        json.dump(get_snp_names(), open(FILE_SNPS, 'w'))


def load_snps():
    return json.load(open(FILE_SNPS))


def save_genotypes(force=False):
    if not os.path.exists(FILE_GENOTYPES) or force:
        json.dump(get_genotypes_names(), open(FILE_GENOTYPES, 'w'))


def load_genotypes():
    return json.load(open(FILE_GENOTYPES))


def save_medical_conditions(force=False):
    if not os.path.exists(FILE_MEDICAL_CONDITIONS) or force:
        json.dump(get_medical_conditions_names(), open(FILE_MEDICAL_CONDITIONS, 'w'))


def load_medical_conditions():
    return json.load(open(FILE_MEDICAL_CONDITIONS))


def get_page(page):
    params = {
        'action': 'parse',
        'format': 'json',
        'page': page,
        'prop': 'text',
        'redirects': ''
    }
    response = requests.get(URL_ENDPOINT, params=params)
    try:
        response_json = response.json()
    except requests.exceptions.JSONDecodeError:
        print("Strange", page)
        print(response.content.decode("utf-8"))
        response_json = json.loads(response.content.decode("utf-8"))
    raw_html = response_json['parse']['text']['*']
    return raw_html


def parse_genotype(genotype):
    raw_html = get_page(genotype)
    soup = BeautifulSoup(raw_html, 'html.parser')
    rows = soup.find_all('tr')
    res = {"others": [], "summary": ""}
    for row in rows:
        columns = row.find_all('td')
        if len(columns) == 0:
            pass
        elif len(columns) == 1:
            res["summary"] = columns[0].text.strip()
        elif len(columns) == 2:
            if columns[0].text.strip() != "mentioned":
                res[columns[0].text.strip()] = columns[1].text.strip()
        elif len(columns) == 3:
            res["others"].append((columns[0].text.strip(), columns[1].text.strip(), columns[2].text.strip()))
        else:
            print(columns)
    res["text"] = [str(x) for x in soup.find_all("p") if len(x.text.strip()) > 0]
    return res


def separate_snpedia_variants(dna, genotypes, snps):
    found, remaining = set(), set()
    genotypes = set(genotypes)
    snps = set(snps)
    counter = 0
    for key, value in tqdm(dna.items(), total=len(dna), desc="Separating SNPedia variants"):
        genotype = get_full_genotype(key, value["forward"])
        forward_genotype = genotype
        backward_genotype = get_full_genotype(key, value["backward"])
        if "Rs" + key[2:] not in snps or \
                (forward_genotype not in genotypes and backward_genotype not in genotypes):
            remaining.add(forward_genotype)
            continue
        orientation, exists = get_orientation(key, autodump=False)
        if orientation is None:
            remaining.add(forward_genotype)
            continue
        if orientation == "plus":
            genotype = forward_genotype
        elif orientation == "minus":
            genotype = backward_genotype
        else:
            DATA_GENOTYPES.dump()
            ValueError("Unknown orientation " + str(orientation))
        if genotype in genotypes:
            found.add((genotype, forward_genotype))
        else:
            remaining.add(forward_genotype)
        if not exists:
            counter += 1
            if counter % 100 == 0:
                DATA_GENOTYPES.dump()
    if counter > 0:
        DATA_GENOTYPES.dump()
    return found, remaining


def get_info_genotype(genotype, autodump=True):
    if DATA_GENOTYPES.exists(genotype):
        return DATA_GENOTYPES.get(genotype), True
    info = parse_genotype(genotype)
    DATA_GENOTYPES.set(genotype, info)
    if autodump:
        DATA_GENOTYPES.dump()
    return info, False


def get_orientation(genotype, autodump=True):
    if "(" in genotype:
        genotype = genotype.split("(")[0]
    info, exists = get_info_genotype(genotype, autodump)
    if "Orientation" not in info:
        return None, exists
    return info["Orientation"], exists


def download_all():
    genotypes = load_genotypes()
    genotypes = [x for x in genotypes if x.startswith("Rs")]
    counter = 0
    for genotype in tqdm(genotypes):
        if genotype.startswith("Rs") and not DATA_GENOTYPES.exists(genotype):
            counter += 1
            info = parse_genotype(genotype)
            DATA_GENOTYPES.set(genotype, info)
            if counter % 100 == 0:
                print("Saving...")
                DATA_GENOTYPES.dump()
                print("Saved")
    if counter > 0:
        print("Saving...")
        DATA_GENOTYPES.dump()
        print("Saved")


def get_all_snpedia_match_genotypes(dna, genotypes, snps):
    known_genotypes, remaining = separate_snpedia_variants(dna, genotypes, snps)
    res = {}
    counter = 0
    for genotype in tqdm(known_genotypes, desc="Gathering SNPedia information"):
        res[genotype[1]], exists = get_info_genotype(genotype[0], autodump=False)
        if not exists:
            counter += 1
            if counter % 100 == 0:
                DATA_GENOTYPES.dump()
    if counter > 0:
        DATA_GENOTYPES.dump()
    for r in remaining:
        res[r] = dict()
    return res


def get_snpedia_link(text):
    if text and text.lower().startswith("rs"):
        return "<a  class=\"link-dark\" href=\"https://www.snpedia.com/index.php/" + text + \
            "\">" + text + "</a>"
    return "None"


def initialize_snpedia(force=False):
    save_genotypes(force)
    save_snps(force)
    save_medical_conditions(force)


if __name__ == '__main__':
    download_all()
