import json

import requests
from bs4 import BeautifulSoup
import pickledb
from tqdm import tqdm

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
    return members


def get_snp_names():
    return get_all_category("Is_a_snp")


def get_genotypes_names():
    return get_all_category("Is_a_genotype")


def get_medical_conditions_names():
    return get_all_category("Is_a_medical_condition")


def save_snps():
    json.dump(get_snp_names(), open(FILE_SNPS, 'w'))


def load_snps():
    return json.load(open(FILE_SNPS))


def save_genotypes():
    json.dump(get_genotypes_names(), open(FILE_GENOTYPES, 'w'))


def load_genotypes():
    return json.load(open(FILE_GENOTYPES))


def save_medical_conditions():
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
        print("Bizarre")
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


def separate_snpedia_variants(dna, genotypes):
    my_genotypes = set(key[0].upper() + key[1:].lower() + "(" + value["genotype"][0].upper() + ";" +
                       value["genotype"][1].upper() + ")"
                       for key, value in dna.items())
    intersection = set(genotypes).intersection(my_genotypes)
    return intersection, my_genotypes.difference(intersection)


def get_info_genotype(genotype, autodump=True):
    if DATA_GENOTYPES.exists(genotype):
        return DATA_GENOTYPES.get(genotype), True
    info = parse_genotype(genotype)
    DATA_GENOTYPES.set(genotype, info)
    if autodump:
        DATA_GENOTYPES.dump()
    return info, False


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


def get_all_snpedia_match_genotypes(dna, genotypes):
    known_genotypes, remaining = separate_snpedia_variants(dna, genotypes)
    res = {}
    counter = 0
    for genotype in known_genotypes:
        res[genotype], exists = get_info_genotype(genotype, autodump=False)
        if exists:
            counter += 1
        if counter % 100 == 0:
            DATA_GENOTYPES.dump()
    if counter > 0:
        DATA_GENOTYPES.dump()
    for r in remaining:
        res[r] = dict()
    return res


if __name__ == '__main__':
    # process_clinvar_release()
    # main()
    download_all()
