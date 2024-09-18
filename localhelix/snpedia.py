import json
import os.path
import re

import requests
from bs4 import BeautifulSoup
import pickledb
from tqdm import tqdm
import wikitextparser as wtp

from localhelix.dna_parsers import get_full_genotype

FILE_SNPS = 'snps.json'
FILE_GENOTYPES = "genotypes.json"
FILE_MEDICAL_CONDITIONS = "medical_conditions.json"
DATA_GENOTYPES_HTML = None
DATA_GENOTYPES_WIKITEXT = None
DATA_GENOTYPES = None
USE_WIKITEXT = True

URL_ENDPOINT = "https://bots.snpedia.com/api.php"
REGEX_PMID = re.compile(r"\[PMID (?P<id>\d*)\]")


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


def save_snps(force=False, data_dir="data/"):
    filename = os.path.join(data_dir, FILE_SNPS)
    if not os.path.exists(filename) or force:
        json.dump(get_snp_names(), open(filename, 'w'))


def load_snps(data_dir):
    filename = os.path.join(data_dir, FILE_SNPS)
    return json.load(open(filename))


def save_genotypes(force=False, data_dir="data/"):
    filename = os.path.join(data_dir, FILE_GENOTYPES)
    print(filename)
    if not os.path.exists(filename) or force:
        json.dump(get_genotypes_names(), open(filename, 'w'))


def load_genotypes(data_dir):
    filename = os.path.join(data_dir, FILE_GENOTYPES)
    return json.load(open(filename))


def save_medical_conditions(force=False, data_dir="data/"):
    filename = os.path.join(data_dir, FILE_MEDICAL_CONDITIONS)
    if not os.path.exists(filename) or force:
        json.dump(get_medical_conditions_names(), open(filename, 'w'))


def load_medical_conditions(data_dir):
    filename = os.path.join(data_dir, FILE_MEDICAL_CONDITIONS)
    return json.load(open(filename))


def get_wikitexts(pages):
    if type(pages) is str:
        pages = pages.split("|")
    if len(pages) == 0:
        return dict()
    if len(pages) > 50:
        res = dict()
        i = 0
        for i in range(len(pages) // 50):
            res.update(get_wikitexts(pages[i * 50:(i + 1) * 50]))
        res.update(get_wikitexts(pages[(i + 1) * 50:]))
        return res
    pages = "|".join(pages)
    params = {
        'action': 'query',
        'titles': pages,
        'prop': 'revisions',
        "format": "json",
        "rvprop": "content",
        "rvslots": "main"
    }
    response = requests.get(URL_ENDPOINT, params=params)
    response_json = response.json()
    res = dict()
    for page in response_json["query"]["pages"].values():
        title = page["title"]
        if "revisions" in page:
            wikitext = page["revisions"][0]["slots"]["main"]["*"]
        else:
            wikitext = ""
        res[title] = parse_wikitext(wikitext)
    return res


def parse_wikitext(wikitext):
    parsed = wtp.parse(wikitext)
    text = wikitext
    res = dict()
    for template in parsed.templates:
        t_name = template.name.strip()
        if t_name == "on chip":
            text = text.replace(template.string, "")
        elif t_name == "PMID":
            text = text.replace(template.string,
                                "[<a href=\"https://pubmed.ncbi.nlm.nih.gov/" + template.arguments[
                                    0].value + "\"\>PMID "
                                + template.arguments[0].value + "</a>]")
        elif t_name == "PMID Auto":
            arguments = {x.name.strip(): x.value.strip() for x in template.arguments}
            title = arguments["PMID"]
            if "Title" in arguments:
                title = arguments["Title"]
            text = text.replace(template.string,
                                "[<a href=\"https://pubmed.ncbi.nlm.nih.gov/" + arguments["PMID"] + "\"\>PMID "
                                + arguments["PMID"] + "</a>] " + title)
        else:
            if t_name not in res:
                res[t_name] = []
            res[t_name].append(dict())
            for arg in template.arguments:
                res[t_name][-1][arg.name.strip()] = arg.value.strip()
            text = text.replace(template.string, "")
    for wikilink in parsed.wikilinks:
        text = text.replace(wikilink.string, "<a href=\"https://www.snpedia.com/index.php/" + wikilink.target + "\">" +
                            wikilink.target + "</a>")
    for external_link in parsed.external_links:
        text_url = external_link.text
        if external_link.text is None:
            text_url = external_link.url
        text = text.replace(external_link.string, "<a href=\"" + external_link.url + "\">" + text_url + "</a>")
    for match in REGEX_PMID.finditer(text):
        pmid_id = match.group("id")
        text = text.replace(match.string,
                            "[<a href=\"https://pubmed.ncbi.nlm.nih.gov/" + pmid_id +"\"\>PMID "
                            + pmid_id + "</a>]")
    res["html"] = text.strip().replace("\n", "<br>")
    res["from"] = "wikitext"
    return res


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


def parse_genotype_from_html(genotype):
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
    res["from"] = "html"
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
            DATA_GENOTYPES_HTML.dump()
            ValueError("Unknown orientation " + str(orientation))
        if genotype in genotypes:
            found.add((genotype, forward_genotype))
        else:
            remaining.add(forward_genotype)
        if not exists:
            counter += 1
            if counter % 100 == 0:
                DATA_GENOTYPES_HTML.dump()
    if counter > 0:
        DATA_GENOTYPES_HTML.dump()
    return found, remaining


def get_info_genotype(genotype, autodump=True):
    genotype = genotype.replace("rs", "Rs")
    if DATA_GENOTYPES.exists(genotype):
        return DATA_GENOTYPES.get(genotype), True
    if USE_WIKITEXT:
        info = get_wikitexts(genotype)[genotype]
    else:
        info = parse_genotype_from_html(genotype)
    DATA_GENOTYPES.set(genotype, info)
    if autodump:
        DATA_GENOTYPES.dump()
    return info, False


def get_orientation(genotype, autodump=True):
    genotype = get_root_name(genotype)
    info, exists = get_info_genotype(genotype, autodump)
    if (("from" not in info or info["from"] == "html") and "Orientation" not in info) or \
            ("from" in info and info["from"] == "wikitext" and (
                    "Rsnum" not in info or "Orientation" not in info["Rsnum"][0])):
        return None, exists
    if "from" not in info or info["from"] == "html":
        return info["Orientation"], exists
    else:
        return info["Rsnum"][0]["Orientation"], exists


def get_root_name(genotype):
    if "(" in genotype:
        return genotype.split("(")[0]
    return genotype


def set_snpedia_info(info, res_dict):
    if "from" not in info or info["from"] == "html":
        res_dict["Magnitude"] = str(info.get("Magnitude", "Unknown"))
        res_dict["Repute"] = info.get("Repute", "Unknown")
        res_dict["summary"] = info.get("summary", "")
        res_dict["text"] = info.get("text", "")
        res_dict["was_on_snpedia"] = len(info) != 0
    else:
        genotype = info.get("Genotype", dict())
        if genotype:
            res_dict["Magnitude"] = str(genotype[0].get("magnitude", "Unknown"))
            res_dict["Repute"] = str(genotype[0].get("repute", "Unknown"))
            res_dict["summary"] = str(genotype[0].get("summary", ""))
            res_dict["text"] = info.get("html", "")
            res_dict["was_on_snpedia"] = len(info) != 0
        else:
            res_dict["Magnitude"] = "Unknown"
            res_dict["Repute"] = "Unknown"
            res_dict["summary"] = ""
            res_dict["text"] = info.get("html", "")
            res_dict["was_on_snpedia"] = len(info) != 0


def download_all(data_dir):
    genotypes = load_genotypes(data_dir) + load_snps(data_dir)
    genotypes = [x for x in genotypes if x.startswith("Rs")]
    if USE_WIKITEXT:
        try:
            download_genotypes_wikitext(genotypes)
        except:
            DATA_GENOTYPES_WIKITEXT.dump()
            raise
    else:
        download_genotypes_html(genotypes)


def download_genotypes_html(genotypes):
    counter = 0
    for genotype in tqdm(genotypes):
        if genotype.startswith("Rs") and not DATA_GENOTYPES_HTML.exists(genotype):
            counter += 1
            info = parse_genotype_from_html(genotype)
            DATA_GENOTYPES_HTML.set(genotype, info)
            if counter % 100 == 0:
                print("Saving...")
                DATA_GENOTYPES_HTML.dump()
                print("Saved")
    if counter > 0:
        print("Saving...")
        DATA_GENOTYPES_HTML.dump()
        print("Saved")


def download_genotypes_wikitext(genotypes):
    counter = 0
    genotypes = [genotype for genotype in genotypes if not DATA_GENOTYPES_WIKITEXT.exists(genotype)]
    for counter in tqdm(range(len(genotypes) // 100), total=len(genotypes) // 100,
                        desc="Predownloading relevant pages"):
        for key, value in get_wikitexts(genotypes[counter * 100:(counter + 1) * 100]).items():
            DATA_GENOTYPES_WIKITEXT.set(key, value)
        if counter % 50 == 49:
            try:
                DATA_GENOTYPES_WIKITEXT.dump()
            except KeyboardInterrupt:
                DATA_GENOTYPES_WIKITEXT.dump()
                raise
    for key, value in get_wikitexts(genotypes[(counter + 1) * 100:]).items():
        DATA_GENOTYPES_WIKITEXT.set(key, value)
    try:
        DATA_GENOTYPES_WIKITEXT.dump()
    except KeyboardInterrupt:
        DATA_GENOTYPES_WIKITEXT.dump()
        raise


def get_all_snpedia_entities(dna, genotypes, snps):
    res = set()
    genotypes = set(genotypes)
    snps = set(snps)
    for key, value in dna.items():
        genotype = get_full_genotype(key, value["forward"])
        forward_genotype = genotype
        backward_genotype = get_full_genotype(key, value["backward"])
        if "Rs" + key[2:] in snps:
            res.add("Rs" + key[2:])
        if forward_genotype in genotypes:
            res.add(forward_genotype)
        if backward_genotype in genotypes:
            res.add(backward_genotype)
    return list(res)


def get_all_snpedia_match_genotypes(dna, genotypes, snps):
    if USE_WIKITEXT:
        download_genotypes_wikitext(get_all_snpedia_entities(dna, genotypes, snps))
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


def initialize_snpedia(force=False, data_dir="data/"):
    global DATA_GENOTYPES_HTML, DATA_GENOTYPES_WIKITEXT, USE_WIKITEXT, DATA_GENOTYPES
    save_genotypes(force, data_dir)
    save_snps(force, data_dir)
    save_medical_conditions(force, data_dir)
    DATA_GENOTYPES_HTML = pickledb.load(os.path.join(data_dir, 'data_genotypes.db'), False)
    DATA_GENOTYPES_WIKITEXT = pickledb.load(os.path.join(data_dir, 'data_wikitext_genotypes.db'), False)
    USE_WIKITEXT = True
    if USE_WIKITEXT:
        DATA_GENOTYPES = DATA_GENOTYPES_WIKITEXT
    else:
        DATA_GENOTYPES = DATA_GENOTYPES_HTML


if __name__ == '__main__':
    initialize_snpedia(data_dir="data/")
    download_all("data/")
