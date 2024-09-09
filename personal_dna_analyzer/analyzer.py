import argparse
from collections import Counter

import pandas as pd
from tqdm import tqdm

from personal_dna_analyzer.dna_parsers import auto_load_dna
from personal_dna_analyzer.gwas import get_gwas_html, get_gwas_traits, initialize_gwas
from personal_dna_analyzer.snpedia import get_all_snpedia_match_genotypes, load_genotypes, initialize_snpedia, \
    get_snpedia_link, load_snps, set_snpedia_info
from personal_dna_analyzer.clinvar import get_clinvar_variant_from_rs, get_clinvar_variant_pathologies, \
    get_clinvar_rs_pathologies, get_clinvar_variants, print_pathologies_html, initialize_clinvar, get_clinvar_var_link, \
    get_haplotypes, find_all_haplotypes


def get_summaries_dict(dna, genotypes, snps, pathologies, variants_mapping, gwas_traits, haplotypes, var_to_pathology):
    genotype_infos = get_all_snpedia_match_genotypes(dna, genotypes, snps)
    res = []
    all_clinvar_variants = []
    for key, value in tqdm(sorted(genotype_infos.items(),
                                  key=lambda x: -(
                                          float(x[1].get("Magnitude", 1)) + len(x[1].get("summary", "")) / 1000.0)),
                           total=len(genotype_infos), desc="Generating report"):
        base_rs = key.split("(")[0].lower()
        al1, al2 = key.split("(")[1].replace(")", "").split(";")
        variant_list, pathos_list = get_clinvar_variant_from_rs(key, variants_mapping)
        for variant in variant_list:
            all_clinvar_variants.append(variant)
        all_pathologies = list(sorted(pathologies.get(base_rs, []), key=lambda x: x.name.lower()))
        variant_pathologies_list = []
        variant_pathology_names = set()
        for pathos in pathos_list:
            variant_pathologies_list.append(list(sorted(pathos, key=lambda x: x.name.lower())))
            variant_pathology_names.update({x.name for x in variant_pathologies_list[-1]})
        all_pathologies = [x for x in all_pathologies if x.name not in variant_pathology_names]
        if not value and not all_pathologies and all(len(x) == 0 for x in variant_pathologies_list) and \
                (base_rs, al1) not in gwas_traits and (base_rs, al2) not in gwas_traits:
            continue
        temp = dict()
        temp["rs"] = key
        set_snpedia_info(value, temp)
        temp["ClinVarAllPathologies"] = all_pathologies
        temp["ClinVarVariants"] = variant_list
        temp["ClinVarVariantPathologies"] = variant_pathologies_list
        temp["GWAS"] = dict()
        if (base_rs, al1) in gwas_traits:
            temp["GWAS"][al1] = gwas_traits[(base_rs, al1)]
        if (base_rs, al2) in gwas_traits:
            temp["GWAS"][al2] = gwas_traits[(base_rs, al2)]
        res.append(temp)
    h_to_v, v_to_h = haplotypes
    found_haplotypes = find_all_haplotypes(all_clinvar_variants, h_to_v, v_to_h)
    for found_haplotype in found_haplotypes:
        variant_pathologies = list(sorted(var_to_pathology.get(found_haplotype, []), key=lambda x: x.name.lower()))
        temp = dict()
        temp["Magnitude"] = "Unknown"
        temp["Repute"] = "Unknown"
        temp["summary"] = ""
        temp["rs"] = "Haplotype: " + str(found_haplotype) + " = " + " + ".join(str(x) for x in h_to_v[found_haplotype])
        temp["text"] = ""
        temp["was_on_snpedia"] = False
        temp["ClinVarAllPathologies"] = []
        temp["ClinVarVariants"] = [found_haplotype]
        temp["ClinVarVariantPathologies"] = [variant_pathologies]
        temp["GWAS"] = dict()
        res.append(temp)
    res = sorted(res, key=lambda x: -score_summary_entry(x))
    return res


def summaries_results_in_html(all_rs):
    res = ["<div>",
           "<h2>Summary of the results</h2>",
           "<b>Do not use for medical advice and always consult your doctor.</b>",
           "<p>We found matches in our databases for " + str(len(all_rs)) + " variants.</p>"]
    snpedia_counter = 0
    pathologies_all = []
    pathologies_all_pathogenic = []
    pathologies_variants = []
    pathologies_variants_pathogenic = []
    gwas_traits = []
    for rs in all_rs:
        if rs["was_on_snpedia"]:
            snpedia_counter += 1
        for pathology_list in rs["ClinVarVariantPathologies"]:
            for pathology in pathology_list:
                _add_pathology_to_counters(pathologies_variants, pathologies_variants_pathogenic, pathology)
        for pathology in rs["ClinVarAllPathologies"]:
            _add_pathology_to_counters(pathologies_all, pathologies_all_pathogenic, pathology)
        gwas = rs["GWAS"]
        for al, value in gwas.items():
            for key, count_or_beta in value.items():
                _, or_betas = count_or_beta
                mapped_trait, _, text = key
                if not text:
                    for or_beta in or_betas:
                        if or_beta:
                            if float(or_beta) > 1.0:
                                text = "Odd ratio >1"
                            else:
                                text = "Odd ratio <1"
                gwas_traits.append(mapped_trait + " (" + text + ")")
    pathologies_all = Counter(pathologies_all)
    pathologies_all_pathogenic = Counter(pathologies_all_pathogenic)
    pathologies_variants = Counter(pathologies_variants)
    pathologies_variants_pathogenic = Counter(pathologies_variants_pathogenic)
    gwas_traits = Counter(gwas_traits)
    res.append("<p>Among theses results, " + str(snpedia_counter) + " variants were on SNPedia.</p>")
    res.append("<p>On ClinVar, we found " + str(len(pathologies_variants)) + " different pathologies or traits " +
               "associated with " +
               "your variants. Among them, " + str(len(pathologies_variants_pathogenic)) + " were classified as " +
               "pathogenic or likely pathogenic. Here are the most frequent pathogenic traits:</p>")
    add_most_commons_list(pathologies_variants_pathogenic, res)
    res.append("<p>On GWAS, we found " + str(len(gwas_traits)) + " different pathologies or traits " +
               "associated with your variants. Here are the most frequent traits:</p>")
    add_most_commons_list(gwas_traits, res)
    res.append("</div>")
    return "".join(res)


def _add_pathology_to_counters(pathologies_variants, pathologies_variants_pathogenic, pathology):
    if pathology.name and not pd.isna(pathology.name) and \
            pathology.name not in ["not provided", "not specified"]:
        pathologies_variants.append(pathology.name)
        if not pd.isna(pathology.is_pathogenic) and "pathogenic" in pathology.is_pathogenic.lower():
            pathologies_variants_pathogenic.append(pathology.name)


def add_most_commons_list(pathologies_variants_pathogenic, res):
    res.append("<ul>")
    for key, value in pathologies_variants_pathogenic.most_common(10):
        res.append("<li>")
        res.append(key + " (" + str(value) + ")")
        res.append("</li>")
    res.append("</ul>")


def get_card_header(rs):
    if rs["Repute"] == "Bad":
        content = '<div class="card text-white bg-danger">'
    elif rs["Repute"] == "Good":
        content = '<div class="card text-white bg-success">'
    else:
        content = '<div class="card bg-light">'
    return content


def rs_to_html(rs):
    clinvar_all_pathologies = print_pathologies_html(rs["ClinVarAllPathologies"])
    clinvar_variant_pathologies = [print_pathologies_html(x) for x in rs["ClinVarVariantPathologies"]]
    var_links = [get_clinvar_var_link(x) for x in rs["ClinVarVariants"]]
    content = get_card_header(rs)
    if type(rs["text"]) is str:
        rs_text = [rs["text"]]
    else:
        rs_text = rs["text"]
    snpedia_text = [x.replace("href=\"/index.php",
                              "href=\"https://www.snpedia.com/index.php")
                    .replace("href=\"//www.ncbi.nlm.nih.gov/pubmed",
                             "href=\"https://www.ncbi.nlm.nih.gov/pubmed")
                    .replace("href=",
                             "class=\"link-dark\" href=")
                    for x in rs_text]
    content += """
      <h5 class="card-header">""" + rs["rs"] + """</h5>
      <div class="card-body">
        <h5 class="card-title">""" + rs["summary"] + """</h5>
        <p class="card-text">""" + " ".join(snpedia_text) + """
        <br> <span title=\"Importance annotated by SNPedia community\"><b>Magnitude</b></span>: """ + rs["Magnitude"] + \
               ", <span title=\"Good or bad SNP, annotaed by SNPedia community\"><b>Repute</b></span>: " + \
               rs["Repute"] + \
               "<br>" + \
               "<b>SNPedia Variant</b>: " + get_snpedia_link(rs["rs"]) + "<br>" + \
               "<b>SNPedia Base SNP</b>: " + get_snpedia_link(rs["rs"].split("(")[0]) + "<br>"
    for var_link, clinvar_var_pathos in zip(var_links, clinvar_variant_pathologies):
        content += ("<span title=\"ClinVar page of your variant\"><b>ClinVar "
                    "Variant</b></span>: ") + \
                   var_link + "<br>" + \
                   ("<span title=\"Associated traits/pathologies with your variants on ClinVar\"><b>Variants "
                    "pathologies for " + var_link + " (ClinVar)</b></span>: ") + \
                   clinvar_var_pathos + "<br>"
    content += ("<span title=\"Other pathologies associated with this SNP. "
                "Interesting to know if you have a good variant\"><b>Pathologies/traits you avoided (ClinVar)"
                "</b></span>: ") + clinvar_all_pathologies + "<br>" + get_gwas_html(rs) + "<br>" + """</p>
      </div>
    </div>
    """
    return content


def get_html_page(all_rss):
    cards = [rs_to_html(rs) for rs in all_rss]
    content = summaries_results_in_html(all_rss) + "<h2>All variants</h2>" + "<br>".join(cards)
    html = """
    <!DOCTYPE html>
        <html lang="en">
        <head>
          <title>DNA Analyser</title>
          <meta charset="utf-8">
          <meta name="viewport" content="width=device-width, initial-scale=1">
          <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
          <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js"></script>
        </head>
        <body>
        
        <div class="container-fluid p-5 bg-primary text-white text-center"> <h1>DNA Analyser</h1> <p>Your results are 
        ready!</p> <form action="https://www.paypal.com/donate" method="post" target="_top"> <input type="hidden" 
        name="hosted_button_id" value="GR3D64Y7S7TU2" /> <input type="image" 
        src="https://www.paypalobjects.com/en_US/i/btn/btn_donate_LG.gif" border="0" name="submit" title="PayPal - 
        The safer, easier way to pay online!" alt="Donate with PayPal button" /> <img alt="" border="0" 
        src="https://www.paypal.com/en_FR/i/scr/pixel.gif" width="1" height="1" /> </form>

        </div>
          
        <div class="container mt-5">
        """ + content + """  
        </div>
        
        </body>
        </html>
    """
    return html


def score_summary_entry(entry):
    magnitude = 0
    if entry["Magnitude"] != "Unknown" and entry["Magnitude"] != "":
        magnitude = float(entry["Magnitude"])
    n_gwas = 0
    for value in entry["GWAS"]:
        n_gwas += len(value)
    n_pathologies = len(entry["ClinVarAllPathologies"])
    n_variant_pathologies = sum(len(x) for x in entry["ClinVarVariantPathologies"])
    n_pathogenic = len([x for y in entry["ClinVarVariantPathologies"] for x in y if not pd.isna(x.is_pathogenic) and
                        "pathogenic" in x.is_pathogenic.lower() and "Conflicting" not in x.is_pathogenic])
    return float(magnitude) / 2.0 + n_gwas / 20.0 + min(n_pathologies, 30) / 60.0 + \
        n_pathogenic / 5.0 + min((n_variant_pathologies - n_pathogenic), 10) / 20.0


def initialize_all(force=False, data_dir="data/"):
    initialize_gwas(force, data_dir)
    initialize_snpedia(force, data_dir)
    initialize_clinvar(force, data_dir)


def main(input_filename, output_filename, force_reload=False, data_dir="data/"):
    initialize_all(force=force_reload, data_dir=data_dir)
    dna = auto_load_dna(input_filename)
    genotypes = load_genotypes(data_dir)
    snps = load_snps(data_dir)
    gwas_traits = get_gwas_traits(data_dir)
    pathology_mapping = get_clinvar_variant_pathologies(data_dir)
    rs_pathologies = get_clinvar_rs_pathologies(pathology_mapping, data_dir)
    variants_mapping = get_clinvar_variants(pathology_mapping, data_dir)
    haplotypes = get_haplotypes(data_dir)
    all_rss = get_summaries_dict(dna, genotypes, snps, rs_pathologies, variants_mapping, gwas_traits, haplotypes,
                                 pathology_mapping)
    html = get_html_page(all_rss)
    with open(output_filename, "w") as f:
        f.write(html)
    print("Report written to " + output_filename)


def get_arguments():
    parser = argparse.ArgumentParser(
        prog="Personal DNA Analyzer",
        description="This program analyzes your DNA to find interesting insights. Do not use for medical advice and"
                    "always consult your doctor."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input file. Accepted formats: MyHeritage, 23andMe, AncestryDNA.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output html file.")
    parser.add_argument("-f", "--force_reload", action="store_true",
                        help="Force reload all data sources (time consuming)")
    parser.add_argument("-d", "--data_dir", default="data/",
                        help="Directory where data files are located.")
    return parser.parse_args()


if __name__ == '__main__':
    args = get_arguments()
    main(args.input, args.output, args.force_reload, args.data_dir)
