import argparse
from collections import Counter

import pandas as pd
from tqdm import tqdm

from personal_dna_analyzer.dna_parsers import auto_load_dna
from personal_dna_analyzer.gwas import get_gwas_html, get_gwas_traits
from personal_dna_analyzer.snpedia import get_all_snpedia_match_genotypes, load_genotypes
from personal_dna_analyzer.clinvar import get_clinvar_variant_from_rs, get_clinvar_variant_pathologies, \
    get_clinvar_rs_pathologies, get_clinvar_variants, print_pathologies_html


def get_summaries_dict(dna, genotypes, pathologies, variants_mapping, gwas_traits):
    genotype_infos = get_all_snpedia_match_genotypes(dna, genotypes)
    res = []
    for key, value in tqdm(sorted(genotype_infos.items(),
                                  key=lambda x: -(
                                          float(x[1].get("Magnitude", 1)) + len(x[1].get("summary", "")) / 1000.0)),
                           total=len(genotype_infos), desc="Generating report"):
        base_rs = key.split("(")[0].lower()
        al1, al2 = key.split("(")[1].replace(")", "").split(";")
        variant, pathos = get_clinvar_variant_from_rs(key, variants_mapping)
        all_pathologies = list(sorted(pathologies.get(base_rs, []), key=lambda x: x.name.lower()))
        variant_pathologies = list(sorted(pathos, key=lambda x: x.name.lower()))
        if not value and not all_pathologies and not variant_pathologies and \
                (base_rs, al1) not in gwas_traits and (base_rs, al2) not in gwas_traits:
            continue
        temp = dict()
        temp["Magnitude"] = str(value.get("Magnitude", "Unknown"))
        temp["Repute"] = value.get("Repute", "Unknown")
        temp["summary"] = value.get("summary", "")
        temp["rs"] = key
        temp["text"] = value.get("text", "")
        temp["was_on_snpedia"] = len(value) != 0
        temp["ClinVarAllPathologies"] = all_pathologies
        temp["ClinVarVariant"] = variant
        temp["ClinVarVariantPathologies"] = variant_pathologies
        temp["GWAS"] = dict()
        if (base_rs, al1) in gwas_traits:
            temp["GWAS"][al1] = gwas_traits[(base_rs, al1)]
        if (base_rs, al2) in gwas_traits:
            temp["GWAS"][al2] = gwas_traits[(base_rs, al2)]
        res.append(temp)
    res = sorted(res, key=lambda x: -score_summary_entry(x))
    return res


def summaries_results_in_html(all_rs):
    res = ["<div>",
           "<h2>Summary of the results</h2>",
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
        for pathology in rs["ClinVarVariantPathologies"]:
            if pathology.name and not pd.isna(pathology.name) and \
                    pathology.name not in ["not provided", "not specified"]:
                pathologies_variants.append(pathology.name)
                if not pd.isna(pathology.is_pathogenic) and "pathogenic" in pathology.is_pathogenic.lower():
                    pathologies_variants_pathogenic.append(pathology.name)
        for pathology in rs["ClinVarAllPathologies"]:
            if pathology.name and not pd.isna(pathology.name) and \
                    pathology.name not in ["not provided", "not specified"]:
                pathologies_all.append(pathology.name)
                if not pd.isna(pathology.is_pathogenic) and "pathogenic" in pathology.is_pathogenic.lower():
                    pathologies_all_pathogenic.append(pathology.name)
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
    res.append("<ul>")
    for key, value in pathologies_variants_pathogenic.most_common(10):
        res.append("<li>")
        res.append(key + " (" + str(value) + ")")
        res.append("</li>")
    res.append("</ul>")
    res.append("<p>On GWAS, we found " + str(len(gwas_traits)) + " different pathologies or traits " +
               "associated with your variants. Here are the most frequent traits:</p>")
    res.append("<ul>")
    for key, value in gwas_traits.most_common(10):
        res.append("<li>")
        res.append(key + " (" + str(value) + ")")
        res.append("</li>")
    res.append("</ul>")
    res.append("</div>")
    return "".join(res)


def rs_to_html(rs):
    clinvar_all_pathologies = print_pathologies_html(rs["ClinVarAllPathologies"])
    clinvar_variant_pathologies = print_pathologies_html(rs["ClinVarVariantPathologies"])
    var_link = "Unknown" if rs["ClinVarVariant"] is None \
        else ('<a class=\"link-dark\" href="https://www.ncbi.nlm.nih.gov/clinvar/variation/' +
              str(rs["ClinVarVariant"]) + '">' +
              str(rs["ClinVarVariant"]) + '</a>')
    if rs["Repute"] == "Bad":
        content = '<div class="card text-white bg-danger">'
    elif rs["Repute"] == "Good":
        content = '<div class="card text-white bg-success">'
    else:
        content = '<div class="card bg-light">'
    content += """
      <h5 class="card-header">""" + rs["rs"] + """</h5>
      <div class="card-body">
        <h5 class="card-title">""" + rs["summary"] + """</h5>
        <p class="card-text">""" + "<br>".join(rs["text"]) + """
        <br> <span title=\"Importance annotated by SNPedia community\"><b>Magnitude</b></span>: """ + rs["Magnitude"] + \
               ", <span title=\"Good or bad SNP, annotaed by SNPedia community\"><b>Repute</b></span>: " + \
               rs["Repute"] + \
               "<br>" + \
               "<b>SNPedia Variant</b>: <a  class=\"link-dark\" href=\"https://www.snpedia.com/index.php/" + rs["rs"] + \
               "\">" + \
               rs["rs"] + "</a><br>" + \
               "<b>SNPedia Base SNP</b>: <a  class=\"link-dark\" href=\"https://www.snpedia.com/index.php/" + \
               rs["rs"].split("(")[0] + "\">" + \
               rs["rs"].split("(")[0] + "</a><br>" + \
               ("<span title=\"ClinVar page of your variant\"><b>ClinVar "
                "Variant</b></span>: ") + \
               var_link + "<br>" + \
               ("<span title=\"Associated traits/pathologies with your variants on ClinVar\"><b>Your variants "
                "pathologies (ClinVar)</b></span>: ") + \
               clinvar_variant_pathologies + "<br>" + \
               ("<span title=\"All pathologies associated with this SNP, not necessary what you have but what others "
                "have. Interesting to know if you have a good variant\"><b>Other possible variant pathologies, "
                "not necessary yours (ClinVar)</b></span>: ") + \
               clinvar_all_pathologies + "<br>" + \
               get_gwas_html(rs) + "<br>" \
               """</p>
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
    if entry["Magnitude"] != "Unknown":
        magnitude = entry["Magnitude"]
    n_gwas = 0
    for value in entry["GWAS"]:
        n_gwas += len(value)
    n_pathologies = len(entry["ClinVarAllPathologies"])
    n_variant_pathologies = len(entry["ClinVarVariantPathologies"])
    n_pathogenic = len([x for x in entry["ClinVarVariantPathologies"] if not pd.isna(x.is_pathogenic) and
                        "pathogenic" in x.is_pathogenic.lower()])
    return float(magnitude) / 5.0 + n_gwas / 20.0 + min((n_pathologies - n_variant_pathologies), 30) / 60.0 + \
        n_pathogenic / 5.0 + min((n_variant_pathologies - n_pathogenic), 10) / 20.0


def main(input_filename, output_filename):
    dna = auto_load_dna(input_filename)
    genotypes = load_genotypes()
    gwas_traits = get_gwas_traits()
    pathology_mapping = get_clinvar_variant_pathologies()
    rs_pathologies = get_clinvar_rs_pathologies(pathology_mapping)
    variants_mapping = get_clinvar_variants(pathology_mapping)
    all_rss = get_summaries_dict(dna, genotypes, rs_pathologies, variants_mapping, gwas_traits)
    html = get_html_page(all_rss)
    with open(output_filename, "w") as f:
        f.write(html)
    print("Report written to " + output_filename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog="Personal DNA Analyzer",
        description="This program analysis your DNA to find interesting insights. Do not use for medical advice and"
                    "always consult your doctor."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input file. Accepted formats: MyHeritage, 23andMe, AncestryDNA.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output html file.")
    args = parser.parse_args()
    main(args.input, args.output)
