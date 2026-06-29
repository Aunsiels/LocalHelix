import argparse
from collections import Counter

import pandas as pd
import numpy as np
from tqdm import tqdm
import json

from localhelix.dna_parsers import auto_load_dna
from localhelix.gwas import get_gwas_html, get_gwas_traits, initialize_gwas
from localhelix.snpedia import get_all_snpedia_match_genotypes, load_genotypes, initialize_snpedia, \
    get_snpedia_link, load_snps, set_snpedia_info, load_genosets
from localhelix.clinvar import get_clinvar_variant_from_rs, get_clinvar_variant_pathologies, \
    get_clinvar_rs_pathologies, get_clinvar_variants, print_pathologies_html, initialize_clinvar, get_clinvar_var_link, \
    get_haplotypes, find_all_haplotypes


def get_summaries_dict(dna, genotypes, snps, pathologies, variants_mapping, gwas_traits, haplotypes, var_to_pathology,
                       genosets):
    genotype_infos = get_all_snpedia_match_genotypes(dna, genotypes, snps, genosets)
    res = []
    all_clinvar_variants = []
    for key, value in tqdm(genotype_infos.items(),
                           total=len(genotype_infos), desc="Generating report"):
        variant_pathologies_list = []
        variant_pathology_names = set()
        all_pathologies = []
        variant_list = []
        if key.lower().startswith("rs"):
            base_rs = key.split("(")[0].lower()
            variant_list, pathos_list = get_clinvar_variant_from_rs(key, variants_mapping)
            for variant in variant_list:
                all_clinvar_variants.append(variant)
            all_pathologies = list(sorted(pathologies.get(base_rs, []), key=lambda x: x.name.lower()))
            for pathos in pathos_list:
                variant_pathologies_list.append(list(sorted(pathos, key=lambda x: x.name.lower())))
                variant_pathology_names.update({x.name for x in variant_pathologies_list[-1]})
            all_pathologies = [x for x in all_pathologies if x.name not in variant_pathology_names]
        temp = dict()
        temp["rs"] = key
        set_snpedia_info(value, temp)
        temp["ClinVarAllPathologies"] = all_pathologies
        temp["ClinVarVariants"] = variant_list
        temp["ClinVarVariantPathologies"] = variant_pathologies_list
        temp["GWAS"] = dict()
        if key.lower().startswith("rs"):
            base_rs = key.split("(")[0].lower()
            al1, al2 = key.split("(")[1].replace(")", "").split(";")
            for al in [al1, al2]:
                gwas_data = gwas_traits.get((base_rs, al))
                if gwas_data:
                    if al not in temp["GWAS"]:
                        temp["GWAS"][al] = []
                    # Structure as a list of objects for valid JSON
                    for item in gwas_data:
                        trait_info, count, or_betas = item[:3], item[3], item[4]
                        temp["GWAS"][al].append({"trait": trait_info, "count": count, "or_betas": or_betas})

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
        if rs.get("was_on_snpedia"):
            snpedia_counter += 1
        for pathology_list in rs.get("ClinVarVariantPathologies", []):
            for pathology in pathology_list:
                _add_pathology_to_counters_from_json(pathologies_variants, pathologies_variants_pathogenic, pathology)
        for pathology in rs.get("ClinVarAllPathologies", []):
            _add_pathology_to_counters_from_json(pathologies_all, pathologies_all_pathogenic, pathology)
        gwas = rs.get("GWAS", {})
        for al, value in gwas.items():
            for gwas_item in value:
                mapped_trait, _, text = gwas_item["trait"]
                or_betas = gwas_item["or_betas"]
                if not text:
                    for or_beta in or_betas:
                        if or_beta:
                            text = "Odd ratio >1" if float(or_beta) > 1.0 else "Odd ratio <1"
                            break  # Only need to set it once
                gwas_traits.append(mapped_trait + " (" + text + ")")

    pathologies_all = Counter(pathologies_all)
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


def _add_pathology_to_counters_from_json(pathologies_variants, pathologies_variants_pathogenic, pathology_json):
    name = pathology_json[1]
    is_pathogenic = pathology_json[2]
    if name and name not in ["not provided", "not specified"]:
        pathologies_variants.append(name)
        if is_pathogenic and "pathogenic" in is_pathogenic.lower():
            pathologies_variants_pathogenic.append(name)

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
    if rs.get("Repute", "") == "Bad":
        content = '<div class="card text-white bg-danger">'
    elif rs.get("Repute", "") == "Good":
        content = '<div class="card text-white bg-success">'
    else:
        content = '<div class="card bg-light">'
    return content


def _clean_for_json(obj):
    """Recursively clean an object for JSON serialization."""
    if isinstance(obj, dict):
        return {str(k): _clean_for_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clean_for_json(elem) for elem in obj]
    if pd.isna(obj):
        return None  # Handles np.nan, pd.NA, None, etc.
    if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
        return None
    if isinstance(obj, pd.Series):
        return _clean_for_json(obj.to_list())
    if isinstance(obj, pd.DataFrame):
        return _clean_for_json(obj.to_dict(orient='records'))
    return obj


def get_html_page(all_rss=None, json_data_url=None):
    summary_html = ""
    if all_rss:
        summary_html = summaries_results_in_html(all_rss)
    content = summary_html + r"""
        <h2>All variants</h2>
        <div class="row mb-3">
            <div class="col-md-8">
                <input type="text" id="search-bar" placeholder="Search variants by name, summary, or pathology..." class="form-control">
            </div>
            <div class="col-md-4">
                <select id="sort-options" class="form-select">
                    <option value="default" selected>Sort by: Relevance</option>
                    <option value="magnitude">Sort by: Magnitude (High to Low)</option>
                    <option value="repute">Sort by: Repute (Bad first)</option>
                </select>
            </div>
        </div>
        <div id="variants-container"></div>
        <nav aria-label="Page navigation">
          <ul class="pagination" id="pagination">
            <!-- Pagination links will be inserted here by JavaScript -->
          </ul>
        </nav>

        <script>
            const variantsContainer = document.getElementById('variants-container');
            const paginationContainer = document.getElementById('pagination');
            const searchBar = document.getElementById('search-bar');
            const sortOptions = document.getElementById('sort-options');
            const itemsPerPage = 50;
            let currentPage = 1;
            let all_variants_data = [];
            let filtered_variants_data = [];

            async function load_data() {
                const json_data_url = '""" + (json_data_url or '') + r"""';
                if (json_data_url) {
                    const response = await fetch(json_data_url);
                    all_variants_data = await response.json();
                } else {
                    all_variants_data = """ + (json.dumps(_clean_for_json(all_rss)) if all_rss else '[]') + r""";
                }
                filtered_variants_data = all_variants_data;
                // Initial render
                renderItems(1);
                setupPagination();
                searchBar.addEventListener('keyup', search);
                sortOptions.addEventListener('change', sortData);
            }

            // This is a trick: we create a JS function from an HTML template string.
            // All HTML generation is now done in JavaScript for clarity and to avoid escaping issues.
            function get_card_header(rs) {
                if (rs.Repute === "Bad") return '<div class="card text-white bg-danger">';
                if (rs.Repute === "Good") return '<div class="card text-white bg-success">';
                return '<div class="card bg-light">';
            }

            function get_snpedia_link(text) {
                if (text && (text.toLowerCase().startsWith("rs") || text.toLowerCase().startsWith("gs"))) {
                    return `<a class="link-dark" href="https://www.snpedia.com/index.php/${text}">${text}</a>`;
                }
                return "None";
            }

            function get_clinvar_var_link(text) {
                return !text ? "Unknown" : `<a class="link-dark" href="https://www.ncbi.nlm.nih.gov/clinvar/variation/${text}">${text}</a>`;
            }

            function get_clinvar_pathology_link(pathology) {
                // Defensive check to prevent crashes on bad data
                if (!pathology) {
                    return "Invalid Pathology Data";
                }

                let link;
                const split = pathology[0] ? pathology[0].split(":") : [];
                if (split.length === 2) {
                    link = `<a class="link-dark" href="https://www.ncbi.nlm.nih.gov/medgen/${split[1]}">${pathology[1]}</a>`;
                } else {
                    link = `<a class="link-dark" href="https://www.ncbi.nlm.nih.gov/medgen/?term=${pathology[1].replace(/ /g, "+")}">${pathology[1]}</a>`;
                }
                if (pathology[2]) {
                    if (pathology[4]) {
                        link += ` [<span title="${pathology[4]}">${pathology[3]}</span>]`;
                    } else {
                        link += ` [<span>${pathology[3]}</span>]`;
                    }
                }
                return link;
            }

            function print_pathologies_html(pathos) {
                if (!pathos || pathos.length === 0) return "None";
                const groups = {};
                for (const pathology of pathos) {
                    const is_pathogenic = pathology[2] || "Unknown";
                    if (!groups[is_pathogenic]) {
                        groups[is_pathogenic] = [];
                    }
                    groups[is_pathogenic].push(pathology);
                }
                if (Object.keys(groups).length === 0) return "None";
                let res = "<ul>";
                for (const group in groups) {
                    const pathologies_links = groups[group].map(p => get_clinvar_pathology_link(p)).join(", ");
                    res += `<li>${group}: ${pathologies_links}</li>`;
                }
                res += "</ul>";
                return res;
            }

            function get_gwas_html(rs) {
                if (!rs.GWAS || Object.keys(rs.GWAS).length === 0) return "<b>GWAS Traits</b>: None";
                let res = "<b>GWAS Traits</b>: <br><ul>";
                for (const al in rs.GWAS) {
                    res += `<li>${al}: `;
                    const temp = [];
                    (rs.GWAS[al] || []).forEach(gwas_item => {
                        const [mapped_trait, trait_uri, text_val] = gwas_item.trait;
                        const count = gwas_item.count;
                        const or_beta_list = gwas_item.or_betas;
                        let text = text_val ? `${text_val}, ` : text_val;
                        let or_beta_text = text ? "Reported Betas:" : "Reported <span title=\">1 accentuate the trait, <1 reduces the trait\">Odd Ratio</span>:";
                        or_beta_text += ` [${or_beta_list.join(", ")}]`;
                        temp.push(`<a class="link-dark" href="${trait_uri}">${mapped_trait}</a> (${text}${count} study(s), ${or_beta_text})`);
                    });
                    res += temp.join(", ");
                    res += "<br>";
                }
                res += `</ul><br><b>GWAS page</b>: <a class="link-dark" href="https://www.ebi.ac.uk/gwas/variants/${rs.rs.split("(")[0]}">${rs.rs.split("(")[0]}</a>`;
                return res;
            }

            function rs_to_html_js(rs) {
                if (!rs) return "";

                const clinvar_all_pathologies = print_pathologies_html(rs.ClinVarAllPathologies || []);
                const clinvar_variant_pathologies = (rs.ClinVarVariantPathologies || []).map(p => print_pathologies_html(p));
                const var_links = (rs.ClinVarVariants || []).map(v => get_clinvar_var_link(v));

                let rs_text = Array.isArray(rs.text) ? rs.text : [rs.text || ""];
                const snpedia_text = rs_text.map(x => x.replace(/href="\/index.php/g, 'href="https://www.snpedia.com/index.php')
                                                       .replace(/href="\/\/www.ncbi.nlm.nih.gov\/pubmed/g, 'href="https://www.ncbi.nlm.nih.gov/pubmed')
                                                       .replace(/href=/g, 'class="link-dark" href=')).join(" ");

                let content = get_card_header(rs);
                content += `
                  <h5 class="card-header">${rs.rs || "N/A"}</h5>
                  <div class="card-body">
                    <h5 class="card-title">${rs.summary}</h5>
                    <p class="card-text">${snpedia_text}
                    <br> <span title="Importance annotated by SNPedia community"><b>Magnitude</b></span>: ${rs.Magnitude}, <span title="Good or bad SNP, annotaed by SNPedia community"><b>Repute</b></span>: ${rs.Repute}<br>`;

                if (get_snpedia_link(rs.rs) !== "None") content += `<b>SNPedia Variant</b>: ${get_snpedia_link(rs.rs)}<br>`;
                if (rs.rs.includes("(") && get_snpedia_link(rs.rs.split("(")[0]) !== "None") content += `<b>SNPedia Base SNP</b>: ${get_snpedia_link(rs.rs.split("(")[0])}<br>`;
                if (rs.rs.toLowerCase().startsWith("r")) content += `<b>OpenSNP Link</b>: <a href="https://opensnp.org/snps/${rs.rs.split("(")[0]}">${rs.rs.split("(")[0]}</a><br>`;

                var_links.forEach((var_link, i) => {
                    content += `<span title="ClinVar page of your variant"><b>ClinVar Variant</b></span>: ${var_link}<br>`;
                    content += `<span title="Associated traits/pathologies with your variants on ClinVar"><b>Variants pathologies for ${var_link} (ClinVar)</b></span>: ${clinvar_variant_pathologies[i]}<br>`;
                });

                if (clinvar_all_pathologies !== "None") {
                    content += `<span title="Other pathologies associated with this SNP. Interesting to know if you have a good variant"><b>Pathologies/traits you avoided (ClinVar)</b></span>: ${clinvar_all_pathologies}<br>`;
                }
                if (rs.GWAS && Object.keys(rs.GWAS).length > 0) {
                    content += get_gwas_html(rs) + "<br>";
                }
                content += `</p></div></div>`;
                return content;
            }

            function sortData() {
                const sortBy = sortOptions.value;
                if (sortBy === 'default') {
                    // The default sort is the original order from the server (by relevance score)
                    // We just need to re-filter based on the current search term.
                    const searchTerm = searchBar.value.toLowerCase();
                    if (!searchTerm) {
                        filtered_variants_data = [...all_variants_data];
                    } else {
                        // Re-run search to get correctly ordered filtered data
                        search({target: {value: searchTerm}});
                    }
                } else if (sortBy === 'magnitude') {
                    filtered_variants_data.sort((a, b) => {
                        const magA = parseFloat(a.Magnitude) || 0;
                        const magB = parseFloat(b.Magnitude) || 0;
                        return magB - magA;
                    });
                } else if (sortBy === 'repute') {
                    const reputeOrder = { 'Bad': 0, 'Good': 1, 'Unknown': 2 };
                    filtered_variants_data.sort((a, b) => {
                        const reputeA = reputeOrder[a.Repute] ?? 2;
                        const reputeB = reputeOrder[b.Repute] ?? 2;
                        return reputeA - reputeB;
                    });
                }
                renderItems(1);
                setupPagination();
            }

            function search(event) {
                const searchTerm = event.target.value.toLowerCase();
                if (!searchTerm) {
                    filtered_variants_data = all_variants_data;
                } else {
                    filtered_variants_data = all_variants_data.filter(rs => {
                        // Search in variant name (rs), summary, and text
                        if (rs.rs && rs.rs.toLowerCase().includes(searchTerm)) return true;
                        if (rs.summary && rs.summary.toLowerCase().includes(searchTerm)) return true;
                        if (rs.text && rs.text.toLowerCase().includes(searchTerm)) return true;

                        // Search in ClinVar pathologies
                        if (rs.ClinVarAllPathologies) {
                            for (const p of rs.ClinVarAllPathologies) {
                                if (p[1] && p[1].toLowerCase().includes(searchTerm)) return true;
                            }
                        }
                        if (rs.ClinVarVariantPathologies) {
                            for (const p_list of rs.ClinVarVariantPathologies) {
                                for (const p of p_list) {
                                    if (p[1] && p[1].toLowerCase().includes(searchTerm)) return true;
                                }
                            }
                        }
                        return false;
                    });
                }
                sortData(); // Apply current sort to new search results
            }

            function renderItems(page) {
                currentPage = page;
                variantsContainer.innerHTML = '';
                const start = (page - 1) * itemsPerPage;
                const end = start + itemsPerPage;
                const paginatedItems = filtered_variants_data.slice(start, end);

                for (const item of paginatedItems) {
                    const card = document.createElement('div');
                    card.innerHTML = rs_to_html_js(item) + '<br>';
                    variantsContainer.appendChild(card);
                }
            }

            function setupPagination() {
                paginationContainer.innerHTML = '';
                const pageCount = Math.ceil(filtered_variants_data.length / itemsPerPage);
                for (let i = 1; i <= pageCount; i++) {
                    const li = document.createElement('li');
                    li.className = 'page-item' + (i === currentPage ? ' active' : '');
                    const a = document.createElement('a');
                    a.className = 'page-link';
                    a.href = '#';
                    a.innerText = i;
                    a.onclick = (e) => { e.preventDefault(); renderItems(i); setupPagination(); };
                    li.appendChild(a);
                    paginationContainer.appendChild(li);
                }
            }

            load_data();
        </script>
    """
    body = """<div class="container-fluid p-5 bg-primary text-white text-center"> <h1>LocalHelix - Your Personal DNA 
        Analyzer</h1> <p>Your results are ready!</p>
        <p><a href="/" class="btn btn-light btn-sm">Back to Upload Page</a></p>
        <form action="https://www.paypal.com/donate" method="post" 
        target="_top"> <input type="hidden" name="hosted_button_id" value="GR3D64Y7S7TU2" /> <input type="image" 
        src="https://www.paypalobjects.com/en_US/i/btn/btn_donate_LG.gif" border="0" name="submit" title="PayPal - 
        The safer, easier way to pay online!" alt="Donate with PayPal button" /> <img alt="" border="0" 
        src="https://www.paypal.com/en_FR/i/scr/pixel.gif" width="1" height="1" /> </form>

        </div>
          
        <div class="container mt-5">
        """ + content + """  
        </div>"""
    html = create_page_from_body(body)
    return html


def create_page_from_body(body):
    html = """
    <!DOCTYPE html>
        <html lang="en">
        <head>
          <title>LocalHelix - Your Personal DNA Analyzer</title>
          <meta charset="utf-8">
          <meta name="viewport" content="width=device-width, initial-scale=1">
          <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
          <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js"></script>
        </head>
        <body>
        """ + body + \
           """
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
        n_pathogenic / 5.0 + min((n_variant_pathologies - n_pathogenic), 10) / 20.0 or 0.0


def initialize_all(force=False, data_dir="data/"):
    initialize_gwas(force, data_dir)
    initialize_snpedia(force, data_dir)
    initialize_clinvar(force, data_dir)


def analyze_dna_to_json(input_filename, output_json_filename, force_reload=False, data_dir="data/", initialize=True):
    if initialize:
        initialize_all(force=force_reload, data_dir=data_dir)
    dna = auto_load_dna(input_filename)
    genotypes = load_genotypes(data_dir)
    snps = load_snps(data_dir)
    genosets = load_genosets(data_dir)
    gwas_traits = get_gwas_traits(data_dir)
    pathology_mapping = get_clinvar_variant_pathologies(data_dir)
    rs_pathologies = get_clinvar_rs_pathologies(pathology_mapping, data_dir)
    variants_mapping = get_clinvar_variants(pathology_mapping, data_dir)
    haplotypes = get_haplotypes(data_dir)
    all_rss = get_summaries_dict(dna, genotypes, snps, rs_pathologies, variants_mapping, gwas_traits, haplotypes,
                                 pathology_mapping, genosets)
    with open(output_json_filename, "w") as f:
        json.dump(_clean_for_json(all_rss), f)
    print("Report data written to " + output_json_filename)

def main(input_filename, output_filename, force_reload=False, data_dir="data/", initialize=True):
    if initialize:
        initialize_all(force=force_reload, data_dir=data_dir)
    dna = auto_load_dna(input_filename)
    genotypes = load_genotypes(data_dir)
    snps = load_snps(data_dir)
    genosets = load_genosets(data_dir)
    gwas_traits = get_gwas_traits(data_dir)
    pathology_mapping = get_clinvar_variant_pathologies(data_dir)
    rs_pathologies = get_clinvar_rs_pathologies(pathology_mapping, data_dir)
    variants_mapping = get_clinvar_variants(pathology_mapping, data_dir)
    haplotypes = get_haplotypes(data_dir)
    all_rss = get_summaries_dict(dna, genotypes, snps, rs_pathologies, variants_mapping, gwas_traits, haplotypes,
                                 pathology_mapping, genosets) # This returns a list of dicts
    html = get_html_page(all_rss)
    with open(output_filename, "w") as f:
        f.write(html)
    print("Report written to " + output_filename)


def get_arguments():
    parser = argparse.ArgumentParser(
        prog="LocalHelix - Your Personal DNA Analyzer",
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
