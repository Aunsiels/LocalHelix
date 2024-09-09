# Personal DNA Analyzer

This program analyzes your DNA to find interesting insights.
**Do not use for medical advice and always consult your doctor**.

The main goal of this project is to get insights on our DNA without having to
send our data to unknown companies.

## How to run the program

To run the program, you need to have Python installed and a package installer like pip.

To install the dependencies, run the command (we recommend using a virtual environment):

```bash
pip install -r requirements.txt
```

Then, run the command:

```bash
python -m personal_dna_analyzer -i INPUT_FILE -o OUTPUT_FILE -d DATA_DIR
```

The input file is an SNP description file. The supported formats are 23andMe,
Ancestry.com, and MyHeritage. To include more format, please open an issue and
copy-paste the header of you file (until the names of the columns for a CSV/TSV).


The data directory will contain all the data required to generate your report.
If an empty directory is provided, the program will download and preprocess the
necessary documents. The process takes a lot of time. Otherwise, you can use an
existing data directory. We provide such a directory on our GitHub. You can
download it, unzip it, and link to it when you start your program.

The output file is a html file that you can open with your favorite web browser.

## Data Sources

Currently, we include the following data sources:
* SNPedia
* ClinVar
* GWAS

When possible, we link our results to the original source. We strongly recommend
to check it.

If you want us to include additional sources or information, please open an issue
on our GitHub.

## Support our work

[![](https://www.paypalobjects.com/en_US/i/btn/btn_donateCC_LG.gif)](https://www.paypal.com/donate/?hosted_button_id=GR3D64Y7S7TU2)