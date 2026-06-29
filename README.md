# LocalHelix - Your Personal DNA Analyzer


>In a small village not so far away, all the Gauls shared their data with Roman companies... All except one! 💪
Meet **LocalHelix**, the brave software powered by privacy potion, keeping your DNA safe from the prying eyes of
Rome’s data-mining legionaries. 🛡️ Unlike Cacofonix, it won’t share your secrets with anyone! 😉


This program analyzes your DNA to find interesting insights.
**Do not use for medical advice and always consult your doctor**.

The main goal of this project is to get insights on our DNA without having to
send our data to unknown companies.

## How to run the program

Before anything, you need to download this project locally, either with Git or
directly using the download button on the GitHub page.

### With Docker

You can run the project with Docker. First, install Docker following
the instructions on the [official website](https://docs.docker.com/engine/install/).
Then, you can optionally download the data folder provided on the release page of
our GitHub project. You must place this folder in the localhelix folder
(localhelix/data). Then run at the root of the project:

```bash
docker build --tag dna_analyzer .
```

and

```bash
docker run dna_analyzer
```

If everything went smoothly, you should see a URL appear in the log that looks like:

```
Running on http://172.17.0.3:5000
```

Open your browser and go to the mentioned URL (http://172.17.0.3:5000 here). On
this page, you have to upload your DNA results (this website is on your machine,
so your data is not sent to anybody). After clicking on the upload button, wait a
few minutes and your results will appear!

### With Python and Command Line

To run the program, you need to have Python installed and a package installer like pip.

To install the dependencies, run the command (we recommend using a virtual environment):

```bash
pip install -r requirements.txt
```

Then, run the command:

```bash
python -m localhelix -i INPUT_FILE -o OUTPUT_FILE -d DATA_DIR
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

### With Python and Flask

Install the dependencies with:

```bash
pip install -r requirements.txt
```

Then, run the commands

```bash
export FLASK_APP=localhelix/web_server.py
flask run
```

You can now open your browser, go to http://127.0.0.1:5000 and upload your DNA.

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