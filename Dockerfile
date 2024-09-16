FROM python:3.10-slim-buster

WORKDIR /analyzer

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY personal_dna_analyzer/*py personal_dna_analyzer/
COPY personal_dna_analyzer/data/*tsv personal_dna_analyzer/data/
COPY personal_dna_analyzer/data/*db personal_dna_analyzer/data/
COPY personal_dna_analyzer/data/*json personal_dna_analyzer/data/

CMD ["flask", "--app", "personal_dna_analyzer/web_server.py", "run", "--host=0.0.0.0"]
