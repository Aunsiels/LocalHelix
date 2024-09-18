FROM python:3.10-slim-buster

WORKDIR /analyzer

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY localhelix/*py localhelix/
COPY localhelix/data/*tsv localhelix/data/
COPY localhelix/data/*db localhelix/data/
COPY localhelix/data/*json localhelix/data/

CMD ["flask", "--app", "localhelix/web_server.py", "run", "--host=0.0.0.0"]
