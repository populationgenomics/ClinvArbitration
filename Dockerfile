FROM hailgenetics/hail:0.2.127-py3.11

COPY scripts /scripts
COPY requirements.txt /scripts/

RUN pip install --no-cache-dir -r /scripts/requirements.txt

WORKDIR /scripts
