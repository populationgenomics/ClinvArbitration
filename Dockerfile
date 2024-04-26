FROM hailgenetics/hail:0.2.127-py3.11

COPY clinvarbitration /clinvarbitration
COPY requirements.txt /clinvarbitration/

RUN pip install --no-cache-dir -r /clinvarbitration/requirements.txt

WORKDIR /clinvarbitration
