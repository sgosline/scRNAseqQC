FROM ubuntu:18.04

RUN apt-get update \
    && apt-get install --no-install-recommends -y \
       python3 \
       python3-setuptools \
       python3-pip \
       python3-pandas 


RUN pip3 install --upgrade pip

RUN pip3 install numpy \
    matplotlib \
    sklearn \
    umap \
    pydpc \
    scipy

COPY *.py /usr/local/bin/
RUN chmod 755 /usr/local/bin/run-qc.py
