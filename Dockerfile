FROM amancevice/pandas

RUN apt-get update \
    && apt-get install --no-install-recommends -y

RUN pip3 install --upgrade pip

RUN pip3 install numpy \
    matplotlib \
    sklearn \
    umap-learn \
    pydpc \
    scipy \
    scanpy

COPY *.py /usr/local/bin/
RUN chmod 755 /usr/local/bin/run-qc.py
