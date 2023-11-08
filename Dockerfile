FROM python:3.6

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y software-properties-common && \
    apt-get install -y r-base && \
    ln -sf /usr/bin/Rscript /usr/local/bin/Rscript

WORKDIR /

RUN pip3 install matplotlib && \
    R -e 'install.packages("getopt")' && \
    R -e 'install.packages("BiocManager")' && \
    R -e 'BiocManager::install("dada2")' && \
    cd /opt && \
    mkdir miqscore16s

COPY . /opt/miqscore16s

RUN cd /opt/miqscore16s &&\
    pip3 install -r requirements.txt

CMD ["python3", "/opt/miqscore16s/analyzeStandardReads.py"]