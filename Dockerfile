FROM python:3.6

RUN apt-get update && \ 
    apt-get upgrade -y && \
    apt-get install -y software-properties-common gnupg2 ca-certificates && \
    echo "deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/" >> /etc/apt/sources.list && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 'B8F25A8A73EACF41' && \
    apt-get update && \
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
