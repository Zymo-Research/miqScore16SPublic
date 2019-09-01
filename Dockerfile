FROM bioconductor/release_core2:R3.5.2_Bioc3.8

#full install of python 3.6.7 on our bioconductor docker
RUN apt-get update &&\
    apt-get upgrade -y &&\
    apt-get install -y software-properties-common &&\
    mkdir python3.6 &&\
    cd python3.6 &&\
    wget https://www.python.org/ftp/python/3.6.7/Python-3.6.7.tgz &&\
    tar xvf Python-3.6.7.tgz &&\
    cd Python-3.6.7 &&\
    /python3.6/Python-3.6.7/configure --enable-optimizations --with-ensurepip=install &&\
    make -j8 &&\
    make altinstall &&\
    rm /usr/bin/lsb_release &&\
    ln -sf /usr/local/bin/python3.6 /usr/bin/python3 &&\
    ln -sf /usr/local/bin/pip3.6 /usr/bin/pip3 &&\
    pip3 install --upgrade pip

WORKDIR /
#move lsb_release and two following lines up to top unit during next major build
RUN pip3 install matplotlib &&\
    R -e 'install("getopt")' &&\
    R -e 'source("https://bioconductor.org/biocLite.R"); biocLite("dada2")' &&\
    cd /opt &&\
    mkdir miqscore16s

COPY . /opt/miqscore16s
RUN cd /opt/miqscore16s &&\
    pip3 install -r requirements.txt


CMD ["python3", "/opt/miqscore16s/analyzeStandardReads.py"]