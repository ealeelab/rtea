FROM r-base:3.6.2
MAINTAINER Boram Lee

ENV LD_LIBRARY_PATH=/usr/local/lib
RUN apt-get update && apt-get install -y \
    cmake \
    libxml2-dev \
    libcurl4-openssl-dev \
    libboost-dev \
    gawk \
    libssl-dev \
    pigz \
    htop
WORKDIR /app

# fastp
RUN wget -P /usr/local/bin http://opengene.org/fastp/fastp \
    && chmod a+x /usr/local/bin/fastp

# hisat2
RUN wget -P /app --content-disposition https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download \
    && unzip /app/hisat2-2.1.0-Linux_x86_64.zip \
    && rm /app/hisat2-2.1.0-Linux_x86_64.zip
ENV PATH=/app/hisat2-2.1.0:$PATH

# samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar -xjf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && ./configure \
    && make \
    && make install \
    && cd /app \
    && rm -r samtools-1.9.tar.bz2 samtools-1.9

# htslib for scallop
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
    && tar -xjf htslib-1.9.tar.bz2 \
    && cd htslib-1.9 \
    && ./configure --disable-bz2 --disable-lzma --disable-gcs --disable-s3 --enable-libcurl=no \
    && make \
    && make install \
    && cd /app \
    && rm -r htslib-1.9.tar.bz2 htslib-1.9

# scallop
RUN wget https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4_linux_x86_64.tar.gz \
    && tar -xzf scallop-0.10.4_linux_x86_64.tar.gz \
    && mv scallop-0.10.4_linux_x86_64/scallop /usr/local/bin \
    && rm -r scallop-0.10.4_linux_x86_64.tar.gz scallop-0.10.4_linux_x86_64

# bamtools
RUN wget https://github.com/pezmaster31/bamtools/archive/v2.5.1.tar.gz \
    && tar -xzf /app/v2.5.1.tar.gz \
    && mkdir bamtools-2.5.1/build \
    && cd bamtools-2.5.1/build \
    && cmake -DBUILD_SHARED_LIBS=ON .. \
    && make \
    && make install \
    && cd /app \
    && rm -r v2.5.1.tar.gz bamtools-2.5.1
ENV PKG_CXXFLAGS="-I/usr/local/include/bamtools"
ENV PKG_LIBS="-L/usr/local/lib -lbamtools"

# bwa
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
    && tar -xjf bwa-0.7.17.tar.bz2 \
    && cd bwa-0.7.17 \
    && make \
    && cd /app \
    && rm bwa-0.7.17.tar.bz2
ENV PATH=/app/bwa-0.7.17:$PATH

# R packages
# R packages
RUN R -e "install.packages(c( \
            'magrittr', \
            'data.table', \
            'stringr', \
            'optparse', \
            'Rcpp', \
            'BiocManager' \
          ))"

RUN R -e "BiocManager::install(c( \
            'GenomicAlignments', \
            'BSgenome.Hsapiens.UCSC.hg19', \
            'BSgenome.Hsapiens.UCSC.hg38', \
            'EnsDb.Hsapiens.v75', \
            'EnsDb.Hsapiens.v86', \
          ))"

RUN apt-get update && apt-get install -y libssl-dev
#RUN R -e 'BiocManager::install("GenomicFiles")'

# rtea files from Boram
COPY ref rtea/
COPY ctea rtea/
COPY etc rtea/
COPY frombam rtea/
COPY ctea2rtea.R rtea/
COPY HISAT2_TE.sh rtea/
COPY isb_cgc_api_v3.py rtea/
COPY rnatea_pipeline_from_bam rtea/
COPY rnatea_pipeline_noscallop_from_bam rtea/
COPY rtea_functions.R rtea/
COPY rtea_pipeline rtea/
COPY rtea_pipeline_without_scallop rtea/
COPY run_ctea.R rtea/

ENV PATH=/app/rtea:$PATH

RUN apt-get update && apt-get install -y openjdk-8-jdk
RUN wget https://github.com/broadinstitute/picard/releases/download/2.21.7/picard.jar

# COPY System Logger Files
COPY kubernetes/persistentVolumes/checkHardwareUsage.sh rtea/
RUN chmod +x /app/rtea/checkHardwareUsage.sh

COPY kubernetes/persistentVolumes/checkMemoryUsage.sh rtea/
RUN chmod +x /app/rtea/checkMemoryUsage.sh

COPY kubernetes/persistentVolumes/checkDiskSpace.sh rtea/
RUN chmod +x /app/rtea/checkDiskSpace.sh

# Copy additional files for Kubernetes
COPY kubernetes/DBLogger.py rtea/
COPY kubernetes/Config.py rtea/
COPY kubernetes/aleelab-ten-key.json rtea/
COPY kubernetes/gdc-user-token.txt rtea/
COPY kubernetes/gdc-user-token-js.txt rtea/
RUN chmod 600 /app/rtea/gdc-user-token-js.txt

COPY kubernetes/gdc-client rtea/
COPY kubernetes/worker.py rtea/

# Set Environment variable for Worker

# Update Packages
RUN apt-get install -y software-properties-common
RUN apt-get update
RUN apt-get install -y emacs wget bzip2

# Anaconda installing
#RUN adduser --disabled-password -gecos '' debian
#RUN adduser debian sudo
#RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
#USER debian
#WORKDIR /home/debian/
#RUN chmod a+rwx /home/debian/

# Ref: http://www.science.smith.edu/dftwiki/index.php/Tutorial:_Docker_Anaconda_Python_--_4
RUN chmod a+rwx /opt
RUN wget https://repo.anaconda.com/archive/Anaconda3-5.3.1-Linux-x86_64.sh
RUN bash Anaconda3-5.3.1-Linux-x86_64.sh -b -p /opt/anaconda3
RUN rm Anaconda3-5.3.1-Linux-x86_64.sh

# Set path
ENV PATH=/opt/anaconda3/bin:$PATH

# Updating Anaconda Packages
RUN conda update conda
RUN conda update anaconda
RUN conda update --all

# Install required python package
RUN pip install mysql-connector-python
RUN pip install google-cloud-storage

## Install AMQP-Tools
RUN apt-get update && \
    apt-get install -y curl ca-certificates amqp-tools python \
       --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*

## NFS setting
RUN apt-get -y update
RUN apt-get -y install nfs-common
RUN pip install pika

# Packages for Google SDK
#RUN apt-get install gnupg gnupg2 gnupg1
# RUN apt-get install libssl1.0.0 libssl-dev

## Google Cloud SDK
#FROM alpine:3.6 RUN apk add --update \ python \ curl \ which \ bash RUN curl -sSL https://sdk.cloud.google.com | bash ENV PATH $PATH:/root/google-cloud-sdk/bin
# RUN curl -sSL https://sdk.cloud.google.com | bash
#RUN apt install sudo
#RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && apt-get update -y && apt-get install google-cloud-sdk -y

#CMD /usr/bin/amqp-consume --url=$BROKER_URL -q $QUEUE -c 1 /app/rtea/worker.py
CMD python /app/rtea/worker.py