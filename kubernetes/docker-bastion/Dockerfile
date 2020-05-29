FROM ubuntu:18.04

RUN apt-get update -y && apt-get install -y \
    bzip2 \
    gcc \
    git \
    less \
    libncurses-dev \
    make \
    time \
    unzip \
    vim \
    wget \
    zlib1g-dev \
    liblz4-tool \
    zlib1g \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev

# installing python3.6 & pip
RUN apt-get update -y && apt-get install -y \
    python3.6-dev \
    python3-setuptools \
    && wget https://bootstrap.pypa.io/get-pip.py \
    && python3.6 get-pip.py

# supporting UTF-8
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# Install required python package
RUN pip install mysql-connector-python

## NFS setting
RUN apt-get -y update
RUN apt-get -y install nfs-common

RUN apt-get update && \
    apt-get install -y curl ca-certificates amqp-tools python \
       --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y gnupg2 gnupg
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && apt-get update -y && apt-get install google-cloud-sdk -y

WORKDIR /root

COPY sampleFiles/ sampleFiles/
COPY insertSamplesToRabbitMQ.sh .

CMD ["/bin/bash"]
