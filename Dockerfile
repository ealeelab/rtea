FROM r-base:3.6.2
MAINTAINER junseok.park@childrens.harvard.edu

WORKDIR /app

RUN wget https://github.com/curl/curl/releases/download/curl-7_60_0/curl-7.60.0.tar.gz \
    && tar -xzf /app/curl-7.60.0.tar.gz \
    && cd curl-7.60.0 \
    && ./configure --with-openssl \
    && make \
    && make install \
    && cd /app \
    && rm -r curl-7.60.0.tar.gz curl-7.60.0

RUN wget --no-check-certificate https://www.gnupg.org/ftp/gcrypt/gnupg/gnupg-1.4.23.tar.gz \
    && tar -xzf /app/gnupg-1.4.23.tar.gz \
    && cd gnupg-1.4.23 \ 
    && ./configure \
    && make \
    && make install \
    && cd /app \
    && rm -r gnupg-1.4.23.tar.gz gnupg-1.4.23 

RUN apt-mark hold r-base

ENV LD_LIBRARY_PATH=/usr/local/lib

RUN gpg --keyserver keyserver.ubuntu.com --recv-keys 0E98404D386FA1D9 6ED0E7B82643E131
RUN gpg --export --armor 0E98404D386FA1D9 | apt-key add - \
    && gpg --export --armor 6ED0E7B82643E131 | apt-key add -

RUN apt-get update && apt-get install -y \
    cmake \
    libxml2-dev \
    libcurl4-openssl-dev \
    libboost-dev \
    gawk \
    libssl-dev \
    pigz \
    htop \
    iputils-ping


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
            'EnsDb.Hsapiens.v86' \
          ))"

# rtea
COPY . rtea/
ENV PATH=/app/rtea:$PATH
