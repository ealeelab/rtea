FROM boramlee/rtea
MAINTAINER Boram Lee


RUN apt-get update && apt-get install -y openjdk-8-jdk pigz
RUN wget https://github.com/broadinstitute/picard/releases/download/2.21.7/picard.jar
COPY rnatea_pipeline_noscallop_from_bam rtea/
