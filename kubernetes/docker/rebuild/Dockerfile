FROM rtea-tcga:1.8.6
MAINTAINER Boram Lee

WORKDIR /app

#COPY rnatea_pipeline_noscallop_from_bam rtea/

#RUN chmod +x /app/rtea/checkDiskSpace.sh

# Install google stroage driver
#RUN pip install google-cloud-storage

# Updated Boram's Code
#RUN R -e 'BiocManager::uninstall("GenomicFiles")'
#RUN R -e 'remove.packages("GenomicFiles")'
#COPY ref rtea/ref
#COPY run_ctea.R rtea/
#COPY rtea_functions.R rtea/

# COPY System Logger Files
#COPY kubernetes/persistentVolumes/checkHardwareUsage.sh rtea/
#RUN chmod +x /app/rtea/checkHardwareUsage.sh
#
#COPY kubernetes/persistentVolumes/checkMemoryUsage.sh rtea/
#RUN chmod +x /app/rtea/checkMemoryUsage.sh
#
#COPY kubernetes/persistentVolumes/checkDiskSpace.sh rtea/
#RUN chmod +x /app/rtea/checkDiskSpace.sh



# Install pika for RabbitMQ
#RUN pip install pika

#COPY kubernetes/worker.py rtea/
#COPY kubernetes/gdc-user-token-js.txt rtea/
RUN chmod 600 /app/rtea/gdc-user-token-js.txt
#COPY kubernetes/DBLogger.py rtea/
#COPY kubernetes/Config.py rtea/
#COPY kubernetes/aleelab-ten-key.json rtea/

#RUN ls -la rtea/*

# AWS SDK install
#RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
#RUN unzip awscliv2.zip
#RUN sudo ./aws/install
# aws s3 ls s3://1000genomes/1000G_2504_high_coverage/data/ERR3243161/NA18874.final.cram

#CMD /usr/bin/amqp-consume --url=$BROKER_URL -q $QUEUE -c 1 /app/rtea/worker.py
#CMD python /app/rtea/worker.py