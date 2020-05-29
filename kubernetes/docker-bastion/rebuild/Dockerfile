FROM rtea-tcga-bastion:1.0

WORKDIR /root

COPY sampleFiles/ sampleFiles/
COPY insertSamplesToRabbitMQ.sh .
COPY consumeAllMessage.sh .

CMD ["/bin/bash"]
