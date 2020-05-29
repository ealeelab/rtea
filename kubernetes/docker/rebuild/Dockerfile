FROM rtea-tcga:1.2
MAINTAINER Boram Lee

WORKDIR /app

COPY ctea2rtea.R rtea/
COPY rtea_functions.R rtea/

CMD /usr/bin/amqp-consume --url=$BROKER_URL -q $QUEUE -c 1 /app/rtea/worker.py