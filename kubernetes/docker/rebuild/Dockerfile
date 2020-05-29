FROM rtea-tcga:1.2
MAINTAINER Boram Lee

WORKDIR /app

COPY boramFiles/ rtea/
COPY worker.py rtea/
COPY DBLogger.py rtea/

CMD /usr/bin/amqp-consume --url=$BROKER_URL -q $QUEUE -c 1 /worker/worker.py