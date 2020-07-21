#!/usr/bin/env bash

export TARGET_POD_GROUP="rtea-tcga-test-deployment"
CONTAINER_NAME="c-rtea-tcga-test"

PodName=`kubectl get pod | grep $TARGET_POD_GROUP | awk '{print $1}'`

kubectl exec -it $PodName --container $CONTAINER_NAME -- /bin/bash

