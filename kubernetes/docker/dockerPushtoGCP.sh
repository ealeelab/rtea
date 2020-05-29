#!/usr/bin/env bash

VERSION="1.2"
IMAGES_NAME="rtea-tcga"

PROJECT_NAME="aleelab-ten"

docker tag $IMAGES_NAME:$VERSION gcr.io/$PROJECT_NAME/$IMAGES_NAME:$VERSION
gcloud docker -- push gcr.io/$PROJECT_NAME/$IMAGES_NAME:$VERSION
