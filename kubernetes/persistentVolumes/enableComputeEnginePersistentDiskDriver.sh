#!/usr/bin/env bash

gcloud beta container clusters update ten-tcga \
      --update-addons=GcePersistentDiskCsiDriver=ENABLED