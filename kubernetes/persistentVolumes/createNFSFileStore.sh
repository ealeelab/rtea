#!/usr/bin/env bash


gcloud filestore instances create nfs-server \
    --project=aleelab-ten \
    --zone=us-west1-b \
    --tier=PREMIUM \
    --file-share=name="pv",capacity=8TB \
    --network=name="ten-vpc"