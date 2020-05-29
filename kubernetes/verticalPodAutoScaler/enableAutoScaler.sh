#!/usr/bin/env bash

# Create
#gcloud container clusters create cluster-name \
#        --enable-vertical-pod-autoscaling --cluster-version=1.14.7

# Enable
gcloud container clusters update ten-tcga --enable-vertical-pod-autoscaling --region us-east1

# Create
kubectl create -f my-opt-vpa.yaml

# Check
kubectl get pod pod-name --output yaml
kubectl get pods
kubectl get vpa my-opt-vpa --output yaml

# Delete
kubectl delete vpa rtea-tcga-test-vpa2