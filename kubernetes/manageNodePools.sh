#!/usr/bin/env bash

# Ref
# https://cloud.google.com/kubernetes-engine/docs/how-to/node-pools

#  Show nodes
gcloud container node-pools list --cluster ten-tcga --zone us-east1

# Create nodes
gcloud container node-pools create test-node --cluster ten-tcga --zone us-east1

# Resize nodes
gcloud container clusters resize cluster-name --node-pool pool-name \
        --num-nodes num-nodes

# Upgrade nodes
gcloud container clusters upgrade cluster-name --node-pool pool-name

# Delete node pools
gcloud container node-pools delete pool-name --cluster cluster-name

## Node Selector
# https://kubernetes.io/ko/docs/concepts/configuration/assign-pod-node/

# list nodes
kubectl get nodes

# Check executed node
kubectl get pods -o wide