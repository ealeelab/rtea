#!/usr/bin/env bash

# Ref: https://kubernetes.io/docs/tasks/debug-application-cluster/get-shell-running-container/

# Check Container
kubectl get pod {pod-name}
# If a Pod has more than one container
kubectl exec -it {pod-name} --container {container-name} -- /bin/bash

# Get Shell to the running Container
kubectl exec -it {pod-name} -- /bin/bash

# Running individual commands in a Container
kubectl exec {pod-name} {shell commands}

