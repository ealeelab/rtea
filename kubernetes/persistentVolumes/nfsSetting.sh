#!/usr/bin/env bash

## NFS setting
apt-get -y update
apt-get -y install nfs-common
mkdir /nfs-data
mount 10.21.191.186:/pv /nfs-data
chmod go+rw /nfs-data