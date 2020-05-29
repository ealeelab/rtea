#!/usr/bin/env python
# Kubernetes Pod worker for TCGA data processing
# Author: Junseok Park (junseok.park@childrens.harvard.edu)

import sys
import time
import os.path
from os import path
import glob
import DBLogger as logger
import mysql.connector

# 0. Configure default system prameter
workPath = os.getenv('WORK_PATH')
outputPath = os.getenv('OUTPUT_PATH')
tempPath = os.getenv('TEMP_PATH')
samplePath = os.getenv('SAMPLE_PATH')
targetFileExt = os.getenv('TARGET_FILE_EXT')
refDir = os.getenv('REF_PATH')
queue = os.getenv("QUEUE")

# 0.0 Database Configuration
dbUsername = os.getenv("DB_USER")
dbPassword = os.getenv("DB_PASSWORD")
dbHost = os.getenv("DB_HOST")
dbName = os.getenv("DB_NAME")

os.chdir(workPath)

credentialFileName = "gdc-user-token.txt"
#refDir = "/worker/ref/hg38"

dbConfig = {
    'user' : dbUsername,
    'password': dbPassword,
    'host': dbHost,
    'database': dbName,
    'raise_on_warnings': True
}

db = mysql.connector.connect(**dbConfig)
cursor = db.cursor()

def mkPath(directory):
    if path.exists(directory):
        cmd = "rm -rf " + directory
        os.system(cmd)

    os.mkdir(directory)

def copyFile(src,dest,fileName):
    cmd = "cp " + src + "/" + fileName + " " + dest + "/"
    os.system(cmd)

def mkPathWithoutReplace(directory):
    if path.exists(directory):
        print("Path exists:"+directory)
    else:
        os.mkdir(directory)

try:
    # 1. Get sample Information
    sampleName = sys.stdin.readlines()[0]
    #sampleName = "88ee67cf-584b-44df-8204-a7e737f44e4e"
    mkPath(tempPath)
    mkPathWithoutReplace(outputPath)

    # 2. Download sample to the designated location of this container
    print("Downloading " + sampleName)

    # 2.1 File exist check
    fileDownload = True
    #if path.exists(os.path.join(samplePath,sampleName)):
    #    if glob.glob(os.path.join(samplePath,sampleName)+"/*."+targetFileExt):
    #        # We do not check the Hash of the file. (Should be updated)
    #        print("File already exists")
    #        fileDownload = False
    # Due to the performance result of FileStore(NFS), we will not use this function

    cmd = "rm -rf " + os.path.join(workPath,sampleName)
    os.system(cmd)

    if fileDownload:
        #cmd = "rm -rf " + os.path.join(samplePath,sampleName)
        #os.system(cmd)
        cmd = 'gdc-client download '+ sampleName + ' -t ' + credentialFileName
        os.system(cmd)
        #print("Copy downloaded file to NFS")
        #cmd = "cp -R " + os.path.join(workPath,sampleName) + " " + os.path.join(samplePath,sampleName)
    else:
        cmd = "cp -R " + os.path.join(samplePath,sampleName) + " " + os.path.join(workPath)
    #os.system(cmd)
    fileNames = glob.glob(os.path.join(workPath,sampleName)+"/*."+targetFileExt)

    if (len(fileNames) == 0):
        raise Exception("Sample download failed. Target file not exist.")
    fileName = fileNames[0]
    fileName = fileName.replace("./"+sampleName+"/","")

    # 3. Initiate data processing
    print("Processing " + fileName)

    # 3.1 File Process logging (Init)
    sno = logger.openLog(db,cursor,sampleName,fileName, queue)

    # 3.2 File process command
    currentPath=os.getcwd()
    targetFileWithPath = os.path.join(currentPath,sampleName,fileName)
    sampleOutputPath = os.path.join(tempPath,sampleName)
    mkPath(sampleOutputPath)

    print(targetFileWithPath)
    cmd = 'rnatea_pipeline_noscallop_from_bam ' \
          ' ' + targetFileWithPath + \
          ' ' + sampleName + \
          ' ' + refDir + \
          ' 4 ' \
          + sampleOutputPath +  \
          ' hg38 2>&1 | tee ' + sampleOutputPath + '/job.log'

    print("Following command will be executed:")
    print(cmd)
    os.system(cmd)

    print("Copy result file to NFS storage")
    storePath = os.path.join(outputPath,sampleName)
    mkPath(storePath)

    copyFile(sampleOutputPath,storePath,"bam/"+sampleName+".fastp.json")
    copyFile(sampleOutputPath,storePath,"bam/"+sampleName+".hisat2.bam")
    copyFile(sampleOutputPath,storePath,"bam/"+sampleName+".hisat2.bam.bai")
    copyFile(sampleOutputPath,storePath,"rtea/"+sampleName+".rtea.tx*")
    copyFile(sampleOutputPath,storePath,"job.log")

    # 3.3 File process logging (End)
    logger.closeLog(db,cursor,sno)

    # Remove Samplefile
    cmd = "rm -rf " + os.path.join(samplePath,sampleName)
    os.system(cmd)

    # Distribute the result to GCP Bucket
    # docker run -v ~/.config/gcloud:/root/.config/gcloud your_docker_image

except mysql.connector.Error as err:
    print("MySQL data connector error: {}".format(err))

except:
    print("Unexpected error:", sys.exc_info()[0])

finally:
    db.close()
    time.sleep(10)