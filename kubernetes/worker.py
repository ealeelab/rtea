#!/usr/bin/env python
# Kubernetes Pod worker for TCGA data processing
# Author: Junseok Park (junseok.park@childrens.harvard.edu)

import sys
import time
import os.path
from os import path
import glob
from DBLogger import Logger
import mysql.connector
from google.cloud import storage
from Config import MySQLConfig
import pika
import random

# 0. Configure default system prameters
# Job Information
projectName = os.getenv("PROJECT_NAME")
totalTaskCounts = os.getenv("TOTAL_TASK_COUNTS")

# RabbitMQ Parameters
queue = os.getenv("QUEUE")
brokerURL = os.getenv("BROKER_URL")

# Local Paths
workPath = os.getenv('WORK_PATH')
tempPath = os.getenv('TEMP_PATH')
targetFileExt = os.getenv('TARGET_FILE_EXT')

# System Parameters
cores = os.getenv("MULTI_THREADS")

# Reference GCP bucket setting
refDir = os.getenv('REF_PATH')
ref_bucket = os.getenv("REF_BUCKET")
ref_prefix = os.getenv("REF_PREFIX")

# Output GCP bucket setting
bucketName = os.getenv("OUTPUT_BUCKET")

# 0.0 Database Configuration
dbUsername = os.getenv("DB_USER")
dbPassword = os.getenv("DB_PASSWORD")
dbHost = os.getenv("DB_HOST")
dbName = os.getenv("DB_NAME")

os.chdir(workPath)

credentialFileName = "gdc-user-token-js.txt"
#refDir = "/worker/ref/hg38"

MySQLConfig.DATABASE_CONFIG = {
    'user' : dbUsername,
    'password': dbPassword,
    'host': dbHost,
    'database': dbName,
    'raise_on_warnings': True
}


def getOneTaskFromQueue(queue, channel, timeoutSec=5):

    body = None
    for method_frame, properties, body in channel.consume(queue,inactivity_timeout=timeoutSec):

        if body is None:
            print("There is no task on the queue")
            return body
        # Acknowledge the message
        channel.basic_ack(method_frame.delivery_tag)
        # Escape out of the loop after 1 messages
        if method_frame.delivery_tag == 1:
            break
    return body.decode()

def insertOneTaskToQueue(channel,queue, message):
    channel.queue_declare(queue=queue,durable=True)
    channel.basic_publish(
        exchange='',
        routing_key=queue,
        body=message,
        properties=pika.BasicProperties(
            delivery_mode=2,  # make message persistent
        ))
    print("{} is delivered to {}".format(message,queue))


def dumpQueueAndTasksFromFile(channel, queue, filePath):
    channel.queue_declare(queue=queue,durable=True)
    with open(filePath) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            print("Line {}: ".format(cnt))
            insertOneTaskToQueue(channel,queue,line.strip())
            line = fp.readline()
            cnt += 1

def uploadFileToBucket(bucketName,queue,sampleName,fileName,fileNameWithPath):
    client = storage.Client()
    bucket = client.get_bucket(bucketName)
    blob = bucket.blob(queue+"/"+sampleName+"/"+fileName)
    blob.upload_from_filename(fileNameWithPath)

    print(
        "Blob {} uploaded to {}. ".format(
            fileNameWithPath,fileName
        )
    )

def downloadFileFromBucket(bucketName,sourceBlobName,destinationFileName):

    client = storage.Client()
    bucket = client.get_bucket(bucketName)
    blob = bucket.blob(sourceBlobName)
    blob.download_to_filename(destinationFileName)

    print(
        "Blob {} downloaded to {}. ".format(
            sourceBlobName,destinationFileName
        )
    )

def downloadDirectoryFromBucket(bucketName,sourcePrefix,downloadPath):
    client = storage.Client()
    bucket = client.get_bucket(bucketName)
    blobs = bucket.list_blobs(prefix=sourcePrefix)
    if downloadPath[-1] != '/':
        downloadPath=downloadPath+'/'

    mkPath(downloadPath)

    for blob in blobs:
        filename = blob.name.split('/')[-1]
        blob.download_to_filename(downloadPath + filename)
        print(
            "Blob {} downloaded to {}. ".format(
                filename,downloadPath+filename
            )
        )

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

#def checkAndCopyResult(sampleOutputPath,storePath,fileName):
def checkAndCopyResult(sampleOutputPath,sampleName,fileName):
    if path.exists(os.path.join(sampleOutputPath,fileName)):
        uploadFileToBucket(bucketName, queue,sampleName,fileName,os.path.join(sampleOutputPath,fileName))
        #copyFile(sampleOutputPath,storePath,fileName)
        return 1
    else:
        return 0

sampleName = ""
sampleOutputPath = ""
resultCode = 'F'
jobNo = 0
taskNo = 0
finalize = True
try:
    # 1. Get sample Information

    # Random waiting
    time.sleep(random.randint(1,60))
    connection = pika.BlockingConnection(pika.URLParameters(brokerURL))
    channel = connection.channel()
    channel.cancel()
    sampleName = getOneTaskFromQueue(queue,channel)
    #sampleName = "88ee67cf-584b-44df-8204-a7e737f44e4e"
    #sampleName = "3fb7de4b-c3b0-49c9-a70d-acb1662f64b9"
    if sampleName is None or sampleName.strip() is "":
        raise SystemExit
    print("Consume from Queue: ", sampleName)
    channel.close()
    connection.close()

    # Check Job and Task status from Logger
    notCompletionStatus = False

    logger = Logger(MySQLConfig)
    (notCompletionStatus,jobNo,taskNo) = logger.openLog(projectName,queue,totalTaskCounts,sampleName,sampleName,queue)
    if notCompletionStatus is False:
        raise SystemExit

    mkPath(tempPath)

    ## Record Hardware usage for Container
    import subprocess
    print("execute hardware usage monitoring deamon")
    hardwareLogger=subprocess.Popen([workPath+"/checkHardwareUsage.sh"])
    print("execute memory usage monitoring deamon")
    memoryLogger=subprocess.Popen([workPath+"/checkMemoryUsage.sh"])
    print("execute disk space monitoring deamon")
    diskLogger=subprocess.Popen([workPath+"/checkDiskSpace.sh"])


    # 1.1. Download reference file from bucket
    print("Downloading reference files..")
    downloadDirectoryFromBucket(ref_bucket,ref_prefix,refDir)

    # 2. Download sample to the designated location of this container
    print("Downloading " + sampleName)

    # 2.1 File exist check
    # Remove redundant files
    cmd = "rm -rf " + os.path.join(workPath,sampleName)
    os.system(cmd)

    # User custom command for download
    cmd = 'gdc-client download '+ sampleName + ' -t ' + credentialFileName
    os.system(cmd)
    fileNames = glob.glob(os.path.join(workPath,sampleName)+"/*."+targetFileExt)

    if (len(fileNames) == 0):
        raise Exception("Sample download failed. Target file not exist.")
    fileName = fileNames[0]
    fileName = fileName.replace("./"+sampleName+"/","")

    # 3. Initiate data processing
    print("Processing " + fileName)

    # 3.1 File Process logging (Init)
    if sampleName is not fileName:
        print("Update file name in queue")
        logger.updateTaskFileName(taskNo,fileName)
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
          ' ' + cores + \
          ' ' + sampleOutputPath +  \
          ' hg38 2>&1 | tee ' + sampleOutputPath + '/job.log'

    print("Following user defined command will be executed:")
    print(cmd)
    os.system(cmd)
    hardwareLogger.kill()
    memoryLogger.kill()
    diskLogger.kill()


except mysql.connector.Error as err:
    print("MySQL data connector error: {}".format(err))

except SystemExit:
    print("System exit")
    finalize = False

except:
    print("Unexpected error:", sys.exc_info()[0])
    ## Insert Failed Message
    if sampleName.strip() is not "" and sampleName is not None:
        print("Re-insert failed message to the queue: Not doing this time")
        connection = pika.BlockingConnection(pika.URLParameters(brokerURL))
        channel = connection.channel()
        channel.cancel()
        insertOneTaskToQueue(channel,queue,sampleName)
        channel.close()
        connection.close()
finally:
    if finalize is True:
        print("Start to check result files")
        ## Check Result Status
        resultScore = 0
        print("ResultScore:",resultScore)

        ## Upload Results to GCP Bucket
        print("Copy files to GCP bucket")
        resultScore += checkAndCopyResult(sampleOutputPath+"/bam",sampleName,sampleName+".fastp.json")
        resultScore += checkAndCopyResult(sampleOutputPath+"/bam",sampleName,sampleName+".hisat2.bam")
        resultScore += checkAndCopyResult(sampleOutputPath+"/bam",sampleName,sampleName+".hisat2.bam.bai")
        resultScore += checkAndCopyResult(sampleOutputPath+"/rtea",sampleName,sampleName+".rtea.txt")
        print("ResultScore:",resultScore)

        if (len(glob.glob(os.path.join(sampleOutputPath,"rtea",sampleName+".rtea.txt.fail"))) > 0):
            resultScore = 0
            checkAndCopyResult(sampleOutputPath+"/rtea",sampleName,sampleName+".rtea.txt.fail")

        ## Upload Log Files to GCP Bucket
        if (len(glob.glob(os.path.join(sampleOutputPath,"job.log"))) > 0):
            print("Upload log files to the GCP bucket")
            uploadFileToBucket(bucketName,queue,sampleName,"job.log",os.path.join(sampleOutputPath,"job.log"))
            uploadFileToBucket(bucketName,queue,sampleName,"hardwareUsage.log",os.path.join(workPath,"hardwareUsage.log"))
            uploadFileToBucket(bucketName,queue,sampleName,"memoryUsage.log",os.path.join(workPath,"memoryUsage.log"))
            uploadFileToBucket(bucketName,queue,sampleName,"diskUsage.log",os.path.join(workPath,"diskUsage.log"))

            # 3.3 File process logging (End)
            errReason=None
            if resultScore >= 3:
                resultCode = 'C'
                print("Worker succeed the task")
            else:
                cmd='dmesg --level=err > systemError.log'
                os.system(cmd)
                uploadFileToBucket(bucketName,queue,sampleName,"systemError.log",os.path.join(workPath,"systemError.log"))
                print("Woker failed the task")
                errReason = "user pipeline error occurred"
                if os.path.getsize(os.path.join(workPath,"systemError.log")) > 0:
                    print("System error occurred")
                    errReason = "system error occurred. Please see the systemError.log"
            logger = Logger(MySQLConfig)
            logger.closeLog(jobNo,taskNo,resultCode,projectName,queue,"gs://"+bucketName+"/"+queue+"/"+sampleName,errReason)

        print("Worker finished")
    else:
        print("Nothing to do. Worker will be closed")
    time.sleep(10)