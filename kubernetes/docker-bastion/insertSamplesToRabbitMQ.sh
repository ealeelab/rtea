#!/usr/bin/env bash
CANCER_TYPES=("COAD" "READ" "ESCA" "STAD" "LUSC" "LUAD" "HNSC" "UCEC" "OV" "CESC")

for i in "${CANCER_TYPES[@]}"
do
    export QUE_NAME=$i
    echo $QUE_NAME
    FILE_NAME="sampleFiles/${QUE_NAME}.txt"

    # Now create a queue:
    /usr/bin/amqp-declare-queue --url=$BROKER_URL -q $QUE_NAME -d
    while IFS= read -r line
    do
        # Publish one message to it:
        /usr/bin/amqp-publish --url=$BROKER_URL -r $QUE_NAME -p -b $line
        echo "$line"
    done < "$FILE_NAME"

done

# And get it back.
#/usr/bin/amqp-consume --url=$BROKER_URL -q foo -c 1 cat && echo

#kubectl exec shell-demo cat /proc/1/mounts