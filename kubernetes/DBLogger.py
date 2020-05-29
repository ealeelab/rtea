import mysql.connector

def openLog(db,cursor,sampleName,fileName,queue):
    sno = 0
    try:
        # Determine open mode : Create or Update
        cursor.execute(
            "SELECT sno, file_name, file_status,file_access_count FROM logger WHERE file_name = %s",
            (fileName,)
        )
        row = cursor.fetchone()
        row_count = cursor.rowcount

        ## Create Mode
        if row_count == 0:
            sql = "INSERT INTO logger " \
                  " (sample_name,file_name,file_status,file_access_count,queue,processed_datetime)" \
                  " VALUES (%s,%s,%s,%s,%s,NOW())"
            values = (sampleName,fileName,'P',1,queue)
            cursor.execute(sql,values)
            db.commit()
            sno = cursor.lastrowid

        ## Update Mode (FileName is unique)
        elif row_count == 1:
            sno = row[0]
            file_access_count = row[3] + 1
            sql = "UPDATE logger SET file_access_count = %s, processed_datetime = NOW() WHERE sno = %s"
            values = (file_access_count,sno)
            cursor.execute(sql,values)
            db.commit()

    except mysql.connector.Error as err:
        print("Failed open log:{}".format(err))
        exit(1)

    return sno

def closeLog(db, cursor, sno):
    try:
        sql = "UPDATE logger SET file_status = 'C', completed_datetime = now() WHERE sno = %s "
        cursor.execute(sql,(sno,))
        db.commit()

    except mysql.connector.Error as err:
        print("Failed open log:{}".format(err))
        exit(1)


# dbConfig = {
#     'user' : 'aleelab',
#     'password': 'aleelabAprom3!',
#     'host': '34.73.61.94',
#     'database': 'logger_db',
#     'raise_on_warnings': True
# }

# db = mysql.connector.connect(**dbConfig)
# cursor = db.cursor()

# Test Codes
#sno = openLog(db,cursor,"TestSample","TestFile")
#print(sno)

#
# closeLog(db,cursor,1)
#
# db.close()





