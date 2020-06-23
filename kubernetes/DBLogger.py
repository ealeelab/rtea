import mysql.connector

class CompletedJobAndTaskError(Exception):
    def __init__(self, job, message="The job or Task is already Completed"):
        self.job = job
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.job} -> {self.message}'


class Logger:
    queryJobExistCheck = "SELECT sno FROM job WHERE project_name=%s AND job_name=%s"
    queryJobCompletionCheck = "SELECT sno, task_count_total, task_count_completed, task_count_failed, task_count_processing FROM job WHERE sno=%s"
    queryJobCreation = "INSERT INTO job (project_name,job_name,task_count_total) VALUES(%s,%s,%s)"
    queryJobProcessing = "UPDATE job SET task_count_processing = task_count_processing + 1 WHERE sno =%s"
    queryJobProcessingWithFailOver = "UPDATE job SET task_count_processing = task_count_processing + 1 " \
                                     " , task_count_failed = task_count_failed -1 WHERE sno =%s"

    queryJobCompleting = "UPDATE job SET task_count_processing = task_count_processing - 1 " \
                               ", task_count_completed = task_count_completed + 1 " \
                               ", last_processed_date = NOW() WHERE sno=%s"

    queryJobFailing = "UPDATE job SET task_count_processing = task_count_processing - 1 " \
                      ", task_count_failed = task_count_failed + 1 " \
                      ", last_processed_date = NOW() WHERE sno=%s"
    queryJobClosing = "UPDATE job SET ended_date = NOW() WHERE sno=%s"

    queryTaskExistCheck = "SELECT sno FROM logger WHERE sample_name=%s AND queue=%s"
    queryTaskCompletionCheck = "SELECT sno, file_status FROM logger WHERE sample_name = %s AND queue = %s"
    queryTaskStatusCheck = "SELECT file_status FROM logger WHERE sno=%s"
    queryTaskCreation = "INSERT INTO logger " \
                        " (sample_name,file_name,file_status,file_access_count,queue,processed_datetime, job_sno) " \
                        " VALUES (%s,%s,%s,%s,%s,NOW(),%s)"
    queryTaskUpdateFileName = "UPDATE logger SET file_name=%s WHERE sno = %s "
    queryTaskProcessing = "UPDATE logger SET file_access_count = file_access_count + 1, " \
                          "processed_datetime = NOW() WHERE sno = %s"
    queryTaskCompleting = "UPDATE logger SET file_status = %s, completed_datetime = now(), result_path = %s, error_reason = %s WHERE sno = %s "

    config = None

    def __init__(self,config):
        self.config = config
        self._conn = mysql.connector.connect(**config.DATABASE_CONFIG)
        self._conn.autocommit = True
        self._cursor=self._conn.cursor()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @property
    def connection(self):
        return self._conn

    @property
    def cursor(self):
        return self._cursor

    def commit(self):
        self.connection.commit()

    def close(self, commit=True):
        if commit:
            self.commit()
        self.connection.close()

    def execute(self, sql, params=None):
        try:
            self.cursor.execute(sql, params or ())
        except (AttributeError, mysql.connector.OperationalError):
            self._conn = mysql.connector.connet(**self.config.DATABASE_CONFIG)
            self._conn.autocommit = True
            self._cursor=self._conn.cursor()
            self.cursor.execute(sql, params or ())
        finally:
            return self.cursor

    def fetchall(self):
        return self.cursor.fetchall()

    def fetchone(self):
        return self.cursor.fetchone()

    def query(self, sql, params=None):
        self.cursor.execute(sql, params or ())
        return self.fetchall()

    def selectOne(self, sql, params=None):
        self.cursor.execute(sql + " limit 1", params or ())
        return self.fetchone()

    def transaction(self, sqlArray):
        try:
            self.connection.autocommit = False
            for sql in sqlArray:
                self.cursor.execute(sql)

        except mysql.connector.Error as error:
            print("Failed to update record to database rollback: {}".format(error))
            self.connection.rollback()

    def checkJobExist(self, projectName,jobName):
        self.cursor.execute(self.queryJobExistCheck,(projectName,jobName))
        row = self.fetchone()
        row_count = self.cursor.rowcount

        if row_count == 1:
            return row[0]  # Return SNO(Key) value
        else:
            return 0

    def checkJobCompletion(self,projectName,jobName):
        sno = self.checkJobExist(projectName,jobName)
        if sno == 0:
            return False
        else:
            row = self.selectOne(self.queryJobCompletionCheck,(sno,))
            if row[1] <= (row[2]+row[3]+row[4]): # 1:taskTotalCount, 2: taskCompletedCount, 3:taskFailedCount, 4: taskCountProcessing
                return True
            else:
                return False

    def createJob(self,projectName,jobName,taskTotal):
        self.execute(self.queryJobCreation,(projectName,jobName,taskTotal))
        return self.cursor.lastrowid

    def updateJobProcessing(self,sno):
        self.execute(self.queryJobProcessing,(sno,))

    def updateJobProcessingWithFailOver(self,sno):
        self.execute(self.queryJobProcessingWithFailOver,(sno,))

    def updateJobCompleting(self,sno,status):
        if status is "C":
            self.execute(self.queryJobCompleting,(sno,))
        else:
            self.execute(self.queryJobFailing,(sno,))

    def updateJobClosing(self,sno):
        self.execute(self.queryJobClosing,(sno,))


    def checkTaskExist(self,sampleName,queue):
        row = self.selectOne(self.queryTaskExistCheck,(sampleName,queue))
        row_count = self.cursor.rowcount

        if row_count == 1:
            return row[0] # Return SNO(Key) value
        else:
            return 0

    def checkTaskCompletion(self,sampleName,queue):
        sno = self.checkTaskExist(sampleName,queue)
        if sno == 0:
            return False
        else:
            row = self.selectOne(self.queryTaskCompletionCheck,(sampleName,queue))
            if row[1] is "C":
                return True
            else:
                return False

    def checkTaskStatus(self,sno):
        row = self.selectOne(self.queryTaskStatusCheck,(sno,))
        row_count = self.cursor.rowcount
        if row_count == 1:
            return row[0] # Return File Status (C,F)
        else:
            return 'N' # No Record, New

    def createTask(self,sampleName,fileName,queue,jobSno):
        self.execute(self.queryTaskCreation,(sampleName,fileName,'P',1,queue,jobSno))
        return self.cursor.lastrowid

    def updateTaskProcessing(self,sno):
        self.execute(self.queryTaskProcessing,(sno,))

    def updateTaskCompleting(self,sno,status,resultPath,errorReason):
        self.execute(self.queryTaskCompleting,(status,resultPath,errorReason,sno))

    def updateTaskFileName(self,sno,fileName):
        self.execute(self.queryTaskUpdateFileName,(fileName,sno))

    def openLog(self, projectName,jobName,taskTotal,
                sampleName,fileName,queue):
        jobNo = 0
        taskNo = 0
        try:

            taskNo = self.checkTaskExist(sampleName,queue)
            taskStatus = self.checkTaskStatus(taskNo)

            if self.checkTaskCompletion(sampleName,queue):
                raise CompletedJobAndTaskError(sampleName)

            if self.checkJobCompletion(projectName,jobName) and taskStatus is not 'P':
                raise CompletedJobAndTaskError(jobName)

            jobNo = self.checkJobExist(projectName,jobName)
            if jobNo == 0:
                jobNo = self.createJob(projectName,jobName,taskTotal)


            self.connection.autocommit = False
            ## Transaction Start (C:Completed is not working because is above line)
            if taskStatus is 'F': # Failed
                self.updateJobProcessingWithFailOver()
            elif taskStatus is 'N': # New
                self.updateJobProcessing(jobNo)

            if taskNo == 0:
                taskNo = self.createTask(sampleName,fileName,queue,jobNo)
            else:
                self.updateTaskProcessing(taskNo)
            ## Transaction End
            self.connection.commit()
            self.connection.autocommit = True
            return (True,jobNo,taskNo)
        except mysql.connector.Error as error:
            print("Failed to crate Log record to database rollback: {}".format(error))
            self.connection.rollback()
            return (False,jobNo,taskNo)
        except CompletedJobAndTaskError:
            print("Job or Task completed")
            return (False,jobNo,taskNo)

    def closeLog(self,jobNo,taskNo,status,projectName,jobName,resultPath, errorReason):

        try:
            self.connection.autocommit = False
            # Transaction Start
            self.updateJobCompleting(jobNo,status)
            self.updateTaskCompleting(taskNo,status,resultPath, errorReason)
            # Transaction End
            self.connection.commit()

            if self.checkJobCompletion(projectName,jobName):
                self.updateJobClosing(jobNo)

        except mysql.connector.Error as error:
            print("Failed to crate Log record to database rollback: {}".format(error))
            self.connection.rollback()
        finally:
            self.connection.autocommit = True

