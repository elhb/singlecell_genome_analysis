class Database(object):

    def __init__(self, dbPath):
        self.path = dbPath

    def getConnection(self,):
        #
        # Import useful stuff
        #
        import sqlite3
        import sys

        #
        # Create database and set
        #
        try: self.conn = sqlite3.connect(self.path)
        except sqlite3.OperationalError:
            print 'ERROR: Trouble with the database, please check your commandline.'
            sys.exit()
        self.c = self.conn.cursor()

    def commitAndClose(self,):
        #
        # commit changes and close connection
        #
        self.conn.commit()
        self.conn.close()

    def create(self,):
        """ creates the database holding all information used in the analysis """
        
        self.getConnection()
        
        #
        # Create tables
        #
        self.c.execute('''CREATE TABLE runs (startTime,command,commandLine,finishedSuccessfully,masterPid)''')
        self.c.execute('''CREATE TABLE fastqs (filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,sampleId,PRIMARY KEY (filePairId))''');
        self.c.execute('''CREATE TABLE settings (variableName,defaultValue,value,setTime,PRIMARY KEY (variableName))''')
        self.c.execute('''CREATE TABLE results (resultName,defaultValue,value,setTime,PRIMARY KEY (resultName))''')
        self.c.execute('''CREATE TABLE samples (sampleId,sampleName,refType,PRIMARY KEY (sampleId))''')
        
        self.commitAndClose()
        
        import os
        os.chmod(self.path, 0664)

    def addToRunsTable(self, startTime, command, commandLine, finishedSuccessfully, masterPid):
        
        self.getConnection()
        
        #
        # check if pid already in database
        #
        t = (masterPid,)
        data = self.c.execute('SELECT masterPid, startTime FROM runs WHERE masterPid=?',t).fetchall()        
        if data:
            for tmp1,tmp2 in data:
    
        #
        # if pid and startTime matches update the "finishedSuccessfully" entry
        #
            if tmp1 == masterPid and tmp2 == startTime:
                values = (startTime, command, commandLine, finishedSuccessfully, masterPid)
                self.c.execute('UPDATE runs SET finishedSuccessfully=? WHERE masterPid=? AND startTime=?', (finishedSuccessfully,masterPid,startTime))
        
        #
        # if not in the database add a new row
        #
        else:
            values = (startTime, command, commandLine, finishedSuccessfully, masterPid)
            self.c.execute('INSERT INTO runs VALUES (?,?,?,?,?)', values)
        
        self.commitAndClose()
        
        return 0

    def addSample(self, newSampleName,newSampleRefType=None):
        
        #
        # Imports
        #
        import sys
        import time
        
        #
        # open connection to database
        #
        self.getConnection()
        
        sampleNames = []
        sampleIds = []
        
        #
        # check if any of the fastqs already in database
        #
        data = self.c.execute('SELECT sampleId,sampleName,refType FROM samples').fetchall()
        if data:
            for (sampleId,sampleName,sampleRefType) in data:
            #sampleName = sampleName[0]
            sampleNames.append(sampleName)
            sampleIds.append(sampleId)
            if newSampleName in sampleNames:
            msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(AnalysisPipe.masterPid)+'# SampleName must be uniq, there is already a sample with name '+newSampleName+' , exiting.\n'
            AnalysisPipe.logfile.write(msg)
            sys.stderr.write(msg)
            sys.exit(1)
    
        
        if sampleIds:  sampleId = max(sampleIds)+1
        else:          sampleId = 0 
        AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(AnalysisPipe.masterPid)+'# Adding sample '+newSampleName+' to database with id '+str(sampleId)+'.\n')
        values = (sampleId,newSampleName,newSampleRefType)
        self.c.execute('INSERT INTO samples VALUES (?,?,?)', values)
        
        sample = Sample(sampleName=newSampleName, sampleId=sampleId,refType=newSampleRefType)
        sample.createDirs()
        
        self.commitAndClose()
        
        return 0

    def getSamples(self):
        #
        # Imports
        #
        import sys
        import time
        
        #
        # open connection to database
        #
        self.getConnection()
        
        refSamples = []
        samples = []
        
        data = self.c.execute('SELECT sampleId,sampleName,refType FROM samples').fetchall()
        if data:
            for (sampleId,sampleName,sampleRefType) in data:
            if sampleRefType: refSamples.append( Sample(sampleName=sampleName,sampleId=int(sampleId),refType=sampleRefType) )
            else:                samples.append( Sample(sampleName=sampleName,sampleId=int(sampleId),refType=None) )
        
        self.commitAndClose()
        
        return refSamples+samples

    def updateFastqReadCount(self,sample):
    
        self.getConnection()
    
        #
        # check if any of the fastqs already in database
        #
        data = self.c.execute('SELECT filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,sampleId FROM fastqs').fetchall()
        if data:
            for filePair in data:
            if int(sample.id) == int(filePair[-1]):
                filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,sampleId = filePair
                tmp = extractData(infile=sample.logPath+'/rubiconWgaTrimming.'+str(filePairId)+'.r1.log.txt',        pattern="Running wgaAdapterTrimmer.py\nProcessed a total of\t(?P<totalReads>\d+)\treads. \(.+\)\nProcessed a total of\t(?P<totalBases>\d+)\tbases \(.+\).\ntrimmed a total of\t(?P<trimmedBases>\d+)\tbases in the start of reads \(.+\).\nwgaAdapterTrimmer.py done exiting ...\n?")
                if type(tmp) != str:newreadcount = int(tmp['totalReads'])
                else: newreadcount = 'Unknown'
                self.c.execute('UPDATE fastqs SET readCount=? WHERE filePairId=?', (newreadcount,filePairId))
    
        self.commitAndClose()

    def addFastqs(self, sampleNameOrId, fastq1, fastq2):
    
        #
        # Imports
        #
        import sys
        import os
        import time
        
        fastq1 = os.path.abspath(fastq1)
        fastq2 = os.path.abspath(fastq2)
        
        samples = AnalysisPipe.database.getSamples()
        samplesbyName = {}
        samplesbyId = {}
        for sample in samples:
            samplesbyId[sample.id]=sample
            samplesbyName[sample.name]=sample
        sampleName = None
        sampleId = None
        try:
            if sampleNameOrId.isdigit() and (int(sampleNameOrId) in [sample.id for sample in samples]):
            sampleId = int(sampleNameOrId);
            sampleName = str(samplesbyId[int(sampleId)].name)
        except ValueError: pass
        if type(sampleId) == int and type(sampleName) == str: pass
        elif   sampleNameOrId  in samplesbyName.keys():
            sampleName = sampleNameOrId;
            sampleId = samplesbyName[sampleName].id
        else:
            msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(AnalysisPipe.masterPid)+'# SampleName (or id) must be registered in the database, there is no sample with name or id '+str(sampleNameOrId)+' ('+str(type(sampleNameOrId))+') , exiting.\n'
            AnalysisPipe.logfile.write(msg)
            sys.stderr.write(msg)
            sys.exit(1)
    
        #
        # open connection to database
        #
        self.getConnection()
        
        filePairId = None
        filePairIds = []
        
        #
        # check if any of the fastqs already in database
        #
        data = self.c.execute('SELECT filePairId,fastq1,fastq2 FROM fastqs').fetchall()
        if data:
            for filePair in data:
            filePairId = int(filePair[0])
            filePairIds.append(filePairId)
            for fastq in [fastq1, fastq2]:
                if fastq in filePair:
                message = 'ERROR: '+fastq+' already in the database.\nExiting after error.'
                print message
                AnalysisPipe.logfile.write(message+'\n')
                sys.exit(1)
        #
        # if not in the database add a new row
        #
        AnalysisPipe.logfile.write('Getting readcount for file'+fastq1+' ... \n')
        readCount = 'Unknown'#bufcount(fastq1)/4 #one read is four lines
        AnalysisPipe.logfile.write('...done. The file has '+str(readCount)+' reads.\n')
        addedToReadsTable = False#SEAseqPipeLine.startTimeStr
        minReadLength = 'NA'
    
        if filePairIds: filePairId = max(filePairIds)+1
        else: filePairId = 0
        values = (filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,sampleId)
        self.c.execute('INSERT INTO fastqs VALUES (?,?,?,?,?,?,?)', values)
        
        self.commitAndClose()
        
        return 0

    def getFastqs(self,):
        #
        # Imports
        #
        import sys
        
        #
        # open connection to database
        #
        self.getConnection()
            
        #
        # get att data in fastqs table
        #
        filePairs = self.c.execute('SELECT filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,sampleId FROM fastqs').fetchall()
        
        self.commitAndClose()
        
        #return [[readCount,fastq1,fastq2] if (not addedToReadsTable) else None for filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength in filePairs]
        return [[filePairId,readCount,fastq1,fastq2,sampleId] for filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,sampleId in filePairs]

    def getRuns(self, runTypes):
    
        self.getConnection()
    
        runsInfo = []
        data = self.c.execute('SELECT * FROM runs').fetchall()
        for startTime, command, commandLine, finishedSuccessfully, masterPid in data:
            if command in runTypes: runsInfo.append([startTime, command, commandLine, finishedSuccessfully, masterPid])
    
        self.commitAndClose()
    
        return runsInfo

class Settings(object,):
    
    def __init__(self, ):
        """ object holding the settings used for each part of the analysis """
	
        self.defaultValues = {
            'debug':False,
            'uppmaxProject':'b2014005',
            'parallelProcesses':16,
            'mode':'exome',
            'patientAndDonorDPcutoff':20,
            'patientAndDonorGQcutoff':30,
            'sampleIdentificationDPcutoff':1,
            'sampleIdentificationGQcutoff':1,
            'referenceSampleDPcutoff':20,
            'referenceSampleGQcutoff':30,
            'sampleADOestDPcutoff':20,
            'sampleADOestGQcutoff':30,
            'RDoverchromYscaleMax':0,
            'skiprubicon':0,
            'skipmalbac':0,
            'skipampli1':0
        }
        self.explenations = {
            'debug':'Flag for running the scripts in multiprocessing or as single process run [True/False] (default=False), Not functional',
            'uppmaxProject':'Project id used at uppmax for sbatch scripts [bXXXXXXX] (default=b2014005)',
            'parallelProcesses':'Number of process to run when doing multiprocess parts of analysis (defaul=16)',
            'mode':'Type of analysis either Whole genome seguencing (wgs) or exome sequencing (exome) [wgs/exome] (default=exome)',
            'patientAndDonorDPcutoff':'read depth cutoff for variation filtering during identification of informative variants to use when identifing cell origin (defaul=20)',
            'patientAndDonorGQcutoff':'genotype quality cutoff for variation filtering during identification of informative variants to use when identifing cell origin (defaul=30)',
            'sampleIdentificationDPcutoff':'read depth cutoff for variation in each sample filtering during identification of cell origin (defaul=1)',
            'sampleIdentificationGQcutoff':'genotype quality cutoff for variation in each sample filtering during identification of cell origin (defaul=1)',
            'referenceSampleDPcutoff':'read depth cutoff for variation filtering during identification of informative variants to use when estimating ADO (defaul=20)',
            'referenceSampleGQcutoff':'genotype quality cutoff for variation filtering during identification of informative variants to use when estimating ADO (defaul=30)',
            'sampleADOestDPcutoff':'read depth cutoff for variation in each sample filtering during ADO estimation (defaul=20)',
            'sampleADOestGQcutoff':'genotype quality cutoff for variation in each sample filtering during ADO estimation (defaul=30)',
            'RDoverchromYscaleMax':'the max value for the y scale of the RD over chrom graphs, 0 means automatic for each sample and chrom',
            'skiprubicon':'flag for skipping rubicon wga adapter trimming (default 0/False)',
            'skipmalbac':'flag for skipping malbac wga adapter trimming (default 0/False)',
            'skipampli1':'flag for skipping ampliOne wga adapter trimming (default 0/False)'
    
            }
        self.isDefault = {}
        self.setTime = {}
    
        self.debug = None
        self.uppmaxProject = None
        self.parallelProcesses = None
        self.mode = None
        self.patientAndDonorDPcutoff = None
        self.patientAndDonorGQcutoff = None
        self.sampleIdentificationDPcutoff = None
        self.sampleIdentificationGQcutoff = None
        self.referenceSampleDPcutoff = None
        self.referenceSampleGQcutoff = None
        self.sampleADOestDPcutoff = None
        self.sampleADOestGQcutoff = None
        self.RDoverchromYscaleMax = None
        self.skiprubicon = None
        self.skipmalbac = None
        self.skipampli1 = None
    
        self.setDefaults()
    
        def setDefaults(self,):
        for variableName, value in self.defaultValues.iteritems():
            self.__dict__[variableName] = value
            self.isDefault[variableName] = True
            self.setTime[variableName] = None
        return 0

    def loadFromDb(self,):
        
        import time
        
        #
        # Get the connection
        #
        AnalysisPipe.database.getConnection()
        
        #
        # Select data
        #
        gotData = False
        while not gotData:
            try:
            data = AnalysisPipe.database.c.execute('SELECT variableName,defaultValue,value,setTime FROM settings').fetchall()
            gotData = True
            except sqlite3.OperationalError:
            time.sleep(1)
        
        #
        # Parse data and add to object __dict__
        #
        if data:
            for variableName,default,value,setTime in data:
            self.__dict__[variableName]  = value
            self.isDefault[variableName] = default
            self.setTime[variableName]   = setTime
        
        #
        # close connection
        #
        AnalysisPipe.database.commitAndClose()
    
        def setVariable(self,variableName,value):
        import time
        assert variableName in self.explenations,'Error: you are trying to set an undefined variable.\n'
        self.__dict__[variableName]  = value
        self.isDefault[variableName] = False
        self.setTime[variableName]   = time.time()
        return 0

    def saveToDb(self,):
        
        #
        # imports
        #
        import time
        
        #
        # get connection
        #
        AnalysisPipe.database.getConnection()
        
            #
            # Look whats already in database, update it if older or default and set what is not
            #
        AnalysisPipe.logfile.write('checking whats in db.\n')
            alreadyInDb = {}
        data = AnalysisPipe.database.c.execute('SELECT variableName,defaultValue,value,setTime FROM settings').fetchall()
            if data:
                for variableName,default,value,setTime in data:
            #AnalysisPipe.logfile.write('processing variable '+variableName+'')
            alreadyInDb[variableName] = True
            
            if variableName in self.__dict__:
                if default and not self.isDefault[variableName] or setTime < self.setTime[variableName]:
                if type(self.__dict__[variableName]) in [dict,list]: self.__dict__[variableName] = str(self.__dict__[variableName])
                AnalysisPipe.logfile.write('processing variable '+variableName+''+', updating from '+str(value)+' to '+str(self.__dict__[variableName])+', old_setTime '+str(time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime(setTime)))+' new_setTime '+str(time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime(self.setTime[variableName])))+'.\n')
                AnalysisPipe.database.c.execute('UPDATE settings SET defaultValue=?, value=?, setTime=? WHERE variableName=?', (self.isDefault[variableName],self.__dict__[variableName],self.setTime[variableName],variableName))
                else: pass#AnalysisPipe.logfile.write(' no update needed.\n')
            
            #
            # Add new vars to database
            #
        AnalysisPipe.logfile.write('adding new vars to db:\n')
            for variableName in self.__dict__:
            if variableName in ['explenations','defaultValues','isDefault','setTime']:continue
            if variableName not in alreadyInDb:
            if type(self.__dict__[variableName]) in [dict,list]: self.__dict__[variableName] = str(self.__dict__[variableName])
            values = (variableName,self.isDefault[variableName],self.__dict__[variableName],self.setTime[variableName])
            AnalysisPipe.database.c.execute('INSERT INTO settings VALUES (?,?,?,?)', values)
            AnalysisPipe.logfile.write('variable '+variableName+' added to db with value '+str(self.__dict__[variableName])+',')
            if self.isDefault[variableName]:AnalysisPipe.logfile.write(' this is the default value.\n')
            else:AnalysisPipe.logfile.write(' non-default value.\n')
            else: pass#SEAseqPipeLine.logfile.write('variable\t'+variableName+'\talready in db.\n')
            
        AnalysisPipe.logfile.write('commiting changes to database.\n')
            AnalysisPipe.database.commitAndClose()
        
        return 0

