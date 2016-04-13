def extractData(infile=None,pattern=None,checkType='data'):
    import re
    if not pattern: return 'NoPatternGiven'
    if not infile: return 'NoInfileGiven'
    try:
	with open(infile) as infile:
	    data = infile.read()
	    p = re.compile(pattern)
	    m = p.search(data)
	    if m:
		if   checkType=='data':   return m.groupdict()
		elif checkType=='program':return 'FinishedOK'
	    else: return 'NoMatchFound'
    except IOError, e:
	assert e.errno == 2;
	return 'NoFileFound'

class Sample(object):

    def __init__(self, sampleName=None,sampleId=None,refType=None, analysispipe=None):
		
		self.analysispipe = analysispipe
		
		self.name = sampleName
		self.id = int(sampleId)
		self.refType = refType
		self.path = self.analysispipe.path+'/samples/'+self.name
		self.scriptPath = self.analysispipe.path+'/samples/'+self.name+'/script'
		self.dataPath   = self.analysispipe.path+'/samples/'+self.name+'/data'
		self.logPath    = self.analysispipe.path+'/samples/'+self.name+'/logs'
		self.plotsPath  = self.analysispipe.path+'/samples/'+self.name+'/plots'
		self.fastqcPath = self.analysispipe.path+'/samples/'+self.name+'/fastQC'
		self.dependencies = {}
		
		self.tempPath = "$SNIC_TMP"
		#self.tempPath = self.dataPath

    @property
    def readCount(self, ):
	tmpCounter = 0
	for filePairId,readCount,fastq1,fastq2,sampleId in self.analysispipe.database.getFastqs():
	    if int(sampleId) == self.id:
		try: tmpCounter+= readCount
		except TypeError: return 'Unknown'
	return tmpCounter

    def getFastqs(self,):
	self.fastqIds = []
	for filePairId,readCount,fastq1,fastq2,sampleId in self.analysispipe.database.getFastqs():
	    if int(sampleId) == self.id:
		self.fastqIds.append(filePairId)
		yield [filePairId,readCount,fastq1,fastq2,sampleId]

    def createDirs(self):
        import os
        try: os.makedirs(self.path)
        except OSError:pass
        try: os.makedirs(self.scriptPath)
        except OSError:pass
        try: os.makedirs(self.dataPath)
        except OSError:pass
        try: os.makedirs(self.fastqcPath)
        except OSError:pass
        try: os.makedirs(self.logPath)
        except OSError:pass
        try: os.makedirs(self.plotsPath)
        except OSError:pass

    def getStats(self):
	
	import re
	import time
	import sys
	import os
	import glob
	
	self.updateOrgReadCounts()
	
	stats = {}
	stats['orientations'] = self.getReadOrientationStats()
	stats['rubiconWgaTrimming'] = {}
	stats['malbacWgaTrimming']  = {}
	stats['illuminaAndNexteraTrimming']  = {}
	stats['ampliOneTrimming']  = {}
	stats['qualityTrimming']  = {}
	stats['removeEmptyReads']  = {}
	stats['bowtie2'] = {}
	
	for filePairId,readCount,fastq1,fastq2,sampleId in self.getFastqs():
	    stats['rubiconWgaTrimming'][filePairId] = {'r1':None,'r2':None}
	    stats['malbacWgaTrimming'][filePairId]  = {'r1':None,'r2':None}
	    stats['illuminaAndNexteraTrimming'][filePairId]  = {'r1':None,'r2':None}
	    stats['ampliOneTrimming'][filePairId]  = {'r1':None,'r2':None}
	    stats['qualityTrimming'][filePairId]    = {'r1':None,'r2':None}
	    
	    for read in ['r1','r2']:
		stats['rubiconWgaTrimming'][filePairId][read]         = extractData(infile=self.logPath+'/rubiconWgaTrimming.'+str(filePairId)+'.'+read+'.log.txt',        pattern="Running wgaAdapterTrimmer.py\nProcessed a total of\t(?P<totalReads>\d+)\treads. \(.+\)\nProcessed a total of\t(?P<totalBases>\d+)\tbases \(.+\).\ntrimmed a total of\t(?P<trimmedBases>\d+)\tbases in the start of reads \(.+\).\nwgaAdapterTrimmer.py done exiting ...\n?")
		stats['malbacWgaTrimming'][filePairId][read]          = extractData(infile=self.logPath+'/malbacWgaTrimming.'+str(filePairId)+'.'+read+'.log.txt',         pattern="(.+)?cutadapt .+\nCommand line parameters: -n 10 -g GTGAGTGATGGTTGAGGTAGTGTGGAG -a CTCCACACTACCTCAACCATCACTCAC .+\nMaximum error rate\: .+\%\n\s+No. of adapters\: 2\n\s+Processed reads\:\s+(?P<totalReads>\d+)\n\s+Processed bases\:\s+(?P<totalBases>\d+) bp \(.+ Mbp\)\n\s+Trimmed reads\:\s+(?P<trimmedReads>\d+) \(.+\%\)\n\s+Trimmed bases\:\s+(?P<trimmedBases>\d+) bp \(.+ Mbp\) \(.+\% of total\)\n\s+Too short reads\:\s+.+ \(.+\% of processed reads\)\n\s+Too long reads\:\s+.+ \(.+\% of processed reads\)\n\s+Total time\:\s+.+ s\n\s+Time per read\:\s+.+ ms")
		stats['illuminaAndNexteraTrimming'][filePairId][read] = extractData(infile=self.logPath+'/illuminaAndNexteraTrimming.'+str(filePairId)+'.'+read+'.log.txt',pattern="(.+)?cutadapt .+\nCommand line parameters: -n .+\nMaximum error rate\: .+\%\n\s+No. of adapters\: 4\n\s+Processed reads\:\s+(?P<totalReads>\d+)\n\s+Processed bases\:\s+(?P<totalBases>\d+) bp \(.+ Mbp\)\n\s+Trimmed reads\:\s+(?P<trimmedReads>\d+) \(.+\%\)\n\s+Trimmed bases\:\s+(?P<trimmedBases>\d+) bp \(.+ Mbp\) \(.+\% of total\)\n\s+Too short reads\:\s+.+ \(.+\% of processed reads\)\n\s+Too long reads\:\s+.+ \(.+\% of processed reads\)\n\s+Total time\:\s+.+ s\n\s+Time per read\:\s+.+ ms")
		stats['ampliOneTrimming'][filePairId][read]           = extractData(infile=self.logPath+'/ampli1WgaTrimming.'+str(filePairId)+'.'+read+'.log.txt',pattern="(.+)?cutadapt .+\nCommand line parameters: -n .+\nMaximum error rate\: .+\%\n\s+No. of adapters\: \d+\n\s+Processed reads\:\s+(?P<totalReads>\d+)\n\s+Processed bases\:\s+(?P<totalBases>\d+) bp \(.+ Mbp\)\n\s+Trimmed reads\:\s+(?P<trimmedReads>\d+) \(.+\%\)\n\s+Trimmed bases\:\s+(?P<trimmedBases>\d+) bp \(.+ Mbp\) \(.+\% of total\)\n\s+Too short reads\:\s+.+ \(.+\% of processed reads\)\n\s+Too long reads\:\s+.+ \(.+\% of processed reads\)\n\s+Total time\:\s+.+ s\n\s+Time per read\:\s+.+ ms")

		if stats['ampliOneTrimming'][filePairId][read] == 'NoMatchFound':
		    stats['ampliOneTrimming'][filePairId][read]      = extractData(infile=self.logPath+'/ampli1WgaTrimming.'+str(filePairId)+'.'+read+'.log.txt',         pattern="(.+)?cutadapt .+\nCommand line parameters: .+\nTrimming .+\nFinished in .+\n\n\=\=\= Summary \=\=\=\n\nTotal reads processed\:\s+(?P<totalReads>([\d\,]+))\nReads with adapters\:.+\nReads written \(passing filters\)\:.+\n\nTotal basepairs processed\:\s+(?P<totalBases>([\d\,]+)) bp\nTotal written \(filtered\)\:\s+(?P<writtenBases>([\d\,]+)) bp \(\d+(\.\d+)?\%\)\n")
		    try:
			stats['ampliOneTrimming'][filePairId][read]['totalReads']   = stats['ampliOneTrimming'][filePairId][read]['totalReads'].replace(',','')
			stats['ampliOneTrimming'][filePairId][read]['writtenBases'] = stats['ampliOneTrimming'][filePairId][read]['writtenBases'].replace(',','')
			stats['ampliOneTrimming'][filePairId][read]['totalBases']   = stats['ampliOneTrimming'][filePairId][read]['totalBases'].replace(',','')
			stats['ampliOneTrimming'][filePairId][read]['trimmedBases'] = int(stats['ampliOneTrimming'][filePairId][read]['totalBases']) - int(stats['ampliOneTrimming'][filePairId][read]['writtenBases'])
			#print 'match!!!!',stats['ampliOneTrimming'][filePairId][read],filePairId,read,self.id,self.logPath+'/ampli1WgaTrimming.'+str(filePairId)+'.'+read+'.log.txt'
		    except TypeError: pass#print 'PLATS1!!!!',stats['ampliOneTrimming'][filePairId][read],filePairId,read,self.id,self.logPath+'/ampli1WgaTrimming.'+str(filePairId)+'.'+read+'.log.txt'

		if stats['malbacWgaTrimming'][filePairId][read] == 'NoMatchFound':
		    stats['malbacWgaTrimming'][filePairId][read]      = extractData(infile=self.logPath+'/malbacWgaTrimming.'+str(filePairId)+'.'+read+'.log.txt',         pattern="(.+)?cutadapt .+\nCommand line parameters: .+\nTrimming .+\nFinished in .+\n\n\=\=\= Summary \=\=\=\n\nTotal reads processed\:\s+(?P<totalReads>([\d\,]+))\nReads with adapters\:.+\nReads written \(passing filters\)\:.+\n\nTotal basepairs processed\:\s+(?P<totalBases>([\d\,]+)) bp\nTotal written \(filtered\)\:\s+(?P<writtenBases>([\d\,]+)) bp \(\d+(\.\d+)?\%\)\n")
		    try:
			stats['malbacWgaTrimming'][filePairId][read]['totalReads']   = stats['malbacWgaTrimming'][filePairId][read]['totalReads'].replace(',','')
			stats['malbacWgaTrimming'][filePairId][read]['writtenBases'] = stats['malbacWgaTrimming'][filePairId][read]['writtenBases'].replace(',','')
			stats['malbacWgaTrimming'][filePairId][read]['totalBases']   = stats['malbacWgaTrimming'][filePairId][read]['totalBases'].replace(',','')
			stats['malbacWgaTrimming'][filePairId][read]['trimmedBases'] = int(stats['malbacWgaTrimming'][filePairId][read]['totalBases']) - int(stats['malbacWgaTrimming'][filePairId][read]['writtenBases'])
			#print 'match!!!!',stats['malbacWgaTrimming'][filePairId][read],filePairId,read,self.id
		    except TypeError: pass#print 'PLATS2!!!!',stats['malbacWgaTrimming'][filePairId][read],filePairId,read,self.id

		if stats['illuminaAndNexteraTrimming'][filePairId][read] == 'NoMatchFound':
		    stats['illuminaAndNexteraTrimming'][filePairId][read]      = extractData(infile=self.logPath+'/illuminaAndNexteraTrimming.'+str(filePairId)+'.'+read+'.log.txt',pattern="(.+)?cutadapt .+\nCommand line parameters: .+\nTrimming .+\nFinished in .+\n\n\=\=\= Summary \=\=\=\n\nTotal reads processed\:\s+(?P<totalReads>([\d\,]+))\nReads with adapters\:.+\nReads written \(passing filters\)\:.+\n\nTotal basepairs processed\:\s+(?P<totalBases>([\d\,]+)) bp\nTotal written \(filtered\)\:\s+(?P<writtenBases>([\d\,]+)) bp \(\d+(\.\d+)?\%\)\n")
		    try:
			stats['illuminaAndNexteraTrimming'][filePairId][read]['totalReads']   = stats['illuminaAndNexteraTrimming'][filePairId][read]['totalReads'].replace(',','')
			stats['illuminaAndNexteraTrimming'][filePairId][read]['writtenBases'] = stats['illuminaAndNexteraTrimming'][filePairId][read]['writtenBases'].replace(',','')
			stats['illuminaAndNexteraTrimming'][filePairId][read]['totalBases']   = stats['illuminaAndNexteraTrimming'][filePairId][read]['totalBases'].replace(',','')
			stats['illuminaAndNexteraTrimming'][filePairId][read]['trimmedBases'] = int(stats['illuminaAndNexteraTrimming'][filePairId][read]['totalBases']) - int(stats['illuminaAndNexteraTrimming'][filePairId][read]['writtenBases'])
			#print 'match!!!!',stats['illuminaAndNexteraTrimming'][filePairId][read],filePairId,read,self.id
		    except TypeError: pass#print 'PLATS3!!!!',stats['illuminaAndNexteraTrimming'][filePairId][read],filePairId,read,self.id

		stats['qualityTrimming'][filePairId][read]            = extractData(infile=self.logPath+'/qualityTrimming.'+str(filePairId)+'.'+read+'.log.txt',           pattern='(?P<totalBasess>\d+)\tbases\n(?P<trimmedBases>\d+)\ttrimmed')

	    stats['removeEmptyReads'][filePairId] = extractData(infile=self.logPath+'/removeEmptyReads.'+str(filePairId)+'.log.txt',pattern="""Running removeEmptyReads.py:\nHeader one is empty exiting.\n(?P<totalReads>\d+) read pairs processed.\n(?P<pairsOut>\d+) read pairs to outfiles .+.\n(?P<singlets>\d+) single reads to outfile .+.\nremoveEmptyReads Exiting.""")
	    stats['bowtie2'][filePairId]          = extractData(infile=self.logPath+'/stderr.bowtie2.'+str(filePairId)+'.txt',      pattern="""(?P<totalReads>\d+) reads; of these:\n\s+(?P<pairedReads>\d+) \(\d+.\d+\%\) were paired; of these:\n\s+(?P<notPropMapedPair>\d+) \(\d+.\d+\%\) aligned concordantly 0 times\n\s+(?P<properPairs>\d+) \(\d+.\d+\%\) aligned concordantly exactly 1 time\n\s+(?P<properPairsMultiMap>\d+) \(\d+.\d+\%\) aligned concordantly >1 times\n\s+----\n\s+(?P<notPropMapedPair2>\d+) pairs aligned concordantly 0 times; of these:\n\s+(?P<discordantPairs>\d+) \(\d+.\d+\%\) aligned discordantly 1 time\n\s+----\n\s+(?P<unMappedPair>\d+) pairs aligned 0 times concordantly or discordantly; of these:\n\s+(?P<possibleSingletons>\d+) mates make up the pairs; of these:\n\s+(?P<unMappedReads>\d+) \(\d+.\d+\%\) aligned 0 times\n\s+(?P<singleSingleMap>\d+) \(\d+.\d+\%\) aligned exactly 1 time\n\s+(?P<singleMultiMap>\d+) \(\d+.\d+\%\) aligned >1 times\n(?P<overallAlignmentRate>\d+.\d+)\% overall alignment rate""")

	stats['merging'] = extractData(infile=self.logPath+'/stderr.merging.'+self.name+'.txt',pattern="Finished reading inputs.+\n.+picard.sam.MergeSamFiles done. Elapsed time",checkType='program')
	pattern = 'LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n(?P<LIBRARY>.+)\t(?P<UNPAIRED_READS_EXAMINED>\d+)\t(?P<READ_PAIRS_EXAMINED>\d+)\t(?P<UNMAPPED_READS>\d+)\t(?P<UNPAIRED_READ_DUPLICATES>\d+)\t(?P<READ_PAIR_DUPLICATES>\d+)\t(?P<READ_PAIR_OPTICAL_DUPLICATES>\d+)\t(?P<PERCENT_DUPLICATION>\d+\,\d+)\t(?P<ESTIMATED_LIBRARY_SIZE>\d+)'
	oldpattern="""LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n(?P<Library>.+)\s+(?P<unPairedReads>\d+)\s+(?P<totalReads>\d+)\s+(?P<unMapped>\d+)\s+(?P<unPairedDups>\d+)\s+(?P<pairDups>\d+)\s+(?P<opticalDups>\d+)\s+(?P<percentageDuplication>\d+\,\d+)\s+(?P<estLibSize>\d+)"""
	stats['markDuplicatesMetrix'] = extractData(infile=self.logPath+'/markDuplicatesMetrix.'+self.name+'.txt',pattern=pattern)
	stats['fixedBamFlagstat']        = extractData(infile=self.logPath+'/fixedBamFlagstat.'+self.name+'.txt',       pattern="""(?P<totalReads>\d+) \+ 0 in total \(QC-passed reads \+ QC-failed reads\)\n(?P<duplicates>\d+) \+ 0 duplicates\n(?P<mapped>\d+) \+ 0 mapped \(\d+.\d+\%:-nan\%\)\n(?P<paired>\d+) \+ 0 paired in sequencing\n(?P<read1>\d+) \+ 0 read1\n(?P<read2>\d+) \+ 0 read2\n(?P<properlyPaired>\d+) \+ 0 properly paired \(\d+.\d+\%:-nan\%\)\n(?P<bothMapped>\d+) \+ 0 with itself and mate mapped\n(?P<singletons>\d+) \+ 0 singletons \(\d+.\d+\%:-nan\%\)\n(?P<mateOnDiffChr>\d+) \+ 0 with mate mapped to a different chr\n(?P<mateOnDiffChrq5>\d+) \+ 0 with mate mapped to a different chr \(mapQ>=5\)""")
	stats['reCalibratedBamFlagstat'] = extractData(infile=self.logPath+'/reCalibratedBamFlagstat.'+self.name+'.txt',pattern="""(?P<totalReads>\d+) \+ 0 in total \(QC-passed reads \+ QC-failed reads\)\n(?P<duplicates>\d+) \+ 0 duplicates\n(?P<mapped>\d+) \+ 0 mapped \(\d+.\d+\%:-nan\%\)\n(?P<paired>\d+) \+ 0 paired in sequencing\n(?P<read1>\d+) \+ 0 read1\n(?P<read2>\d+) \+ 0 read2\n(?P<properlyPaired>\d+) \+ 0 properly paired \(\d+.\d+\%:-nan\%\)\n(?P<bothMapped>\d+) \+ 0 with itself and mate mapped\n(?P<singletons>\d+) \+ 0 singletons \(\d+.\d+\%:-nan\%\)\n(?P<mateOnDiffChr>\d+) \+ 0 with mate mapped to a different chr\n(?P<mateOnDiffChrq5>\d+) \+ 0 with mate mapped to a different chr \(mapQ>=5\)""")
	stats['unmapRemovedBamFlagstat'] = extractData(infile=self.logPath+'/unmapRemovedBamFlagstat.'+self.name+'.txt',pattern="""(?P<totalReads>\d+) \+ 0 in total \(QC-passed reads \+ QC-failed reads\)\n(?P<duplicates>\d+) \+ 0 duplicates\n(?P<mapped>\d+) \+ 0 mapped \(\d+.\d+\%:-nan\%\)\n(?P<paired>\d+) \+ 0 paired in sequencing\n(?P<read1>\d+) \+ 0 read1\n(?P<read2>\d+) \+ 0 read2\n(?P<properlyPaired>\d+) \+ 0 properly paired \(\d+.\d+\%:-nan\%\)\n(?P<bothMapped>\d+) \+ 0 with itself and mate mapped\n(?P<singletons>\d+) \+ 0 singletons \(\d+.\d+\%:-nan\%\)\n(?P<mateOnDiffChr>\d+) \+ 0 with mate mapped to a different chr\n(?P<mateOnDiffChrq5>\d+) \+ 0 with mate mapped to a different chr \(mapQ>=5\)""")
	stats['qualFilteredBamFlagstat'] = extractData(infile=self.logPath+'/qualFilteredBamFlagstat.'+self.name+'.txt',pattern="""(?P<totalReads>\d+) \+ 0 in total \(QC-passed reads \+ QC-failed reads\)\n(?P<duplicates>\d+) \+ 0 duplicates\n(?P<mapped>\d+) \+ 0 mapped \(\d+.\d+\%:-nan\%\)\n(?P<paired>\d+) \+ 0 paired in sequencing\n(?P<read1>\d+) \+ 0 read1\n(?P<read2>\d+) \+ 0 read2\n(?P<properlyPaired>\d+) \+ 0 properly paired \(\d+.\d+\%:-nan\%\)\n(?P<bothMapped>\d+) \+ 0 with itself and mate mapped\n(?P<singletons>\d+) \+ 0 singletons \(\d+.\d+\%:-nan\%\)\n(?P<mateOnDiffChr>\d+) \+ 0 with mate mapped to a different chr\n(?P<mateOnDiffChrq5>\d+) \+ 0 with mate mapped to a different chr \(mapQ>=5\)""")
	stats['noDuplicatesBamFlagstat'] = extractData(infile=self.logPath+'/noDuplicatesBamFlagstat.'+self.name+'.txt',pattern="""(?P<totalReads>\d+) \+ 0 in total \(QC-passed reads \+ QC-failed reads\)\n(?P<duplicates>\d+) \+ 0 duplicates\n(?P<mapped>\d+) \+ 0 mapped \(\d+.\d+\%:-nan\%\)\n(?P<paired>\d+) \+ 0 paired in sequencing\n(?P<read1>\d+) \+ 0 read1\n(?P<read2>\d+) \+ 0 read2\n(?P<properlyPaired>\d+) \+ 0 properly paired \(\d+.\d+\%:-nan\%\)\n(?P<bothMapped>\d+) \+ 0 with itself and mate mapped\n(?P<singletons>\d+) \+ 0 singletons \(\d+.\d+\%:-nan\%\)\n(?P<mateOnDiffChr>\d+) \+ 0 with mate mapped to a different chr\n(?P<mateOnDiffChrq5>\d+) \+ 0 with mate mapped to a different chr \(mapQ>=5\)""")
	
	#self.logPath+'/'+self.name+'.qacompute.stdout.txt'
	#self.logPath+'/'+self.name+'.qacompute.stderr.txt'
	#self.dataPath+'/'+self.name+'.qacompute.out '
	#self.logPath+'/'+self.name+'.stderr.caluclateHsmetrics.txt'
	pattern = 'READ_GROUP\n(?P<BAIT_SET>.+)\t(?P<GENOME_SIZE>\d+)\t(?P<BAIT_TERRITORY>\d+)\t(?P<TARGET_TERRITORY>\d+)\t(?P<BAIT_DESIGN_EFFICIENCY>\d+(\,\d+)?)\t(?P<TOTAL_READS>\d+)\t(?P<PF_READS>\d+)\t(?P<PF_UNIQUE_READS>\d+)\t(?P<PCT_PF_READS>\d+(\,\d+)?)\t(?P<PCT_PF_UQ_READS>\d+(\,\d+)?)\t(?P<PF_UQ_READS_ALIGNED>\d+)\t(?P<PCT_PF_UQ_READS_ALIGNED>\d+(\,\d+)?)\t(?P<PF_UQ_BASES_ALIGNED>\d+)\t(?P<ON_BAIT_BASES>\d+)\t(?P<NEAR_BAIT_BASES>\d+)\t(?P<OFF_BAIT_BASES>\d+)\t(?P<ON_TARGET_BASES>\d+)\t(?P<PCT_SELECTED_BASES>\d+(\,\d+)?)\t(?P<PCT_OFF_BAIT>\d+(\,\d+)?)\t(?P<ON_BAIT_VS_SELECTED>\d+(\,\d+)?)\t(?P<MEAN_BAIT_COVERAGE>\d+(\,\d+)?)\t(?P<MEAN_TARGET_COVERAGE>\d+(\,\d+)?)\t(?P<PCT_USABLE_BASES_ON_BAIT>\d+(\,\d+)?)\t(?P<PCT_USABLE_BASES_ON_TARGET>\d+(\,\d+)?)\t(?P<FOLD_ENRICHMENT>\d+(\,\d+)?)\t(?P<ZERO_CVG_TARGETS_PCT>\d+(\,\d+)?)\t(?P<FOLD_80_BASE_PENALTY>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_2X>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_10X>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_20X>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_30X>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_40X>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_50X>(\?)|(\d+?(\,\d+)?))\t(?P<PCT_TARGET_BASES_100X>\d+(\,\d+)?)\t(?P<HS_LIBRARY_SIZE>(\s?)|(\d+(\,\d+)?))\t(?P<HS_PENALTY_10X>\d+(\,\d+)?)\t(?P<HS_PENALTY_20X>\d+(\,\d+)?)\t(?P<HS_PENALTY_30X>\d+(\,\d+)?)\t(?P<HS_PENALTY_40X>\d+(\,\d+)?)\t(?P<HS_PENALTY_50X>\d+(\,\d+)?)\t(?P<HS_PENALTY_100X>\d+(\,\d+)?)'
	stats['hs_metrics.summary'] = extractData(infile=self.dataPath+'/'+self.name+'.hs_metrics.summary.txt',pattern=pattern)

	# make sums
	if not list(self.getFastqs()):
	    self.analysispipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# No fastq files found for sample: '+self.name+' continuing with next sample.\n')
	    sys.stderr.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# No fastq files found for sample: '+self.name+' continuing with next sample.\n')
	    self.stats = stats
	    return
	for program in ['illuminaAndNexteraTrimming','malbacWgaTrimming','qualityTrimming','rubiconWgaTrimming','ampliOneTrimming']:
	    try:
    		sums = {variable:0 for variable in stats[program][self.getFastqs().next()[0]]['r1'].keys()}
		for filePairId,readCount,fastq1,fastq2,sampleId in self.getFastqs():
		    for read in ['r1','r2']:
			for variable, value in stats[program][filePairId][read].iteritems(): sums[variable]+=float(value)
		stats[program]['sum']= sums
	    except AttributeError: pass
	for program in ['removeEmptyReads','bowtie2']:
	    try:
		sums = {variable:0 for variable in stats[program][self.getFastqs().next()[0]].keys()}
		for filePairId,readCount,fastq1,fastq2,sampleId in self.getFastqs():
		    for variable, value in stats[program][filePairId].iteritems(): sums[variable]+=float(value)
		stats[program]['sum']= sums
	    except AttributeError: pass

	# classification
	if os.path.exists(self.dataPath+'/sampleClassification.txt'):
	    data = open(self.dataPath+'/sampleClassification.txt').read().split('\n')
	    stats['sampleIs'] = data[0]
	    if data[0] == 'mix' or data[0] =='Unknown': refsample = None
	    else: refsample = data[1]
	    if 'LowConf' in data: stats['sampleIsLowConf'] = True
	    else:stats['sampleIsLowConf'] = False
	    data = open(self.dataPath+'/identificationVariantsSummary.txt').read().split('\n')
	    stats['pat%'] = data[0].split('=')[-1]
	    stats['don%'] = data[1].split('=')[-1]
	    stats['mix%'] = data[2].split('=')[-1]
	    stats['other%'] = data[3].split('=')[-1]
	    stats['totalIdVar'] = data[4].split('=')[-1]
	else:
	    stats['sampleIsLowConf'] = 'NA'
	    stats['sampleIs'] = 'NA'
	    stats['pat%'] = 'NA'
	    stats['don%'] = 'NA'
	    stats['mix%'] = 'NA'
	    stats['other%'] = 'NA'
	    stats['totalIdVar'] = 'NA'

	# ado summary:
	#self.dataPath+'/adoVariantsSummary.ref='+mainReferenceSample.name+'.txt'
	potentialInFiles = list(glob.iglob( self.dataPath+'/adoVariantsSummary.ref=*.txt' ))
	eraseAdo = False
	if len(potentialInFiles) >= 1:
	    if len(potentialInFiles) == 1: infile = potentialInFiles[0]
	    if len(potentialInFiles) > 1:
		if refsample: infile = self.dataPath+'/adoVariantsSummary.ref='+refsample+'.txt'
		else:  infile = potentialInFiles[0]; eraseAdo = True
	    with open(infile) as infile:
		stats['refsample'] = infile.name.split('=')[-1][:-4]
		data = infile.read().split('\n')
		stats['correct%'] = data[0].split('=')[-1]
		stats['dropout%'] = data[1].split('=')[-1]
		stats['other%'] = data[2].split('=')[-1]
		stats['totalAdoVar'] = data[3].split('=')[-1]
	else:
		stats['refsample']= 'NA'
		stats['correct%'] = 'NA'
		stats['dropout%'] = 'NA'
		stats['other%'] = 'NA'
		stats['totalAdoVar'] = 'NA'
	if eraseAdo:
		stats['refsample']= 'NA'
		stats['correct%'] = 'NA'
		stats['dropout%'] = 'NA'
		stats['other%'] = 'NA'
		stats['totalAdoVar'] = 'NA'

	# target coverage
	if os.path.exists(self.dataPath+'/targetCoveragestat.txt'):
	    with open(self.dataPath+'/targetCoveragestat.txt') as data: stats['averageTartgetCoverage'] = data.read().rstrip().split(' ')[-2]
	else: stats['averageTartgetCoverage'] = 'NA'
	    
	self.stats = stats

	# debug output
	#print "\n######## "+self.name+" ######## "
	#for key,value in stats.iteritems():
	#    print '    ',key
	#    try:
	#	for key2,value2 in value.iteritems():
	#	    assert type(value2) == dict
	#	    print '        ',key2,value2
	#    except:print '        ',value
	
	return 0

    def updateOrgReadCounts(self,):
	self.analysispipe.database.updateFastqReadCount(self)

    def getReadOrientationStats(self,):

	import time
	import operator
	import shutil
	import os
	import sys
	import pysam
	import time
	import gzip
	orientations = {'PP':0,'FR':0,'FRsp':0,'FF':0,'RF':0,'RR':0,'difChrom':0}
	
	if not os.path.exists(self.dataPath+'/orientations.pylist.gz'):	
	
	    try: uppmax_temp = os.environ["SNIC_TMP"]
	    except:
		uppmax_temp = None
		print 'Not on uppmax no temporary directory'
	
	    if not os.path.exists(self.dataPath+'/'+self.name+'.noDuplicates.bam'): self.analysispipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# Skipping oriantations stats for sample '+self.name+' the infile has not been created yet...\n'); return orientations
	
	    if uppmax_temp:
		try:os.mkdir(uppmax_temp+'/fnuttglugg_TMP')
		except OSError:pass
		if not os.path.exists(uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bam'):
		    self.analysispipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# Copying '+self.dataPath+'/'+self.name+'.noDuplicates.bam'+' to temp location for faster reading from disk, '+uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bam'+' \n')
		    shutil.copy(self.dataPath+'/'+self.name+'.noDuplicates.bam',uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bam')
		else:
		    print 'WARNING: rerun of '+uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bam'+' skipping copy!!'
		    self.analysispipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# WARNING: rerun of '+uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bam'+' skipping copy!!\n')
		bamfileName  = uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bam'
	    else:bamfileName = self.dataPath+'/'+self.name+'.noDuplicates.bam'
	    if uppmax_temp:
		try:os.mkdir(uppmax_temp+'/fnuttglugg_TMP')
		except OSError:pass
		if not os.path.exists(uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bai'):
		    self.analysispipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# Copying '+self.dataPath+'/'+self.name+'.noDuplicates.bai'+' to temp location for faster reading from disk, '+uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bai'+' \n')
		    try: shutil.copy(self.dataPath+'/'+self.name+'.noDuplicates.bai',uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bai')
		    except IOError as e: pass
		else:
		    print 'WARNING: rerun of '+uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bai'+' skipping copy!!'
		    self.analysispipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# WARNING: rerun of '+uppmax_temp+'/fnuttglugg_TMP'+'/'+self.name+'.noDuplicates.bai'+' skipping copy!!\n')
	    
    #	bamfileName = 	self.dataPath+'/'+self.name+'.noDuplicates.bam'
	    try: bamfile = pysam.Samfile(bamfileName, "rb")
	    except IOError:
		self.analysispipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# Skipping oriantation stats for sample '+self.name+' the infile has not been created yet...\n')
		return orientations
	    except ValueError:
		self.analysispipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# Skipping oriantation stats for sample '+self.name+' the infile is not finished for processing...\n')
		return orientations
	
	    self.analysispipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# Loading oriantation stats for sample '+self.name+'...\n')
	    try:
		for read in bamfile.fetch():
    
		    orientation = None
		    if read.is_paired and not read.is_unmapped and not read.mate_is_unmapped:
			if read.tid == read.rnext:
			    if read.is_proper_pair:
				orientation = 'PP'
			    elif   read.pos < read.pnext:
				if read.is_reverse: orientation = 'R'
				else: orientation = 'F'
				if read.mate_is_reverse: orientation += 'R'
				else: orientation += 'F'
			    elif read.pnext < read.pos:
				if read.mate_is_reverse: orientation = 'R'
				else: orientation = 'F'
				if read.is_reverse: orientation += 'R'
				else: orientation += 'F'
			    elif   read.pos == read.pnext: orientation = 'FRsp'
			else: orientation = 'difChrom'
			orientations[orientation]+=1
		total = sum(orientations.values())
    
		output = self.name+'\n'
		for thingy,plingy in orientations.iteritems():
		    output+= thingy+' '+str(percentage(plingy,total))+'\n'
		print output
    
	    except ValueError as e:
		if e == 'fetch called on bamfile without index':
		    self.analysispipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# Skipping insert size plot for sample '+self.name+' the bam index is not present...\n')
		    sys.stderr.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# Skipping insert size plot for sample '+self.name+' the bam index is not present...\n')
		    return orientations

	    orientationsOnDisk = gzip.open(self.dataPath+'/orientations.pylist.gz','wb',9)
	    orientationsOnDisk.write(str(orientations))
	    orientationsOnDisk.close()
	else:
	    self.analysispipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# Loading orientations from orintations.pylist-file '+self.name+' ...\n')
	    orientations = eval(gzip.open(self.dataPath+'/orientations.pylist.gz','rb').read())

	self.analysispipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(self.analysispipe.masterPid)+'# Loaded oriantation stats for sample '+self.name+' joining to main...\n')
	return orientations

