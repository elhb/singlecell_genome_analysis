import sys
sys.stderr.write('module loaded\n')

class VcfEntry(object):
    """Class representing a vcf entry or genomic variant recorded in a vcf file
    
    Has the following attributes:
    self.chrom      -- holding the chromosome name
    self.pos        -- position on chromosome
    self.id         -- variation id
    self.refBase    -- the reference base
    self.altBases   -- listr with alternative bases
    self.varQual    -- variant quality
    self.passFilter -- if pass filter
    self.info       -- variant info field
    perSampleInfo   -- all the per sample infos, this is a dict with samplenames as keys and the following info:
	sampleNames -- key in dict, value is a dict with the following keys
	    DP -- sample read depth at variant position
	    GT -- sample genotype
	    AD -- sample allele distribution
	    GQ -- sample genotype quality
    """    

    def __init__(self,line,sampleNames):
	""" The init function takes a line from a vcf file and extracts the info into different suitable contatiners accesible by sample name keys
	"""

	#
	# import functions
	#
	import time
	import sys

	# set some initital values
	self.isIndel = None
	self.perSampleInfo = {sampleName:{} for sampleName in sampleNames}
	
	# save the original line for access later
	self.vcfLine = line
	
	# split the tab delimeted line into a list
	line = self.vcfLine.rstrip().split('\t')
	
	# go through the list by incrementing the position in the list (i)
	for i in range(len(line)):
	    
	    # if we are on first position in list we have the chromosome value
	    if i == 0:   self.chrom = line[i]
	    
	    # second position is position on chromosome
	    elif i == 1: self.pos = line[i]
	    
	    # third position is the variation ID
	    elif i == 2:
		self.id = line[i]
		# 3 following lines NOT USED:
		##if line[i] == '.':pass
		#    #AnalysisPipe.nonDbSNPVarCount+=1;
		#    #self .id = 'novel#'+str(AnalysisPipe.nonDbSNPVarCount)
	    
	    # postition 4 in list is reference base(s)
	    elif i == 3: self.refBase = line[i]
	    # postition 5 in list is alternate base(s)
	    elif i == 4: self.altBases = line[i].split(',') # split the alternate bases on "," so that we get each possible alternate "allele" in the list
	    
	    # position6 = variant quality (note not same as genotype quality which is sample specific)
	    elif i == 5: self.varQual = line[i]
	    
	    # on position 7 in vcf line we find the flag for passing filter os not
	    elif i == 6: self.passFilter = line[i]
	    
	    # then comes some info format is described in the vcf-header and we ignore it for now
	    elif i == 7: self.info = line[i]
	    
	    # here is the format of the per sample information stored in the vcf mostly this field equals "GT:AD:DP:GQ:PL"
	    elif i == 8: self.perSampleFormat = line[i]
	    
	    # from position ten and onwards we find the per sample informations formated as described in last field
	    elif i >= 9: # in other words for each of the fields after position ten do the following:
		
		# stor the per sample info in a new variable
		sampleInfo = line[i]
		
		# get the format and split it based on ":" so we now whats what in the per sample information
		infoTags = self.perSampleFormat.split(':')
		
		#  check that the first information is the genotype also check that the genotype is different from "./." which usually means that no reads are present
		if infoTags[0] == 'GT' and sampleInfo[0:2] != './.' or sampleInfo != './.':
		    
		    # split the per sample info in the same way as the format so we can pair them based on position in the list
		    sampleInfoList = sampleInfo.split(':') #GT:AD:DP:GQ:PL
		    #if sampleInfoList[0] == './.':continue
		    
		    # save all the info values accesible by info name/tag and sample name in a dictionary
		    for i2 in range(len(sampleInfoList)): self.perSampleInfo[sampleNames[i-9]][infoTags[i2]] = sampleInfoList[i2]
		    #    go through all pos inj list ----> the dictionary   [ key 1 is sampleName ] [ key 2 is name of the info field eg GT, DP etc] this equals the value in the same position of the per sample info list
		    
		    # if the DP tag/info is "." then set it to "0" so we can do some mathematics with it later
		    try:
			if self.perSampleInfo[sampleNames[i-9]]['DP'] == '.': self.perSampleInfo[sampleNames[i-9]]['DP'] = 0
			self.perSampleInfo[sampleNames[i-9]]['DP'] = int(self.perSampleInfo[sampleNames[i-9]]['DP'])
		    except KeyError: self.perSampleInfo[sampleNames[i-9]]['DP'] = 0
		    
		    # if the GQ tag/info is "." then set it to "0" so we can do some mathematics with it later
		    try:
			if self.perSampleInfo[sampleNames[i-9]]['GQ'] == '.': self.perSampleInfo[sampleNames[i-9]]['GQ'] = 0
			self.perSampleInfo[sampleNames[i-9]]['GQ'] = int(self.perSampleInfo[sampleNames[i-9]]['GQ'])
		    except KeyError: self.perSampleInfo[sampleNames[i-9]]['GQ'] = 0
		    
		    
		    # if the GT info is not there set it to "unknown"
		    if 'GT' not in self.perSampleInfo[sampleNames[i-9]]: self.perSampleInfo[sampleNames[i-9]]['GT'] = 'Unknown'
		    # same thing for AD
		    if 'AD' not in self.perSampleInfo[sampleNames[i-9]]: self.perSampleInfo[sampleNames[i-9]]['AD'] = 'Unknown'

		    # get angry if we cant find the DP tag and print an error message
		    try:
			if self.perSampleInfo[sampleNames[i-9]]['DP'] == '.': self.perSampleInfo[sampleNames[i-9]]['DP'] = 0
		    except KeyError:
			msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# VcfEntry format error: '+str(infoTags)+' '+str(sampleInfoList)+' '+str(self.perSampleInfo[sampleNames[i-9]])+' '+str(i)+' >'+str(sampleInfo)+'<'+'.\n'
			sys.stderr.write(msg)

	# now we have stored all the info we need finally we chech if the length of alt base(s) is the same as ref bases and set the type of the variation as indel or not
	for altBase in list(self.altBases):
	    if len(self.refBase) != len(altBase):
		self.type = 'INDEL'
		self.isIndel = True
	    elif not self.isIndel:
		self.type = 'SNV'
		self.isIndel = False

def vcfParser(infilename,subsetSize=None,logfile=False):
	"""This is a Generator taking an infilename ie vcf-file as input and yields VcfEntries Objects defined above"""
	
	#
	# import usefull stuff
	#
	import os
	import sys
	import time
	from misc import bufcount
	from misc import Progress
	
	#
	# set some initial values
	#
	nonDbSNPVarCount = 0
	vcfLineCount = 0
	toQcounter = 0
	lastChomosomeRead = None
	
	#
	# write a log message
	#
	if logfile: logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Loading vcf to memory, my id is '+str(os.getpid())+'.\n')
	
	
	#
	# check if we asked for a full file or a subset
	#
	if not subsetSize:
	    totalLineNum = bufcount(infilename)
	    if logfile:
		logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# '+str(totalLineNum)+' lines found in vcf file.\n')
	    	progress = Progress(totalLineNum, verb='full', logfile=logfile, unit='variants read from disk' ,mem=False, printint=1)
	else:
	    if logfile: progress = Progress(subsetSize, verb='full', logfile=logfile, unit='variants read from disk' ,mem=False)

	with open(infilename) as infile:
	    for line in infile:
		if not subsetSize and logfile: progress.update()
		if subsetSize and toQcounter >= subsetSize:
		    if logfile: logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Warning running on a subset of variants debbugging purpose.\n')
		    break

                # check the first letter in the line to know if it is the header line or not
		if line[0] == "#": #THIS LINE IS HEADER
                    if line[1] == 'C': #THIS LINE IS HEADER WITH SAMPLE INFO -- 
                        line = line.rstrip().split('\t') # split the  header line  with sample names in to a list
                        sampleNames = [line[i] for i in range(9,len(line))] # store only the sample names is a list ie skip "chr pos id ref alt ... etc"
                    continue # continue with, ie go to, next line

		# incremeant the line counter and create a vcfEnctry object as templated above in "class VcfEntry"
		vcfLineCount +=1
		variation = VcfEntry(line,sampleNames)

		# check the pass filter value of recently created variation object and if it is an indel and that the variation quality is above 30 otherwise trash the variation
		if variation.passFilter == 'PASS' and not variation.isIndel and variation.varQual >= 30:
		    
		    yield variation # if the variation pass the check above send it back to "the user"
		    
	# some log messages:
		    if subsetSize and logfile: progress.update()
		    toQcounter+=1
		    if variation.chrom != lastChomosomeRead and logfile: logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# started loading variations from '+str(variation.chrom)+'.\n')
		    lastChomosomeRead = variation.chrom
	if logfile: logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# vcf file loaded.\n')

def genotypeAsBases(variant,genotype):
    """ Use this function to translate a "0/1"-type genotype to a "A/T"-genotype """
    
    # get the alternative and refences bases into an array so that the bases can be accessed by integer indexes
    theBases = [variant.refBase]+variant.altBases
    
    # check that the genotype is not "no data" or "unknown"
    if genotype == 'Unknown' or genotype == './.':
        return genotype
    
    # if the genotype is integers sepparated by "/" translate the integers to corresponding base
    else:
        baseInts = genotype.split('/')
        return '/'.join([theBases[int(tmp_integer)] for tmp_integer in baseInts])
    
def parseVcfForADOrate(AnalysisPipe):

    import sys
    import operator
    import time
    from misc import sorted_nicely, percentage, thousandString

    patId = None
    donId = None
    samplesbyName = {}
    samplesNamesNoPD = []
    samplesbyId = {}
    refSamplesFirst = []
    samples = []
    refSamples = []
    for sample in AnalysisPipe.database.getSamples():
        samplesbyId[sample.id]=sample
        samplesbyName[sample.name]=sample
        samples.append(sample)
        if sample.refType:
            refSamplesFirst.append(sample)
            refSamples.append(sample)
        if sample.refType == 'Patient': patId = sample.name
        if sample.refType == 'Donor': donId = sample.name
        if not sample.refType: samplesNamesNoPD.append(sample.name)

    for sampleName in sorted_nicely(samplesNamesNoPD):
        sample = samplesbyName[sampleName]
        if sample not in refSamplesFirst: refSamplesFirst.append(sample)

    mainReferenceSample = 'no sample name here'
    if len(refSamples) == 1:
        mainReferenceSample = refSamples[0].name
    elif len(refSamples) > 1:
        if patId:
            AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# using patient sample ' +patId+ ' as reference sample for ADO rate estimation etc.\n')
            mainReferenceSample = patId
        else:
            mainReferenceSample = raw_input('cannot find which sample to use as reference in ADO rate estimation etc\nfollowing reference samples are around:\n'+'\n'.join([sample.name for sample in refSamples])+'\nplease give a Sample name to use as reference: ')
            AnalysisPipe.settings.setVariable('adoRefSample',mainReferenceSample)
            AnalysisPipe.settings.saveToDb()
    elif len(refSamples) == 0:
            mainReferenceSample = raw_input('cannot find any sample to use as reference in ADO rate estimation etc\nfollowing samples are around:\n'+'\n'.join(sorted_nicely(samplesbyName.keys()))+'\nplease give a Sample name to use as reference: ')
            AnalysisPipe.settings.setVariable('adoRefSample',mainReferenceSample)
            AnalysisPipe.settings.saveToDb()

    try: mainReferenceSample = samplesbyName[mainReferenceSample]
    except KeyError:
        msg = '#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# sample name ' +mainReferenceSample+ ' is not valid exiting.\n'
        sys.stderr.write(msg)
        AnalysisPipe.logfile.write(msg)
        sys.exit()

    referenceSamples = [mainReferenceSample]
    if patId and donId:
        referenceSamples.append(samplesbyName[donId])
    for mainReferenceSample in referenceSamples:
        for sample in refSamplesFirst:
            sample.onDiskVariantsList = open(sample.dataPath+'/adoVariantsList.ref='+mainReferenceSample.name+'.tsv','w')
            sample.onDiskVariantsSummary = open(sample.dataPath+'/adoVariantsSummary.ref='+mainReferenceSample.name+'.txt','w')
        coloredTableOut = open(AnalysisPipe.path+'/adoRate.colored.ref='+mainReferenceSample.name+'.txt','w')
        AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# using sample ' +mainReferenceSample.name+ ' as reference sample for ADO rate estimation etc.\n')
  
  
        identifcationPercentages = {}
        for sample in samples: identifcationPercentages[sample.name] = {'totalVariationCount':0,'correct':0,'dropout':0,'other':0}
  
        totalVariationCount = 0
        informativeVariationCount = 0
  
        AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Parsing variation per chromosome.\n')
        for chrom in sorted_nicely(AnalysisPipe.variations.keys()):#,variations in sorted(AnalysisPipe.variations.iteritems(), key=operator.itemgetter(0)):
            AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# now at '+str(chrom)+'.\n')
            for variation in AnalysisPipe.variations[ chrom ]:
    
                totalVariationCount += 1
    
                # patient and donor filters
                refDP = bool(variation.perSampleInfo[mainReferenceSample.name]['DP'] >= int(AnalysisPipe.settings.referenceSampleDPcutoff))
                refQual = bool(variation.perSampleInfo[mainReferenceSample.name]['GQ'] >= int(AnalysisPipe.settings.referenceSampleGQcutoff))
                refHom = bool(variation.perSampleInfo[mainReferenceSample.name]['GT'].split('/').count(variation.perSampleInfo[mainReferenceSample.name]['GT'].split('/')[0])==2)
    
                if refDP and not refHom and refQual:
                    informativeVariationCount += 1
                    identify = {variation.perSampleInfo[mainReferenceSample.name]['GT']:'\033[34m'+'correct'+'\033[0m',
                              './.':'\033[90m'+'NoData'+'\033[0m',
                              str(variation.perSampleInfo[mainReferenceSample.name]['GT'].split('/')[0])+'/'+variation.perSampleInfo[mainReferenceSample.name]['GT'].split('/')[0]:'\033[93m'+'dropout'+'\033[0m',
                              str(variation.perSampleInfo[mainReferenceSample.name]['GT'].split('/')[1])+'/'+variation.perSampleInfo[mainReferenceSample.name]['GT'].split('/')[1]:'\033[93m'+'dropout'+'\033[0m'}
                    identifyNoColor = {variation.perSampleInfo[mainReferenceSample.name]['GT']:'correct',
                              './.':'NoData',
                              str(variation.perSampleInfo[mainReferenceSample.name]['GT'].split('/')[0])+'/'+variation.perSampleInfo[mainReferenceSample.name]['GT'].split('/')[0]:'dropout',
                              str(variation.perSampleInfo[mainReferenceSample.name]['GT'].split('/')[1])+'/'+variation.perSampleInfo[mainReferenceSample.name]['GT'].split('/')[1]:'dropout'}
                    test=''
                    pfSamples = 0
                    for sample in refSamplesFirst:
                        sampleDP = None
                        sampleQual = None
      
                        # per sample filters
                        sampleDP   = bool(variation.perSampleInfo[sample.name]['DP'] >= int(AnalysisPipe.settings.sampleADOestDPcutoff))
                        sampleQual = bool(variation.perSampleInfo[sample.name]['GQ'] >= int(AnalysisPipe.settings.sampleADOestGQcutoff))
      
                        if sampleDP and sampleQual:
                            try:
                                test+=sample.name+'='+identify[variation.perSampleInfo[sample.name]['GT']]+'('+variation.perSampleInfo[sample.name]['GT']+')\t'
                                identifcationPercentages[sample.name][identifyNoColor[variation.perSampleInfo[sample.name]['GT']]] += 1
                                identifcationPercentages[sample.name]['totalVariationCount'] += 1
                                pfSamples+=1
                            except KeyError:
                                test+=sample.name+'='+'\033[31m''other('+variation.perSampleInfo[sample.name]['GT']+')\033[0m\t'
                                identifcationPercentages[sample.name]['other']
    
                            try:             sample.onDiskVariantsList.write(str(variation.chrom)+'\t'+str(variation.pos)+'\t'+str(variation.id)+'\t'+identifyNoColor[variation.perSampleInfo[sample.name]['GT']]+'\t'+variation.perSampleInfo[sample.name]['GT']+'\t'+str(variation.perSampleInfo[sample.name]['GQ'])+'\t'+str(variation.perSampleInfo[sample.name]['DP'])+'\t'+str(variation.perSampleInfo[sample.name]['AD'])+'\n')
                            except KeyError: sample.onDiskVariantsList.write(str(variation.chrom)+'\t'+str(variation.pos)+'\t'+str(variation.id)+'\t'+'other'                                                    +'\t'+variation.perSampleInfo[sample.name]['GT']+'\t'+str(variation.perSampleInfo[sample.name]['GQ'])+'\t'+str(variation.perSampleInfo[sample.name]['DP'])+'\t'+str(variation.perSampleInfo[sample.name]['AD'])+'\n')
      
                        else: test+=sample.name+'='+'\033[90m'+'LowDPorGQ'+'\033[0m'+'\t'
                    if pfSamples >= 3: coloredTableOut.write( str(variation.chrom)+'\t'+str(variation.pos)+'\t'+str(variation.id)+'\t'+variation.perSampleInfo[mainReferenceSample.name]['GT']+'\t'+test+'\n')
    
        coloredTableOut.write('\nSAMPLESUMS:')
        coloredTableOut.write('\n\t\t\t\t')
        for sample in refSamplesFirst:
            sampleIs = 'Unknown'
            sample.correctPercentage = percentage(identifcationPercentages[sample.name]['correct'],identifcationPercentages[sample.name]['totalVariationCount'])
            sample.adoPercentage = percentage(identifcationPercentages[sample.name]['dropout'],  identifcationPercentages[sample.name]['totalVariationCount'])
            #print sample.name
            #print sample.name,'-correct->',identifcationPercentages[sample.name]['correct']
            #print sample.name,'-dropout->',identifcationPercentages[sample.name]['dropout']
            #print sample.name,'---total->',identifcationPercentages[sample.name]['totalVariationCount']
  
            coloredTableOut.write(sample.name+'\t')
        coloredTableOut.write('\n\t\t\t\t')
        for sample in refSamplesFirst:
            coloredTableOut.write('correct%='+str(sample.correctPercentage)+''+'\t')
            sample.onDiskVariantsSummary.write('correct%='+str(sample.correctPercentage)+''+'\n')
        coloredTableOut.write('\n\t\t\t\t')
        for sample in refSamplesFirst:
            coloredTableOut.write('ADO%='+str(sample.adoPercentage)+''+'\t')
            sample.onDiskVariantsSummary.write('ADO%='+str(sample.adoPercentage)+''+'\n')
        coloredTableOut.write('\n\t\t\t\t')
        for sample in refSamplesFirst:
            coloredTableOut.write('other%='+str(percentage(identifcationPercentages[sample.name]['other'],identifcationPercentages[sample.name]['totalVariationCount']))+''+'\t')
            sample.onDiskVariantsSummary.write('other%='+str(percentage(identifcationPercentages[sample.name]['other'],identifcationPercentages[sample.name]['totalVariationCount']))+''+'\n')
        coloredTableOut.write('\n\t\t\t\t')
        for sample in refSamplesFirst:
            coloredTableOut.write('total='+str(identifcationPercentages[sample.name]['totalVariationCount'])+'\t')
            sample.onDiskVariantsSummary.write('total='+str(identifcationPercentages[sample.name]['totalVariationCount'])+'\n')
        coloredTableOut.write('\n')
  
        sys.stdout.write('ADOcheck: Found '+thousandString(str(informativeVariationCount))+' informative variant out of '+thousandString(str(totalVariationCount))+' parsed variants.\n')
        coloredTableOut.write('Found '+thousandString(str(informativeVariationCount))+' informative variant out of '+thousandString(str(totalVariationCount))+' parsed variants.\n')
        coloredTableOut.close()
    for sample in refSamplesFirst: sample.onDiskVariantsList.close()
    AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# All done.\n')

def readVcfToMem(AnalysisPipe):

    import os
    import sys
    import time
    
    AnalysisPipe.variations = {}
    AnalysisPipe.nonDbSNPVarCount = 0

    AnalysisPipe.vcfLineCount = 0
    if os.path.exists(AnalysisPipe.path+'/indels.recalibrated.vcf'): infilename = AnalysisPipe.path+'/indels.recalibrated.vcf'
    elif os.path.exists(AnalysisPipe.path+'/recalibrated_snps_raw_indels.vcf'): infilename = AnalysisPipe.path+'/recalibrated_snps_raw_indels.vcf'
    else:
        msg = '#ERROR_MSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(AnalysisPipe.masterPid)+'# cannot find the vfc file has the variant calling and recalibration finished without errors?\n'
        AnalysisPipe.logfile.write(msg)
        sys.stderr.write(msg)
        sys.exit()
    AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(AnalysisPipe.masterPid)+'# Loading vcf to memory.\n')
    lastChomosomeRead = None
    with open(infilename) as infile:
        for line in infile:
            #if sum([len(AnalysisPipe.variations[chrom]) for chrom in AnalysisPipe.variations]) >= 1e4:
            #    AnalysisPipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'# Warning running on a subset of variants debbugging purpose.\n')
            #    break

            if line[0] == "#": #THIS LINE IS HEADER
                if line[1] == 'C': #THIS LINE IS HEADER WITH SAMPLE INFO
                    line = line.rstrip().split('\t')
                    sampleNames = [line[i] for i in range(9,len(line))]
                continue

            AnalysisPipe.vcfLineCount +=1

            variation = VcfEntry(line,sampleNames)

            if variation.passFilter == 'PASS' and not variation.isIndel:
                try:             AnalysisPipe.variations[variation.chrom].append(variation)
                except KeyError: AnalysisPipe.variations[variation.chrom] = [variation]
                if variation.chrom != lastChomosomeRead: AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(AnalysisPipe.masterPid)+'# started loading variations from '+str(variation.chrom)+'.\n')
                lastChomosomeRead = variation.chrom
    AnalysisPipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(AnalysisPipe.masterPid)+'# vcf file loaded.\n')