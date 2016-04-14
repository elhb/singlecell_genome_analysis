def generate_index(analysispipe):
    
    outputStr = """
    <!DOCTYPE HTML>
        <html lang="en-US">
            <head>
                <meta charset="UTF-8">
                <script src="http://d3js.org/d3.v3.js"></script>
                <link rel="stylesheet" href="static/css/style.css">
            </head>
            <body>
            <script src="https://rawgit.com/gka/d3-jetpack/master/d3-jetpack.js"></script>
            <script src="/static/js/all_sample_summary_table.js"></script>
            <div id="all_sample_summary_table"></div>"""
    outputStr += """
            </body>
        </html>"""
        
    return outputStr

def generate_oldStyle_index(analysispipe):
    import time
    import operator
    import os
    import glob
    import multiprocessing
    import pysam
    from misc import thousandString, percentage, sorted_nicely
    analysispipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(analysispipe.masterPid)+'# Creating a report located at '+analysispipe.path+'/report.htm'+' \n')
  
    analysispipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(analysispipe.masterPid)+'# Getting stats for all samlpes ....'+' \n')
    #samples = analysispipe.database.getSamples()
    samples = []
    sampleCount = sum(1 for sample in analysispipe.database.getSamples())

    poolOfProcesses = multiprocessing.Pool(min([multiprocessing.cpu_count(),sampleCount]),maxtasksperchild=1)
    parallelResults = poolOfProcesses.imap(reportParallel,analysispipe.database.getSamples(),chunksize=1)

    samplesbyName = {}
    samplesbyId = {}
    refSamplesFirst = []
    lengthOfChromByName = {}
    #for sample in parallelResults:
    for sample in analysispipe.database.getSamples():
        analysispipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(analysispipe.masterPid)+'# Loading stats for '+sample.name+' ....'+' \n')
        samplesbyId[sample.id]=sample
        samplesbyName[sample.name]=sample
        samples.append(sample)
        sample.getStats()
        if sample.refType: refSamplesFirst.append(sample)
    ### get chromosome sizes
        if not os.path.exists(sample.dataPath+'/'+sample.name+'.noDuplicates.bam'):
          analysispipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(analysispipe.masterPid)+'# Cannot fetch chromsizes from bam for sample '+sample.name+' the infile has not been created yet...\n')
        bamfileName = sample.dataPath+'/'+sample.name+'.noDuplicates.bam'
        try: bamfile = pysam.Samfile(bamfileName, "rb")
        except IOError:
          analysispipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(analysispipe.masterPid)+'# Cannot fetch chromsizes from bam for sample '+sample.name+' the infile has not been created yet...\n')
          continue
        except ValueError:
          analysispipe.logfile.write('#WARNING#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(analysispipe.masterPid)+'# Cannot fetch chromsizes from bam for sample '+sample.name+' the infile is not finished for processing...\n')
          continue
        for chrom in bamfile.header['SQ']: lengthOfChromByName[chrom['SN']] = chrom['LN']
    ### done
    for sampleId, sample in sorted(samplesbyId.iteritems(), key=operator.itemgetter(0)):
        if sample not in refSamplesFirst: refSamplesFirst.append(sample)
    analysispipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(analysispipe.masterPid)+'# Sample stats loaded, creating output ....'+' \n')
    
    outputStr = ''
    #reportFile = open(analysispipe.path+'/report.htm','w',1)
    outputStr += '<html>'
    outputStr += """<head><style>
        body {
          background-color:white
          font-family : Arial,"Myriad Web",Verdana,Helvetica,sans-serif;
        }
        h1   {color:black}
        p    {color:green}
        
        table, th, td {
          border: 1px solid black;
          border-collapse: collapse;
          
          font-size : 12;
        }
        th {
          text-align: center;
          background-color:darkgray;
          color:white;
          border: 1px solid black;
        }
        table {border-spacing: 1px;}
        th,td {padding: 5px;}
        td {text-align: center;}
    </style></head>"""

    outputStr += '<body>'
    outputStr += '<h1>Analysis Report '+analysispipe.path+'</h1>\n'

    outputStr += '<h2>List of samples:</h2>'
    outputStr += '<table>'
    outputStr += '<tr>'
    outputStr += '<th>Sample Id</th>'
    outputStr += '<th>Sample Name</th>'
    outputStr += '<th>Read Pair Count</th>'
    outputStr += '<th>Pairs After Trimming</th>'
    outputStr += '<th>Overall Mapping Rate<br>(bowtie)</th>'
    outputStr += '<th>Overall Mapping Rate<br>(flagstat)</th>'
    outputStr += '<th>Duplication Rate<br>(preFilter)</th>'
    outputStr += '<th>Unmapped or<br>Low Mapping Quality</th>'
    outputStr += '<th>Left After<br>All Filtering</th>'
    outputStr += '<th>'+analysispipe.settings.mode+'<br>Coverage</th>'
    outputStr += '</tr>'
    for sample in refSamplesFirst:
        outputStr += '<tr>'
        outputStr += '<td>'
        outputStr += str(sample.id)
        outputStr += '</td>'
        outputStr += '<td>'+sample.name+'</td>'
        outputStr += '<td>'+thousandString(str(sample.readCount))+'</td>'
        try:            outputStr += '<td>'+thousandString(str(int(sample.stats['removeEmptyReads']['sum']['pairsOut'])))+' ('+str(percentage(sample.stats['removeEmptyReads']['sum']['pairsOut'],sample.readCount))+'%)</td>'
        except KeyError:outputStr += '<td><font color="red">NA (NA%)</font></td>'
        try:
          tmpCount  = sum([sample.stats['bowtie2']['sum'][tmp] for tmp in ['singleSingleMap','singleMultiMap']])
          tmpCount += sum([sample.stats['bowtie2']['sum'][tmp] for tmp in ['discordantPairs','properPairsMultiMap','properPairs']])*2
          assert sample.stats['bowtie2']['sum']['totalReads'] == sample.stats['removeEmptyReads']['sum']['pairsOut'], 'Error: pair counts do not match'
          outputStr += '<td>'+str(percentage(tmpCount,sample.stats['bowtie2']['sum']['totalReads']*2))+'%</td>'
        except KeyError:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['fixedBamFlagstat']['mapped']),int(sample.stats['fixedBamFlagstat']['totalReads'])))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['fixedBamFlagstat']['duplicates']),int(sample.stats['fixedBamFlagstat']['totalReads'])))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['fixedBamFlagstat']['totalReads'])-int(sample.stats['qualFilteredBamFlagstat']['totalReads']),int(sample.stats['fixedBamFlagstat']['totalReads'])))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['noDuplicatesBamFlagstat']['mapped']),int(sample.stats['fixedBamFlagstat']['totalReads'])))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        if 'averageTartgetCoverage' in sample.stats and sample.stats['averageTartgetCoverage'] != 'NA':outputStr += '<td>'+str(sample.stats['averageTartgetCoverage'])+'%</td>'
        else: outputStr += '<td><font color="red">NA%</font></td>'
        outputStr += '<tr>'
    outputStr += '</table>'

    outputStr += '<h2>List of fastqs:</h2>'
    outputStr += '<table>'
    outputStr += '<tr>'
    outputStr += '<th>File Id</th>'
    outputStr += '<th># Read Pairs</th> '
    outputStr += '<th>Sample Name</th>'
    outputStr += '<th>File Name (r1)</th>'
    outputStr += '<th>FastQC r1</th>'
    outputStr += '<th>FastQC r2</th>'
    outputStr += '</tr>'
    for filePairId,readCount,fastq1,fastq2,sampleId in sorted(analysispipe.database.getFastqs(), key=operator.itemgetter(0)):
        outputStr += '<tr>'
        outputStr += '<td>'+str(filePairId)+'</td>'
        outputStr += '<td>'+thousandString(str(readCount))+'</td>'
        try:               outputStr += '<td>'+samplesbyId[int(sampleId)].name+'</td>'
        except KeyError:   outputStr += '<td>'+'Unknown'                      +'</td>'
        outputStr += '<td>'+fastq1+'</td>'
        try:
          if os.path.exists(samplesbyId[int(sampleId)].fastqcPath+'/'+str(filePairId)+'.r1.allTrimmed_fastqc.html'): outputStr += '<td><a href="'+'static/runData'+'/per_sample_info/'+samplesbyId[int(sampleId)].name+'/fastQC'+'/'+str(filePairId)+'.r1.allTrimmed_fastqc.html'+'">here</a></td>'
          else:outputStr += '<td><font color="red">NA</font></td>'
          if os.path.exists(samplesbyId[int(sampleId)].fastqcPath+'/'+str(filePairId)+'.r2.allTrimmed_fastqc.html'): outputStr += '<td><a href="'+'static/runData'+'/per_sample_info/'+samplesbyId[int(sampleId)].name+'/fastQC'+'/'+str(filePairId)+'.r2.allTrimmed_fastqc.html'+'">here</a></td>'
          else:outputStr += '<td><font color="red">NA</font></td>'
        except KeyError:
          outputStr += '<td><font color="red">NA</font></td>'
          outputStr += '<td><font color="red">NA</font></td>'
        outputStr += '</tr>'
    outputStr += '</table>'

    outputStr += '<h2>Trimming Details:</h2>'
    outputStr += '<h3>Samples:</h3>'
    outputStr += '<table>'
    outputStr += '<tr>'
    outputStr += '<th>Sample Id</th>'
    outputStr += '<th>Sample Name</th>'
    outputStr += '<th>Original Bases</th>'
    outputStr += '<th>% rubicon adapter</th>'
    outputStr += '<th>% malbac adapter</th>'
    outputStr += '<th>% ampliOne adapter</th>'
    outputStr += '<th>% illumina Adapters</th>'
    outputStr += '<th>% quality trimmed</th>'
    outputStr += '</tr>'
    for sample in refSamplesFirst:
        outputStr += '<tr>'
        outputStr += '<td>'+str(sample.id)+'</td>'
        outputStr += '<td>'+sample.name+'</td>'
        try: outputStr += '<td>'+thousandString(str(int(sample.stats['rubiconWgaTrimming']['sum']['totalBases'])))+'</td>'
        except (KeyError, TypeError) as e: outputStr += '<td><font color="red">NA</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['rubiconWgaTrimming']['sum']['trimmedBases']),int(sample.stats['rubiconWgaTrimming']['sum']['totalBases'])))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['malbacWgaTrimming']['sum']['trimmedBases']),int(sample.stats['rubiconWgaTrimming']['sum']['totalBases'])))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['ampliOneTrimming']['sum']['trimmedBases']),int(sample.stats['rubiconWgaTrimming']['sum']['totalBases'])))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['illuminaAndNexteraTrimming']['sum']['trimmedBases']),int(sample.stats['rubiconWgaTrimming']['sum']['totalBases'])))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['qualityTrimming']['sum']['trimmedBases']),int(sample.stats['rubiconWgaTrimming']['sum']['totalBases'])))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        outputStr += '<tr>'
    outputStr += '</table>'

    outputStr += '<h3>Files (r1 then r2):</h3>'
    outputStr += '<table>'
    outputStr += '<tr>'
    outputStr += '<th>File Id</th>'
    outputStr += '<th>Sample Name</th>'
    outputStr += '<th colspan="2">Original Bases</th>'
    outputStr += '<th colspan="2">% rubicon adapter</th>'
    outputStr += '<th colspan="2">% malbac adapter</th>'
    outputStr += '<th colspan="2">% ampliOne adapter</th>'
    outputStr += '<th colspan="2">% illumina Adapters</th>'
    outputStr += '<th colspan="2">% quality trimmed</th>'
    outputStr += '</tr>'
    for filePairId,readCount,fastq1,fastq2,sampleId in sorted(analysispipe.database.getFastqs(), key=operator.itemgetter(0)):
        try:
          sample = samplesbyId[int(sampleId)]
          outputStr += '<tr>'
          outputStr += '<td>'+str(filePairId)+'</td>'
          outputStr += '<td>'+sample.name+'</td>'
          for read in ['r1','r2']:
              try: outputStr += '<td>'+thousandString(str(int(sample.stats['rubiconWgaTrimming'][filePairId][read]['totalBases'])))+'</td>'
              except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA</font></td>'
          for read in ['r1','r2']:
              try: outputStr += '<td>'+str(percentage(int(sample.stats['rubiconWgaTrimming'][filePairId][read]['trimmedBases']),int(sample.stats['rubiconWgaTrimming'][filePairId][read]['totalBases'])))+'%</td>'
              except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
          for read in ['r1','r2']:
              try: outputStr += '<td>'+str(percentage(int(sample.stats['malbacWgaTrimming'][filePairId][read]['trimmedBases']),int(sample.stats['rubiconWgaTrimming'][filePairId][read]['totalBases'])))+'%</td>'
              except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
          for read in ['r1','r2']:
              try: outputStr += '<td>'+str(percentage(int(sample.stats['ampliOneTrimming'][filePairId][read]['trimmedBases']),int(sample.stats['rubiconWgaTrimming'][filePairId][read]['totalBases'])))+'%</td>'
              except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
          for read in ['r1','r2']:
              try: outputStr += '<td>'+str(percentage(int(sample.stats['illuminaAndNexteraTrimming'][filePairId][read]['trimmedBases']),int(sample.stats['rubiconWgaTrimming'][filePairId][read]['totalBases'])))+'%</td>'
              except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
          for read in ['r1','r2']:
              try: outputStr += '<td>'+str(percentage(int(sample.stats['qualityTrimming'][filePairId][read]['trimmedBases']),int(sample.stats['rubiconWgaTrimming'][filePairId][read]['totalBases'])))+'%</td>'
              except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
          outputStr += '<tr>'
        except KeyError: pass
    outputStr += '</table>'

    outputStr += '<h2>Read Orientations:</h2>'
    outputStr += '<table>'
    outputStr += '<tr>'
    outputStr += '<th>Sample Id</th>'
    outputStr += '<th>Sample Name</th>'
    outputStr += '<th>Left After<br>All Filtering</th>'
    outputStr += '<th>Properly<br>Paired</th>'
    outputStr += '<th>Fwd Rev</th>'
    outputStr += '<th>Fwd Rev<br>full overlapp</th>'
    outputStr += '<th>Fwd Fwd</th>'
    outputStr += '<th>Rev Fwd</th>'
    outputStr += '<th>Rev Rev</th>'
    outputStr += '<th>r1 r2 on<br>diferent chrom</th>'
    outputStr += '</tr>'
    for sample in refSamplesFirst:
        outputStr += '<tr>'
        outputStr += '<td>'
        outputStr += str(sample.id)
        outputStr += '</td>'
        outputStr += '<td>'+sample.name+'</td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['noDuplicatesBamFlagstat']['mapped']),int(sample.stats['fixedBamFlagstat']['totalReads'])))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['orientations']['PP']), int(sum(sample.stats['orientations'].values()))))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['orientations']['FR']), int(sum(sample.stats['orientations'].values()))))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['orientations']['FRsp']), int(sum(sample.stats['orientations'].values()))))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['orientations']['FF']), int(sum(sample.stats['orientations'].values()))))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['orientations']['RF']), int(sum(sample.stats['orientations'].values()))))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['orientations']['RR']), int(sum(sample.stats['orientations'].values()))))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        try: outputStr += '<td>'+str(percentage(int(sample.stats['orientations']['difChrom']), int(sum(sample.stats['orientations'].values()))))+'%</td>'
        except (KeyError, TypeError) as e:outputStr += '<td><font color="red">NA%</font></td>'
        outputStr += '<tr>'
    outputStr += '</table>'

    outputStr += '<h2>Variants:</h2>'
    outputStr += '<table>'
    outputStr += '<tr>'
    outputStr += '<th>Sample Id</th>'
    outputStr += '<th>Sample Name</th>'
    outputStr += '<th>refSample</th>'
    outputStr += '<th>totalVars</th>'
    outputStr += '<th>% Hetero</th>'
    outputStr += '<th>% Dropout</th>'
    outputStr += '<th>% Other</th>'
    if os.path.exists(analysispipe.path+'/SNPidentificationTable.colored.txt'):
        outputStr += '<th>totalVars<br>Identification</th>'
        outputStr += '<th>% Patient</th>'
        outputStr += '<th>% Donor</th>'
        outputStr += '<th>% Mix</th>'
        outputStr += '<th>% Other</th>'
        outputStr += '<th>Classification</th>'
        outputStr += '<th>Low Confidence</th>'
    outputStr += '</tr>'
    for sample in refSamplesFirst:
        outputStr += '<tr>'
        outputStr += '<td>'+str(sample.id)+'</td>'
        outputStr += '<td>'+sample.name+'</td>'
        if 'refsample' in sample.stats and sample.stats['refsample'] != 'NA':outputStr += '<td>'+str(sample.stats['refsample'])+'</td>'
        else: outputStr += '<td><font color="red">NA%</font></td>'
        if 'totalAdoVar' in sample.stats and sample.stats['totalAdoVar'] != 'NA':outputStr += '<td>'+str(sample.stats['totalAdoVar'])+'</td>'
        else: outputStr += '<td><font color="red">NA%</font></td>'
        if 'correct%' in sample.stats and sample.stats['correct%'] != 'NA':outputStr += '<td>'+str(sample.stats['correct%'])+'%</td>'
        else: outputStr += '<td><font color="red">NA%</font></td>'
        if 'dropout%' in sample.stats and sample.stats['dropout%'] != 'NA':outputStr += '<td>'+str(sample.stats['dropout%'])+'%</td>'
        else: outputStr += '<td><font color="red">NA%</font></td>'
        if 'other%' in sample.stats and sample.stats['other%'] != 'NA':outputStr += '<td>'+str(sample.stats['other%'])+'%</td>'
        else: outputStr += '<td><font color="red">NA%</font></td>'
        if os.path.exists(analysispipe.path+'/SNPidentificationTable.colored.txt'):
          if sample.stats['totalIdVar'] != 'NA':outputStr += '<td>'+str(sample.stats['totalIdVar'])+'</td>'
          else: outputStr += '<td><font color="red">NA%</font></td>'
          if sample.stats['pat%'] != 'NA':outputStr += '<td>'+str(sample.stats['pat%'])+'%</td>'
          else: outputStr += '<td><font color="red">NA%</font></td>'
          if sample.stats['don%'] != 'NA':outputStr += '<td>'+str(sample.stats['don%'])+'%</td>'
          else: outputStr += '<td><font color="red">NA%</font></td>'
          if sample.stats['mix%'] != 'NA':outputStr += '<td>'+str(sample.stats['mix%'])+'%</td>'
          else: outputStr += '<td><font color="red">NA%</font></td>'
          if sample.stats['other%'] != 'NA':outputStr += '<td>'+str(sample.stats['other%'])+'%</td>'
          else: outputStr += '<td><font color="red">NA%</font></td>'
          if sample.stats['sampleIs'] != 'NA':outputStr += '<td>'+str(sample.stats['sampleIs'])+'</td>'
          else: outputStr += '<td><font color="red">NA%</font></td>'
          if sample.stats['sampleIsLowConf'] != 'NA':outputStr += '<td>'+str(sample.stats['sampleIsLowConf'])+'</td>'
          else: outputStr += '<td><font color="red">NA%</font></td>'

        outputStr += '<tr>'
    outputStr += '</table>'

    outputStr += '<h2>InsertSize graphs:</h2>'
    if os.path.exists(analysispipe.path+'/graphics/insertSizes.png'):
        outputStr += '<a href="'+'static/runData'+'/all_samples_graphics/insertSizes.pdf'+'">'
        outputStr += '<img src="'+'static/runData'+'/all_samples_graphics/insertSizes.png'+'">'
        outputStr += '</a>'
    else: outputStr += '<font color="red">Insert size graph not created yet.</font><br>'
    outputStr += '<h2>Coverage graphs:</h2>'
    if os.path.exists(analysispipe.path+'/graphics/lorentzCurve.png'):
        outputStr += '<a href="'+'static/runData'+'/all_samples_graphics/lorentzCurve.pdf'+'">'
        outputStr += '<img src="'+'static/runData'+'/all_samples_graphics/lorentzCurve.png'+'">'
        outputStr += '</a>'
    else: outputStr += '<font color="red">Lorentz curve graph not created yet.</font><br>'
    if os.path.exists(analysispipe.path+'/graphics/lorentzCurveCoveredOnly.png'):
        outputStr += '<a href="'+'static/runData'+'/all_samples_graphics/lorentzCurveCoveredOnly.pdf'+'">'
        outputStr += '<img src="'+'static/runData'+'/all_samples_graphics/lorentzCurveCoveredOnly.png'+'">'
        outputStr += '</a>'
    else: outputStr += '<font color="red">Lorentz curve graph for CoveredOnly not created yet.</font><br>'
    if os.path.exists(analysispipe.path+'/graphics/exomecoverage.png'):
        outputStr += '<a href="'+'static/runData'+'/all_samples_graphics/exomecoverage.pdf'+'">'
        outputStr += '<img src="'+'static/runData'+'/all_samples_graphics/exomecoverage.png'+'">'
        outputStr += '</a>'
    else: outputStr += '<font color="red">Exome coverage graph not created yet.</font><br>'
    outputStr += '<br>'
    
    if not list(glob.iglob(analysispipe.path+'/samples/*/plots/readDepth.*.png')):
        outputStr += '<font color="red">No read depth graphs created yet.</font><br>'
    else:
        for chrom in sorted_nicely(lengthOfChromByName.keys()):
          if chrom[0] == 'G':continue
          outputStr += '<br><br>Coverage over chromosome '+chrom+':<br>'
          for sample in refSamplesFirst:
              if os.path.exists(sample.plotsPath+'/readDepth.'+sample.name+'.'+chrom+'.png'):
                outputStr += '<a href="'+'static/runData'+'/per_sample_info/'+samplesbyId[int(sampleId)].name+'/plots'+'/readDepth.'+sample.name+'.'+chrom+'.pdf'+'">'
                #os.path.relpath(sample.plotsPath+'/readDepth.'+sample.name+'.'+chrom+'.pdf','/'.join(reportFile.name.split('/')[:-1]))+'">'
                outputStr += '<img src="'+'static/runData'+'/per_sample_info/'+samplesbyId[int(sampleId)].name+'/plots'+'/readDepth.'+sample.name+'.'+chrom+'.png'+'">'
                outputStr += '</a><br>'

    #for imgFileName in glob.iglob(analysispipe.path+'/graphics/readDepth.*.png'):
    #    relPath = os.path.relpath(imgFileName, '/'.join(reportFile.name.split('/')[:-1])
    #    outputStr += '<br>Coverage over chromosome '+imgFileName.split('/readDepth.')[-1].split('.png')[0]+':<br>'
    #    outputStr += '<img src="'+relPath+'">'

    outputStr += '</body></html>\n'
    analysispipe.logfile.write('#LOGMSG#'+time.strftime("%Y-%m-%d:%H:%M:%S",time.localtime())+'#'+str(analysispipe.masterPid)+'# Full report now written to '+analysispipe.path+'/report.htm'+' \n')
    return outputStr

def reportParallel(sample):
    sample.getStats()
    return sample

def createSampleSummaryCsv(analysispipe):
    csv = "name"+ '\n'#"name"+","+"thing" + '\n'
    for sample in analysispipe.database.getSamples():
        csv += sample.name + '\n'
    # csv += 'kurt'+','+'gurka' + '\n'
    # csv += 'goran'+','+'smurgas' + '\n'
    # csv += 'svenne'+','+'fisk' + '\n'
    return csv
