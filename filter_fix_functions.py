class SampleFilterAndFix():
    
    def __init__(self, analysispipe):
        self.analysispipe = analysispipe # link/connection to the analysispipe object
        
        
    def filter_and_fix(self, sample):
        #for filePairId,readCount,fastq1,fastq2,sampleId in sample.getFastqs():
        output_file = open(sample.scriptPath+'/fnf.'+sample.name+'.sbatch.sh', 'w')
            
        #
        # sbatch header
        #
        output = '#! /bin/bash -l'+'\n'
        output += '#SBATCH -A '+self.analysispipe.settings.uppmaxProject+'\n'
        output += '#SBATCH -n 1 -p core'+'\n'
        output += '#SBATCH -t 48:00:00'+'\n'
        output += '#SBATCH -J fnf.'+sample.name+'\n'
        output += '#SBATCH -e '+sample.logPath+'/stderr.filterAndFix.'+sample.name+'.txt'+'\n'
        output += '#SBATCH -o '+sample.logPath+'/stdout.filterAndFix.'+sample.name+'.txt'+'\n'

        #
        # define variebles and go to path
        #
        output += 'echo "$(date) Running on: $(hostname)"'+'\n'
        output += 'cd '+sample.path+'\n'
        output += 'echo'+'\n'
        output += 'module load bioinfo-tools pysam/0.8.3-py27  FastQC cutadapt/1.8.0 bowtie2 samtools picard/1.114 BEDTools/2.16.2  GATK/3.1.1'+'\n'
        output += 'module load python/2.7'+'\n'
        
        
        #
        # sort the bam file
        #
        output += 'java -Xmx5g -jar '+self.analysispipe.settings.picardLocation+'/SortSam.jar MAX_RECORDS_IN_RAM=2500000 SORT_ORDER=coordinate '
        output += 'INPUT='+ sample.dataPath+'/'+sample.name+'.merged.bam '
        output += 'OUTPUT='+sample.tempPath+'/'+sample.name+'.sorted.bam '
        output += 'CREATE_INDEX=true 1>&2  2> '
        output += sample.logPath+'/stderr.sortBam.'+sample.name+'.txt \n'
        output += 'echo -e "bam2sort Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
        output += 'rm -v '+sample.dataPath+'/'+sample.name+'.merged.bam\n'

        #
        # mark duplicates
        #
        output += 'java -Xmx5g -jar '+self.analysispipe.settings.picardLocation+'/MarkDuplicates.jar MAX_RECORDS_IN_RAM=2500000 VALIDATION_STRINGENCY=LENIENT '
        output += 'INPUT='+ sample.tempPath+'/'+sample.name+'.sorted.bam '
        output += 'OUTPUT='+sample.tempPath+'/'+sample.name+'.marked.bam '
        output += 'METRICS_FILE='+sample.logPath+'/markDuplicatesMetrix.'+sample.name+'.txt '
        output += '1>&2  2> '+sample.logPath+'/stderr.markDuplicates.'+sample.name+'.txt \n'
        output += 'echo -e "mark Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
        output += 'rm -v '+sample.tempPath+'/'+sample.name+'.sorted.bam\n'

        #
        # fix missing information
        #
        output += 'java -Xmx5g -jar '+self.analysispipe.settings.picardLocation+'/AddOrReplaceReadGroups.jar '
        output += 'MAX_RECORDS_IN_RAM=2500000 '
        output += 'INPUT='+ sample.tempPath+'/'+sample.name+'.marked.bam '
        output += 'OUTPUT='+sample.dataPath+'/'+sample.name+'.fixed.bam '
        output += 'CREATE_INDEX=true RGID='+sample.name+' RGLB='+sample.name+' RGPL=ILLUMINA RGSM='+sample.name+' RGCN="NA" RGPU="NA"'+'  '
        output += '1>&2  2> '+sample.logPath+'/stderr.addAndReplaceReadGroups.'+sample.name+'.txt \n'
        output += 'echo "addorreplace Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
        output += 'rm -v '+sample.tempPath+'/'+sample.name+'.marked.bam\n'

        #
        # samtools flagstat
        #
        output += 'samtools flagstat '+sample.dataPath+'/'+sample.name+'.fixed.bam'+' > '+sample.logPath+'/fixedBamFlagstat.'+sample.name+'.txt \n'
        output += 'echo "flagstat Done. $(date) Running on: $(hostname)" 1>&2'+'\n'

        #
        # Final output and write script to file
        #
        #output += '\n'+AnalysisPipe.programPath+' '+AnalysisPipe.path+' report\n'
        output += 'wait'+'\n'
        output += 'echo "$(date) AllDone"'+'\n'
        output += 'echo "$(date) AllDone" >&2'+'\n'
        
        output_file.write(output)
        output_file.close()
        
        

