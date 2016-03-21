class SampleMapper():

    def __init__(self, analysispipe):
       self.analysispipe = analysispipe # link/connection to the analysispipe object

    def Bowtie2_mapping(self, sample):
        for filePairId,readCount,fastq1,fastq2,sampleId in sample.getFastqs():
            output_file = open(sample.scriptPath+'/mapping.'+str(filePairId)+'.sbatch.sh', 'w')
            output = '#! /bin/bash -l'+'\n'
	    output += '#SBATCH -A '+self.analysispipe.settings.uppmaxProject+'\n'
	    output += '#SBATCH -n 16 -p node'+'\n'
	    output += '#SBATCH -t 240:00:00'+'\n'
	    output += '#SBATCH -J map.'+sample.name+'.'+str('test')+'\n'
	    output += '#SBATCH -e '+sample.logPath+'/stderr.mapping.'+sample.name+'.'+str(filePairId)+'.txt'+'\n'
	    output += '#SBATCH -o '+sample.logPath+'/stdout.mapping.'+sample.name+'.'+str(filePairId)+'.txt'+'\n'
            #
	    # define variebles and go to path
	    #
	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'cd '+sample.path+'\n'
	    output += 'echo'+'\n'
	    
	    output += 'module load bioinfo-tools pysam/0.8.3-py27 FastQC cutadapt/1.8.0 bowtie2 samtools picard/1.114 BEDTools/2.16.2  GATK/3.1.1'+'\n'
	    output += 'module load python/2.7'+'\n'

	    #
	    # Bowtie2 mapping
	    #output += 'module load bioinfo-tools pysam bwa/0.7.8\n'
	    #output += 'bwa mem -t 16 /sw/data/uppnex/reference/Homo_sapiens/GRCh37/program_files/bwa/concat.fa '+sample.r1files[0]+' '+sample.r2files[0]+' > '+sample.sam+'\n'
	    #output += 'bowtie2 -1 '+sample.r1files[0]+' -2 '+sample.r2files[0]+' --very-sensitive-local -p16 -x '+sample.reference+' > '+sample.sam+'\n'
	    output += 'bowtie2 --maxins 2000 -p16 '
	    output += '-1 '+sample.dataPath+'/'+str(filePairId)+'.r1.allTrimmed.fq.gz '
	    output += '-2 '+sample.dataPath+'/'+str(filePairId)+'.r2.allTrimmed.fq.gz '
#	    output += '-U '+sample.dataPath+'/'+str(filePairId)+'.singletts.fq.gz '
	    output += '-x '+self.analysispipe.settings.bowtie2Reference+' '
	    #output += '> '+sample.tempPath+'/'+str(filePairId)+'.sam '
	    output += '2> '+sample.logPath+'/stderr.bowtie2.'+str(filePairId)+'.txt |'
	    #
	    # convert to bam file
	    #
	    output += 'java -Xmx5g -jar '+self.analysispipe.settings.picardLocation+'/SamFormatConverter.jar MAX_RECORDS_IN_RAM=2500000 '
	    output += 'INPUT='+ '/dev/stdin '#sample.tempPath+'/'+str(filePairId)+'.sam '
	    output += 'OUTPUT='+sample.dataPath+'/'+str(filePairId)+'.bam '
	    output += '1>&2  2> '+sample.logPath+'/stderr.'+str(filePairId)+'.sam2bam.'+sample.name+'.txt \n'
	    output += 'echo -e "sam2bam Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
	    output += 'rm -v '+ sample.dataPath+'/'+str(filePairId)+'.sam\n'

	    output += 'echo -e "mapping Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
	    output += 'rm -v '+sample.dataPath+'/'+str(filePairId)+'.r1.allTrimmed.fq.gz'+'\n'
	    output += 'rm -v '+sample.dataPath+'/'+str(filePairId)+'.r2.allTrimmed.fq.gz'+'\n'
	    output += 'rm -v '+sample.dataPath+'/'+str(filePairId)+'.singletts.fq.gz'+'\n'

	    #
	    # Final output and write script to file
	    #
	    #output += '\n'+self.analysispipe.programPath+' '+self.analysispipe.path+' report\n'
	    output += 'wait'+'\n'
	    output += 'echo "$(date) AllDone"'+'\n'
	    output += 'echo "$(date) AllDone" >&2'+'\n'

            output_file.write(output)
            output_file.close()
       







