class SampleMerger():
    
    def __init__(self, analysispipe):
        self.analysispipe = analysispipe # link/connection to the analysispipe object
        
    def merge_mapped_reads(self, sample):
        #for filePairId,readCount,fastq1,fastq2,sampleId in sample.getFastqs():
        output_file = open(sample.scriptPath+'/merge.'+sample.name+'.sbatch.sh', 'w')
        #
        # Sbatch header
        #
        output = '#! /bin/bash -l'+'\n'
        output += '#SBATCH -A '+self.analysispipe.settings.uppmaxProject+'\n'
        output += '#SBATCH -n 1 -p core'+'\n'
        output += '#SBATCH -t 240:00:00'+'\n'
        output += '#SBATCH -J merge.'+sample.name+'.'+str('test')+'\n'
        output += '#SBATCH -e '+sample.logPath+'/stderr.merge.'+sample.name+'.txt'+'\n'
        output += '#SBATCH -o '+sample.logPath+'/stdout.merge.'+sample.name+'.txt'+'\n'
        #
        # define variebles and go to path
        #
        output += 'echo "$(date) Running on: $(hostname)"'+'\n'
        output += 'cd '+sample.path+'\n'
        output += 'echo'+'\n'
        output += 'module load bioinfo-tools pysam/0.8.3-py27 FastQC cutadapt/1.8.0 bowtie2 samtools picard/1.114 BEDTools/2.16.2  GATK/3.1.1'+'\n'
        output += 'module load python/2.7'+'\n'
        #
        # merge
        #

        if len(list(sample.getFastqs())) > 1:
            inputFiles = ' INPUT='+' INPUT='.join([sample.dataPath+'/'+str(filePairId)+'.bam' for filePairId,readCount,fastq1,fastq2,sampleId in sample.getFastqs()])
            output += 'java -Xmx5g -jar '+self.analysispipe.settings.picardLocation+'/MergeSamFiles.jar '+inputFiles+' OUTPUT='+sample.dataPath+'/'+sample.name+'.merged.bam '
            output += '1>&2  2>  '+sample.logPath+'/stderr.merging.'+sample.name+'.txt \n'
            output += 'echo -e "mapping Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
            output += 'rm -v '+' '.join([sample.dataPath+'/'+str(filePairId)+'.bam' for filePairId,readCount,fastq1,fastq2,sampleId in sample.getFastqs()])+'\n'
        elif len(list(sample.getFastqs())) == 1:
            output += 'mv -v '+' '.join([sample.dataPath+'/'+str(filePairId)+'.bam' for filePairId,readCount,fastq1,fastq2,sampleId in [list(sample.getFastqs())[0]]])+' '+sample.dataPath+'/'+sample.name+'.merged.bam\n'
        else: print 'ERROR: no fastqfiles for sample '+sample.name

        #
        # Final output and write script to file
        #
        output += 'wait'+'\n'
        output += 'echo "$(date) AllDone"'+'\n'
        output += 'echo "$(date) AllDone" >&2'+'\n'
        output_file.write(output)
        output_file.close()
    



