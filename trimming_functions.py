class SampleTrimmer():
    
    def __init__(self, analysispipe):
        self.analysispipe = analysispipe # link/connection to the analysispipe object
        
    def trimming(self, sample):
        for filePairId,readCount,fastq1,fastq2,sampleId in sample.getFastqs():
            output_file = open(sample.scriptPath+'/trimming.'+sample.name+'.'+str(filePairId)+'.sbatch.sh', 'w')
            output = '#! /bin/bash -l'+'\n'
        output += '#SBATCH -A '+self.analysispipe.settings.uppmaxProject+'\n'
        output += '#SBATCH -n 2 -p node'+'\n'
        output += '#SBATCH -t 240:00:00'+'\n'
        output += '#SBATCH -J trim.'+sample.name+'.'+str('test')+'\n'
        output += '#SBATCH -e '+sample.logPath+'/stderr.trimming.'+sample.name+'.'+str(filePairId)+'.txt'+'\n'
        output += '#SBATCH -o '+sample.logPath+'/stdout.trimming.'+sample.name+'.'+str(filePairId)+'.txt'+'\n'

        #
        # Go to sample path, load tools and define variables
        #
        output += 'echo "$(date) Running on: $(hostname)"'+'\n'
        output += 'cd '+sample.path+'\n'
        output += 'echo'+'\n'
        #output += 'module load bioinfo-tools pysam/0.8.3-py27 FastQC cutadapt/1.8.0'+'\n'
        output += 'module load bioinfo-tools pysam/0.8.3-py27 FastQC cutadapt/1.8.0 bowtie2 samtools picard/1.114 BEDTools/2.16.2  GATK/3.1.1'+'\n'
        output += 'module load python/2.7'+'\n'
        r1_in = fastq1
        r2_in = fastq2
        #
        # Rubicon trimming
        #
        if not self.analysispipe.settings.skiprubicon:
            output += ''+self.analysispipe.settings.trim_scripts+'/wgaAdapterTrimmer.py -i '+r1_in+' > '+sample.tempPath+'/'+str(filePairId)+'.r1.wgaTrimmed.fq 2> '+sample.logPath+'/rubiconWgaTrimming.'+str(filePairId)+'.r1.log.txt &\n'
            output += ''+self.analysispipe.settings.trim_scripts+'/wgaAdapterTrimmer.py -i '+r2_in+' > '+sample.tempPath+'/'+str(filePairId)+'.r2.wgaTrimmed.fq 2> '+sample.logPath+'/rubiconWgaTrimming.'+str(filePairId)+'.r2.log.txt &\n'
            output += 'wait\n'
            
            r1_in = sample.tempPath+'/'+str(filePairId)+'.r1.wgaTrimmed.fq'
            r2_in = sample.tempPath+'/'+str(filePairId)+'.r2.wgaTrimmed.fq'
        #
        # Malbac trimming 
        #
        if not self.analysispipe.settings.skipmalbac:
            output += '\n'
            output += 'cutadapt -n 10 -g GTGAGTGATGGTTGAGGTAGTGTGGAG -a CTCCACACTACCTCAACCATCACTCAC '+r1_in+' > '+sample.tempPath+'/'+str(filePairId)+'.r1.wgaTrimmed2.fq  2> '+sample.logPath+'/malbacWgaTrimming.'+str(filePairId)+'.r1.log.txt &\n'
            output += 'cutadapt -n 10 -g GTGAGTGATGGTTGAGGTAGTGTGGAG -a CTCCACACTACCTCAACCATCACTCAC '+r2_in+' > '+sample.tempPath+'/'+str(filePairId)+'.r2.wgaTrimmed2.fq  2> '+sample.logPath+'/malbacWgaTrimming.'+str(filePairId)+'.r2.log.txt &\n'
            output += 'wait\n'
            output += '\n'
            r1_in = sample.tempPath+'/'+str(filePairId)+'.r1.wgaTrimmed2.fq'
            r2_in = sample.tempPath+'/'+str(filePairId)+'.r2.wgaTrimmed2.fq'
        
            if not self.analysispipe.settings.skiprubicon:output += 'rm -v '+sample.tempPath+'/'+str(filePairId)+'.r1.wgaTrimmed.fq '+sample.tempPath+'/'+str(filePairId)+'.r2.wgaTrimmed.fq\n'
        #
        # Ampli1 trimming
        #
        if not self.analysispipe.settings.skipampli1:
            output += '\n'
            output += 'cutadapt -n 10 -g AGTGGGATTCCTGCTGTCAGT '+r1_in+' > '+sample.tempPath+'/'+str(filePairId)+'.r1.wgaTrimmed3.fq  2> '+sample.logPath+'/ampli1WgaTrimming.'+str(filePairId)+'.r1.log.txt &\n'
            output += 'cutadapt -n 10 -g AGTGGGATTCCTGCTGTCAGT '+r2_in+' > '+sample.tempPath+'/'+str(filePairId)+'.r2.wgaTrimmed3.fq  2> '+sample.logPath+'/ampli1WgaTrimming.'+str(filePairId)+'.r2.log.txt &\n'
            r1_in = sample.tempPath+'/'+str(filePairId)+'.r1.wgaTrimmed3.fq'
            r2_in = sample.tempPath+'/'+str(filePairId)+'.r2.wgaTrimmed3.fq'
            output += 'wait\n'
            output += '\n'
            if not self.analysispipe.settings.skiprubicon and self.analysispipe.settings.skipmalbac:
                output += 'rm -v '+sample.tempPath+'/'+str(filePairId)+'.r1.wgaTrimmed.fq '+sample.tempPath+'/'+str(filePairId)+'.r2.wgaTrimmed.fq\n'
            if not self.analysispipe.settings.skipmalbac:
                output += 'rm -v '+sample.tempPath+'/'+str(filePairId)+'.r1.wgaTrimmed2.fq '+sample.tempPath+'/'+str(filePairId)+'.r2.wgaTrimmed2.fq\n'
        #
        # Illumina  adapter trimming
        #
        adaptersToTrim = '-a CTGTCTCTTATACACATCTGACGCTGCCGACGA -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
        output += 'cutadapt -n 3 '+adaptersToTrim+' '+r1_in+' > '+sample.tempPath+'/'+str(filePairId)+'.r1.wgaAndilluminaTrimmed.fq 2> '+sample.logPath+'/illuminaAndNexteraTrimming.'+str(filePairId)+'.r1.log.txt &\n'
        output += 'cutadapt -n 3 '+adaptersToTrim+' '+r2_in+' > '+sample.tempPath+'/'+str(filePairId)+'.r2.wgaAndilluminaTrimmed.fq 2> '+sample.logPath+'/illuminaAndNexteraTrimming.'+str(filePairId)+'.r2.log.txt &\n'
        output += 'wait\n'
        #
        # Remove temporary files
        #
        output += 'rm -v '+sample.tempPath+'/'+str(filePairId)+'.r1.wgaTrimmed3.fq '+sample.tempPath+'/'+str(filePairId)+'.r2.wgaTrimmed3.fq \n'
        output += 'wait\n'
        #
        # Quality trimmming
        #
        output += ''+self.analysispipe.settings.trim_scripts+'/TrimBWAstyle.pl -q 20 '+sample.tempPath+'/'+str(filePairId)+'.r1.wgaAndilluminaTrimmed.fq > '+sample.tempPath+'/'+str(filePairId)+'.r1.wgaIlluminaAndQualityTrimmed.fq 2> '+sample.logPath+'/qualityTrimming.'+str(filePairId)+'.r1.log.txt &\n'
        output += ''+self.analysispipe.settings.trim_scripts+'/TrimBWAstyle.pl -q 20 '+sample.tempPath+'/'+str(filePairId)+'.r2.wgaAndilluminaTrimmed.fq > '+sample.tempPath+'/'+str(filePairId)+'.r2.wgaIlluminaAndQualityTrimmed.fq 2> '+sample.logPath+'/qualityTrimming.'+str(filePairId)+'.r2.log.txt &\n'
        output += 'wait\n'
        #
        # Remove temporary files
        #
        output += 'rm -v '+sample.tempPath+'/'+str(filePairId)+'.r1.wgaAndilluminaTrimmed.fq '+sample.tempPath+'/'+str(filePairId)+'.r2.wgaAndilluminaTrimmed.fq\n'
        output += 'wait\n'
        #
        # Remove empty or "N" only sequences
        #
        output += 'python '+self.analysispipe.settings.trim_scripts+'/removeEmptyReads.py '
        output += sample.tempPath+'/'+str(filePairId)+'.r1.wgaIlluminaAndQualityTrimmed.fq '
        output += sample.tempPath+'/'+str(filePairId)+'.r2.wgaIlluminaAndQualityTrimmed.fq '
        output += sample.dataPath+'/'+str(filePairId)+'.r1.allTrimmed.fq '
        output += sample.dataPath+'/'+str(filePairId)+'.r2.allTrimmed.fq '
        output += sample.dataPath+'/'+str(filePairId)+'.singletts.fq '
        output += '>&2 2> '+sample.logPath+'/removeEmptyReads.'+str(filePairId)+'.log.txt\n'
        #
        # Compress files
        #
        output += 'gzip -v9 '+sample.dataPath+'/'+str(filePairId)+'.r1.allTrimmed.fq &\n'
        output += 'gzip -v9 '+sample.dataPath+'/'+str(filePairId)+'.r2.allTrimmed.fq  &\n'
        output += 'gzip -v9 '+sample.dataPath+'/'+str(filePairId)+'.singletts.fq &\n'
        output += 'wait\n'
        #
        # FASTQC
        #
        output += 'fastqc '+sample.dataPath+'/'+str(filePairId)+'.r1.allTrimmed.fq.gz &\n'
        output += 'fastqc '+sample.dataPath+'/'+str(filePairId)+'.r2.allTrimmed.fq.gz &\n'
        output += 'fastqc '+sample.dataPath+'/'+str(filePairId)+'.singletts.fq.gz &\n'
        output += 'wait\n'
        output += 'mv -v '+sample.dataPath+'/*fastqc* '+sample.fastqcPath+'/\n'
        #
        # Final output and write script to file
        #
        output += 'echo'+'\n'
        output += 'wait'+'\n'
        output += 'echo "$(date) AllDone"'+'\n'
        output += 'echo "$(date) AllDone" >&2'+'\n'
        
        output_file.write(output)
        output_file.close()
            


    
