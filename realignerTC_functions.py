class SampleRealignTargetCreator():
    
    def __init__(self, analysispipe):
        self.analysispipe = analysispipe # link/connection to the analysispipe object
        
        
    def realign_target_creator(self, sample):
        #for filePairId,readCount,fastq1,fastq2,sampleId in sample.getFastqs():
        #output_file = open(sample.scriptPath+'/realTC.'+sample.name+'.'+str(filePairId)+'.sbatch.sh', 'w')
        output_file = open(sample.scriptPath+'/realTC.'+sample.name+'.sbatch.sh', 'w')
            
        #
        # sbatch header
        #
        output = '#! /bin/bash -l'+'\n'
        output += '#SBATCH -A '+self.analysispipe.settings.uppmaxProject+'\n'
        output += '#SBATCH -n 1 -p core'+'\n'
        output += '#SBATCH -t 240:00:00'+'\n'
        output += '#SBATCH -J realTC.'+sample.name+'\n'
        output += '#SBATCH -e '+sample.logPath+'/stderr.realTC.'+sample.name+'.txt'+'\n'
        output += '#SBATCH -o '+sample.logPath+'/stdout.realTC.'+sample.name+'.txt'+'\n'

        #
        # define variebles and go to path
        #
        output += 'echo "$(date) Running on: $(hostname)"'+'\n'
        output += 'cd '+sample.path+'\n'
        output += 'echo'+'\n'
        output += 'module load bioinfo-tools pysam/0.8.3-py27  FastQC cutadapt/1.8.0 bowtie2 samtools picard/1.114 BEDTools/2.16.2  GATK/3.1.1'+'\n'
        output += 'module load python/2.7'+'\n'
        
        #
        # Find targets for indel realignment
        #
        output += 'echo -e "-> RealignerTargetCreator <-"\n'
        output += 'java -Xmx72g -jar '+self.analysispipe.settings.GATKlocation+' -T RealignerTargetCreator '
        output += '-nt 16 '
        output += '-I '+sample.dataPath+'/'+sample.name+'.fixed.bam'+' '
        output += '-R '+self.analysispipe.settings.bowtie2Reference+' '
        output += '-o '+sample.dataPath+'/'+sample.name+'.reAlignemntTargetIntervals.bed '
        output += ' -known '+self.analysispipe.settings.GATKbundleLocation+'/Mills_and_1000G_gold_standard.indels.b37.vcf'
        output += ' -known '+self.analysispipe.settings.GATKbundleLocation+'/1000G_phase1.indels.b37.vcf '
        output += '1>&2 2> '+sample.logPath+'/stderr.RealignerTargetCreator.'+sample.name+'.txt;'
        output += '\n'
        
        output_file.write(output)
        output_file.close()