class SampleReAlignAndReCalibrator():
    
    def __init__(self, analysispipe):
        self.analysispipe = analysispipe # link/connection to the analysispipe object
        
        
    def realign_and_recalibrate(self, sample):
        for filePairId,readCount,fastq1,fastq2,sampleId in sample.getFastqs():
            output_file = open(sample.scriptPath+'/reAlign.'+sample.name+'.sbatch.sh', 'w')
            
        #
        # sbatch header
        #
        output = '#! /bin/bash -l'+'\n'
        output += '#SBATCH -A '+self.analysispipe.settings.uppmaxProject+'\n'
        output += '#SBATCH -n 1 -p core'+'\n'
        output += '#SBATCH -t 240:00:00'+'\n'
        output += '#SBATCH -J reAlign.'+sample.name+'\n'
        output += '#SBATCH -e '+sample.logPath+'/stderr.reAlign.'+sample.name+'.txt'+'\n'
        output += '#SBATCH -o '+sample.logPath+'/stdout.reAlign.'+sample.name+'.txt'+'\n'

        #
        # define variebles and go to path
        #
        output += 'echo "$(date) Running on: $(hostname)"'+'\n'
        output += 'cd '+sample.path+'\n'
        output += 'echo'+'\n'
        output += 'module load bioinfo-tools pysam/0.8.3-py27  FastQC cutadapt/1.8.0 bowtie2 samtools picard/1.114 BEDTools/2.16.2  GATK/3.1.1'+'\n'
        output += 'module load python/2.7'+'\n'
        
        
        
        #
        # Realign reads around indels
        #
        output += 'echo -e "-> Running IndelRealigner <-"\n'
        output += 'java -Xmx5g -jar '+self.analysispipe.settings.GATKlocation+' -T IndelRealigner '
        output += '-I '+sample.dataPath+'/'+sample.name+'.fixed.bam'+' '
        output += '-R '+self.analysispipe.settings.GATK_reference+' '
        output += '-targetIntervals '+sample.dataPath+'/'+sample.name+'.reAlignemntTargetIntervals.bed '
        output += ' -o '+sample.tempPath+'/'+sample.name+'.reAligned.bam'+' '
        output += ' -known '+self.analysispipe.settings.GATKbundleLocation+'/Mills_and_1000G_gold_standard.indels.b37.vcf'
        output += ' -known '+self.analysispipe.settings.GATKbundleLocation+'/1000G_phase1.indels.b37.vcf  '
        output += '1>&2 2> '+sample.logPath+'/stderr.indelRealigner.'+sample.name+'.txt;'+'\n'
        output += '\n'
        output += 'echo "Done. $(date) Running on: $(hostname)"\n'
        output += 'echo "$(date) Running on: $(hostname)"\n'
        output += 'rm -v '+sample.dataPath+'/'+sample.name+'.fixed.bam'+'\n'
        
        #
        # Quality recalibration
        #
        output += 'echo -e "-> Running BaseRecalibrator <-"\n'
        output += 'java -Xmx5g -jar '+self.analysispipe.settings.GATKlocation+' -T BaseRecalibrator '
        output += '-I '+sample.tempPath+'/'+sample.name+'.reAligned.bam'+' '
        output += '-R '+self.analysispipe.settings.GATK_reference+' '
        output += '-o '+sample.dataPath+'/'+sample.name+'.BQSR.grp'+' '
        output += ' -knownSites '+self.analysispipe.settings.GATKbundleLocation+'/dbsnp_138.b37.vcf '
        output += '1>&2 2> '+sample.logPath+'/stderr.baseRecalibrator.'+sample.name+'.txt;'+'\n'

        output += '\n'
        output += 'echo -e "-> Running PrintReads <-"\n'
        output += 'java -Xmx5g -jar '+self.analysispipe.settings.GATKlocation+' -T PrintReads '
        output += '-I '+sample.tempPath+'/'+sample.name+'.reAligned.bam'+' '
        output += '-R '+self.analysispipe.settings.GATK_reference+' '
        output += '-BQSR '+sample.dataPath+'/'+sample.name+'.BQSR.grp'+' '
        output += '-o '+sample.tempPath+'/'+sample.name+'.reCalibrated.bam'+' '
        output += '1>&2 2> '+sample.logPath+'/stderr.printreads.txt ;\n'
        output += 'rm -v '+sample.tempPath+'/'+sample.name+'.reAligned.bam'+'\n'
        output += 'samtools flagstat '+sample.tempPath+'/'+sample.name+'.reCalibrated.bam > '+sample.logPath+'/reCalibratedBamFlagstat.'+sample.name+'.txt \n'

        output += 'samtools view -b -F 4 '   +sample.tempPath+'/'+sample.name+'.reCalibrated.bam > '+sample.tempPath+'/'+sample.name+'.unmapRemoved.bam  2> '+sample.logPath+'/stderr.samtoolsView.removeUnmap.'+sample.name+'.txt \n'
        output += 'java -Xmx5g -jar '+self.analysispipe.settings.picardLocation+'/BuildBamIndex.jar INPUT='+sample.tempPath+'/'+sample.name+'.unmapRemoved.bam '+'1>&2  2>  '+sample.logPath+'/stderr.buildIndex1.'+sample.name+'.txt \n'
        output += 'samtools flagstat '+sample.tempPath+'/'+sample.name+'.unmapRemoved.bam > '+sample.logPath+'/unmapRemovedBamFlagstat.'+sample.name+'.txt \n'
        output += 'rm -v '+sample.tempPath+'/'+sample.name+'.reCalibrated.bam\n'

        output += 'samtools view -b -q 20 '  +sample.tempPath+'/'+sample.name+'.unmapRemoved.bam > '+sample.tempPath+'/'+sample.name+'.qualFiltered.bam  2> '+sample.logPath+'/stderr.samtoolsView.qualFilter.'+sample.name+'.txt \n'
        output += 'java -Xmx5g -jar '+self.analysispipe.settings.picardLocation+'/BuildBamIndex.jar INPUT='+sample.tempPath+'/'+sample.name+'.qualFiltered.bam '+'1>&2  2>  '+sample.logPath+'/stderr.buildIndex2.'+sample.name+'.txt \n'
        output += 'samtools flagstat '+sample.tempPath+'/'+sample.name+'.qualFiltered.bam > '+sample.logPath+'/qualFilteredBamFlagstat.'+sample.name+'.txt \n'
        output += 'rm -v '+sample.tempPath+'/'+sample.name+'.unmapRemoved.bam\n'

        output += 'samtools view -b -F 1024 '+sample.tempPath+'/'+sample.name+'.qualFiltered.bam > '+sample.dataPath+'/'+sample.name+'.noDuplicates.bam  2> '+sample.logPath+'/stderr.samtoolsView.removeDups.'+sample.name+'.txt \n'
        output += 'java -Xmx5g -jar '+self.analysispipe.settings.picardLocation+'/BuildBamIndex.jar INPUT='+sample.dataPath+'/'+sample.name+'.noDuplicates.bam '+'1>&2  2>  '+sample.logPath+'/stderr.buildIndex3.'+sample.name+'.txt \n'
        output += 'samtools flagstat '+sample.dataPath+'/'+sample.name+'.noDuplicates.bam > '+sample.logPath+'/noDuplicatesBamFlagstat.'+sample.name+'.txt \n'
        output += 'rm -v '+sample.tempPath+'/'+sample.name+'.qualFiltered.bam\n'

        #
        # Final output and write script to file
        #
        #output += '\n'+AnalysisPipe.programPath+' '+AnalysisPipe.path+' report\n'
        output += 'echo "Done. $(date) Running on: $(hostname)"\n'
        output += 'wait\n'
        output += 'echo "$(date) AllDone"\n'
        
        output_file.write(output)
        output_file.close()