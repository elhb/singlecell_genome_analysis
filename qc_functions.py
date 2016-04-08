class SampleQC():
    
    def __init__(self, analysispipe):
        self.analysispipe = analysispipe # link/connection to the analysispipe object
        
        
    def qcSteps(self, sample):
        
        import sys
        
        #for filePairId,readCount,fastq1,fastq2,sampleId in sample.getFastqs():
        output_file = open(sample.scriptPath+'/qcSteps.'+sample.name+'.sbatch.sh', 'w')
            
        #
        # sbatch header
        #
        output = '#! /bin/bash -l'+'\n'
        output += '#SBATCH -A '+self.analysispipe.settings.uppmaxProject+'\n'
        output += '#SBATCH -n 1 -p core'+'\n'
        output += '#SBATCH -t 240:00:00'+'\n'
        output += '#SBATCH -J qcSteps.'+sample.name+'\n'
        output += '#SBATCH -e '+sample.logPath+'/stderr.qcSteps.'+sample.name+'.txt'+'\n'
        output += '#SBATCH -o '+sample.logPath+'/stdout.qcSteps.'+sample.name+'.txt'+'\n'

        #
        # define variebles and go to path
        #
        output += 'echo "$(date) Running on: $(hostname)"'+'\n'
        output += 'cd '+sample.path+'\n'
        output += 'echo'+'\n'
        output += 'module load bioinfo-tools pysam/0.8.3-py27  FastQC cutadapt/1.8.0 bowtie2 samtools picard/1.114 BEDTools/2.16.2  GATK/3.1.1'+'\n'
        output += 'module load python/2.7'+'\n'
        
        #
        # GATK callable Loci
        #
        output += 'echo "$(date) Running on: $(hostname)"'+'\n'
        output += 'echo -e "-> CallableLoci <-"'+'\n'
        output += 'java -Xmx5g -jar '+self.analysispipe.settings.GATKlocation+' -T CallableLoci '
        output +='-I '+sample.dataPath+'/'+sample.name+'.noDuplicates.bam '
        output +='-summary '+sample.dataPath+'/'+sample.name+'.callableLociSummary.txt '
        output +='-o '+sample.dataPath+'/'+sample.name+'.callableLoci.bed '
        output +='-R '+self.analysispipe.settings.GATK_reference+' '+'\n'
        output += 'echo "Done. $(date) Running on: $(hostname)"'+'\n'
        output += 'echo'+'\n'
        output += 'echo "-----"'+'\n'

        #
        # qacompute
        #
        output += 'echo "$(date) Running on: $(hostname)"'+'\n'
        output += 'echo -e "-> Pauls qacompute <-"'+'\n'
        output += '/proj/b2010052/scripts/qaCompute -d -q 10 '
        output += '-m '+sample.dataPath+'/'+sample.name+'.noDuplicates.bam '
        output += sample.dataPath+'/'+sample.name+'.qacompute.out '
        output += '> '+sample.logPath+'/'+sample.name+'.qacompute.stdout.txt '
        output += '2> '+sample.logPath+'/'+sample.name+'.qacompute.stderr.txt '+'\n'
        output += 'echo "Done. $(date) Running on: $(hostname)"'+'\n'

        #
        # picard HS metrics
        #
        output += 'java -Xmx3g -jar '+self.analysispipe.settings.picardLocation+'/CalculateHsMetrics.jar '
        
        #if wgsOrExome == 'exome':output += 'BAIT_INTERVALS='  +self.analysispipe.settings.targets+'/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem '
        #output += 'BAIT_INTERVALS='  +self.analysispipe.settings.targets+'/wgs '
        #if wgsOrExome == 'exome':output += 'TARGET_INTERVALS='+self.analysispipe.settings.targets+'/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem '
        #output += 'TARGET_INTERVALS='+self.analysispipe.settings.targets+'/wgs '
        
        output += 'INPUT='+sample.dataPath+'/'+sample.name+'.noDuplicates.bam '
        output += 'OUTPUT='+sample.dataPath+'/'+sample.name+'.hs_metrics.summary.txt '
        output += 'PER_TARGET_COVERAGE='+sample.dataPath+'/'+sample.name+'.hs_metrics.perTargetCoverage.txt '
        output += 'REFERENCE_SEQUENCE='+self.analysispipe.settings.bowtie2Reference+'  '
        output += '1>&2 2> '+sample.logPath+'/'+sample.name+'.stderr.caluclateHsmetrics.txt \n'
        
        #
        #make files for coverage checks
        #
        #if wgsOrExome == 'exome':
            #output += "bedtools coverage -abam "+sample.dataPath+'/'+sample.name+'.noDuplicates.bam'+" -b "+self.analysispipe.settings.targets+"/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem.bed -d | sort -k 1,1 -k2,2n -k3,6n > "+sample.dataPath+'/'+sample.name+'.bedtools.coverage.bed\n'
            #output += "bedtools coverage -abam "+sample.dataPath+'/'+sample.name+'.noDuplicates.bam'+" -b "+self.analysispipe.settings.targets+"/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem.bed -hist > "+sample.dataPath+'/'+sample.name+'.bedtools.exomecoverage.histogram\n'
        output += "bedtools genomecov -ibam "+sample.dataPath+'/'+sample.name+'.noDuplicates.bam'+" > "+sample.dataPath+'/'+sample.name+'.bedtools.genomecov.histogram\n'
        output += "bedtools genomecov -bga -ibam "+sample.dataPath+'/'+sample.name+'.noDuplicates.bam'+" > "+sample.dataPath+'/'+sample.name+'.bedtools.genomecov.bedgraph\n'
    
        #if wgsOrExome == 'exome':
            #output += "bedtools coverage -split -abam "+sample.dataPath+'/'+sample.name+'.noDuplicates.bam'+" -b "+self.analysispipe.settings.targets+"/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem.bed -d | sort -k 1,1 -k2,2n -k3,6n > "+sample.dataPath+'/'+sample.name+'.bedtools.coverage.nonPhysical.bed\n'
            #output += "bedtools coverage -split -abam "+sample.dataPath+'/'+sample.name+'.noDuplicates.bam'+" -b "+self.analysispipe.settings.targets+"/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem.bed -hist > "+sample.dataPath+'/'+sample.name+'.bedtools.exomecoverage.nonPhysical.histogram\n'
        output += "bedtools genomecov -split -ibam "+sample.dataPath+'/'+sample.name+'.noDuplicates.bam'+" > "+sample.dataPath+'/'+sample.name+'.bedtools.genomecov.nonPhysical.histogram\n'
        output += "bedtools genomecov -split -bga -ibam "+sample.dataPath+'/'+sample.name+'.noDuplicates.bam'+" > "+sample.dataPath+'/'+sample.name+'.bedtools.genomecov.nonPhysical.bedgraph\n'
    
        output += "awk '{print $7}' "+sample.dataPath+'/'+sample.name+'.bedtools.coverage.bed'+" | sort -n | uniq -c | awk '{print $1\"\\t\"$2}'> "+sample.dataPath+'/'+sample.name+".coverageDistribution.tsv\n"
        output += "awk '{print $1 \"\\t\" $2+$6-1 \"\\t\" $7}' "+sample.dataPath+'/'+sample.name+'.bedtools.coverage.bed'+" > "+sample.dataPath+'/'+sample.name+".depthPerPosition.tsv\n"
    
        #if wgsOrExome == 'exome':
            #output += 'gzip -v9 '+sample.dataPath+'/'+sample.name+'.bedtools.coverage.bed\n'
            #output += 'gzip -v9 '+sample.dataPath+'/'+sample.name+'.bedtools.exomecoverage.histogram\n'
            #output += 'gzip -v9 '+sample.dataPath+'/'+sample.name+'.bedtools.coverage.nonPhysical.bed\n'
            #output += 'gzip -v9 '+sample.dataPath+'/'+sample.name+'.bedtools.exomecoverage.nonPhysical.histogram\n'
        output += 'gzip -v9 '+sample.dataPath+'/'+sample.name+'.bedtools.genomecov.histogram\n'
        output += 'gzip -v9 '+sample.dataPath+'/'+sample.name+'.bedtools.genomecov.bedgraph\n'
        output += 'gzip -v9 '+sample.dataPath+'/'+sample.name+'.bedtools.genomecov.nonPhysical.histogram\n'
        output += 'gzip -v9 '+sample.dataPath+'/'+sample.name+'.bedtools.genomecov.nonPhysical.bedgraph\n'
        output += 'gzip -v9 '+sample.dataPath+'/'+sample.name+".coverageDistribution.tsv\n"
        output += 'gzip -v9 '+sample.dataPath+'/'+sample.name+".depthPerPosition.tsv\n"
        
        #
        # create bed for cnv stuff
        #
        output += "bamToBed -i "+sample.dataPath+'/'+sample.name+'.noDuplicates.bam'+" > "+sample.dataPath+'/'+sample.name+'.bed'+"\n"
        output += "gzip -v9 "+sample.dataPath+'/'+sample.name+'.bed'+"\n"
    
        #
        # Final output and write script to file
        #
        output += '\n'+self.analysispipe.settings.program_path+' '+self.analysispipe.path+' report\n'
        output += 'echo'+'\n'
        output += 'wait'+'\n'
        output += 'echo "$(date) AllDone"'+'\n'
        output += 'echo "$(date) AllDone" >&2'+'\n'
        
        output_file.write(output)
        output_file.close()