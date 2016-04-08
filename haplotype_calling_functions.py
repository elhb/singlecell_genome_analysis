class SampleHaplotypeCaller():
    
    def __init__(self, analysispipe):
        self.analysispipe = analysispipe # link/connection to the analysispipe object
        
        
    def haplotype_calling(self, sample):
        #for filePairId,readCount,fastq1,fastq2,sampleId in sample.getFastqs():
        output_file = open(sample.scriptPath+'/haplotypeCalling.'+sample.name+'.sbatch.sh', 'w')
            
        #
        # sbatch header
        #
        output = '#! /bin/bash -l'+'\n'
        output += '#SBATCH -A '+self.analysispipe.settings.uppmaxProject+'\n'
        output += '#SBATCH -n 16 -p node'+'\n'
        output += '#SBATCH -t 240:00:00'+'\n'
        output += '#SBATCH -J haplotypeCalling.'+sample.name+'\n'
        output += '#SBATCH -e '+sample.logPath+'/stderr.haplotypeCalling.'+sample.name+'.txt'+'\n'
        output += '#SBATCH -o '+sample.logPath+'/stdout.haplotypeCalling.'+sample.name+'.txt'+'\n'

        #
        # define variebles and go to path
        #
        output += 'echo "$(date) Running on: $(hostname)"'+'\n'
        output += 'cd '+sample.path+'\n'
        output += 'echo'+'\n'
        output += 'module load bioinfo-tools pysam/0.8.3-py27  FastQC cutadapt/1.8.0 bowtie2 samtools picard/1.114 BEDTools/2.16.2  GATK/3.1.1'+'\n'
        output += 'module load python/2.7'+'\n'
        
        output += 'echo "Haplotype Calling" '+'\n'
        
        output += 'cp -v '+sample.dataPath+'/'+sample.name+'.noDuplicates.bam $SNIC_TMP\n'
        output += 'cp -v '+sample.dataPath+'/'+sample.name+'.noDuplicates.bai $SNIC_TMP\n'
        output += 'java -Xmx50g -jar '+self.analysispipe.settings.GATKlocation+' '
        output += '-T HaplotypeCaller -nct 16 '
        output += '-R '+self.analysispipe.settings.bowtie2Reference+' '
        output += '-I $SNIC_TMP'+'/'+sample.name+'.noDuplicates.bam '
        output += '--genotyping_mode DISCOVERY '
        output += '-stand_emit_conf 10 '
        output += '-stand_call_conf 30 '
        
        #
        # Keep these?
        #
        #if wgsOrExome == 'exome': output += '-L '+self.AnalysisPipe.referencePath+'/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem.bed '
        #output += '-L '+self.analysispipe.settings.targets+'/wgs.bed '
        
        output += '--dbsnp '+self.analysispipe.settings.GATKbundleLocation+'/dbsnp_138.b37.vcf '
        output += '--annotation AlleleBalance --annotation AlleleBalanceBySample --annotation BaseCounts --annotation BaseQualityRankSumTest '
        output += '--annotation ChromosomeCounts --annotation ClippingRankSumTest --annotation Coverage --annotation DepthPerAlleleBySample '
        output += '--annotation DepthPerSampleHC --annotation FisherStrand --annotation GCContent --annotation HaplotypeScore --annotation HardyWeinberg '
        output += '--annotation HomopolymerRun --annotation InbreedingCoeff --annotation LikelihoodRankSumTest --annotation LowMQ '
        output += '--annotation MVLikelihoodRatio --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation MappingQualityZeroBySample '
        output += '--annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation ReadPosRankSumTest --annotation SampleList '
        output += '--annotation SnpEff --annotation SpanningDeletions --annotation StrandBiasBySample --annotation TandemRepeatAnnotator '
        output += '--annotation TransmissionDisequilibriumTest --annotation VariantType '#--annotation StrandOddsRatio 
        output += '--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 '
        output += '-o '+sample.dataPath+'/'+sample.name+'.gvcf '
        output += '1>&2 2> '+sample.logPath+'/stderr.haplotypeCallerGatk.'+sample.name+'.txt &'+'\n'
        output += 'wait'+'\n'

        #
        # Final output and write script to file
        #
        output += 'echo "Done. $(date) Running on: $(hostname)"\n'
        output += 'wait\n'
        output += 'echo "$(date) AllDone"\n'
        
        output_file.write(output)
        output_file.close()