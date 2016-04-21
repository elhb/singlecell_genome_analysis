import vcf
import sys

AllHet=open('/proj/b2014159/fateMapping/parsedTEST', 'w')


for variant in vcf.vcfParser('/proj/b2014159/fateMapping/results/150702.TcellClones/vcfParsing/with_bulk_indels_recalibrated.vcf'):
    #for samplename in variant.perSampleInfo:
    samplenames = sorted(variant.perSampleInfo.keys())
    if variant.perSampleInfo["1002Bulk"]["DP"] >= 10 and variant.perSampleInfo["1002Bulk"]["GQ"] >=30 and len(set(variant.perSampleInfo["1002Bulk"]['GT']))== 3 and variant.chrom:
    # Skriv all info till filen
        AllHet.write( str( variant.chrom )+'\t'+str( variant.pos ) +'\t' + str(variant.refBase) + '/' + str(variant.altBases) + '\t')
        for samplename2 in samplenames:
            AllHet.write( str( samplename2 ) +'\t')
            AllHet.write( vcf.genotypeAsBases(variant, variant.perSampleInfo[samplename2]['GT'] ) +'\t')
            AllHet.write( str( variant.perSampleInfo[samplename2]['DP'] ) +'\t')
            AllHet.write( str( variant.perSampleInfo[samplename2]['AD'] ) +'\t')
        AllHet.write( '\n' )
