# singlecell_genome_analysis

This is a pipeline for generation and submission of sbatch scripts for genomic analysis of single cells


##The generated scripts

1. **Read trimming**
This is a script containing optional WGA adapter trimming, Illumina adapter trimming, Quality trimming, removal of empty reads and FastQC
2. **Mapping**
This is a script for bowtie2 mapping
3. **Mergeing**
This is a script for merging mapped file pairs
4. **Filter and Fix**
This is a script for sorting bam files, mark duplicates, add or replace read groups and samtools flagstat
5. **Realign target creator**
This is a script for finding targets for indel realignment
6. **Realign and Calibrate**
This is a script for realigning reads around indels and quality recalibration
7. **Haplotype Calling**
This is a script for haplotype calling
8. **QC Steps**
This is a script for GATK callable loci, qa compute, picard HS metrics and generation of bed files

##Reproducing the analysis
###Note on Uppmax and Slurm
All scripts were run at Uppmax, Uppsala University's resource of high-performance computers in Uppsala, Sweden (http://www.uppmax.uu.se/). The Uppmax resource uses the SLURM workload manager (http://slurm.schedmd.com/) and the scripts are therefore configure to use this system. it also loads binaries such as bowtie2, fastqc etc into the PATH through "module load" commands.  If these type of systems is not present on you machine you should be able to run the automatically generated sbatch scripts using bash (possibly after some manual editing).

###Dependencies
To run the analysis some additional software is requiered:

1. TrimBWAstyle.pl (from http://wiki.bioinformatics.ucdavis.edu/index.php/TrimBWAstyle.pl)
1. samtools (http://samtools.sourceforge.net/)
1. picard tools (http://picard.sourceforge.net/)
1. bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
1. cutadapt (https://code.google.com/p/cutadapt/)
1. fastqc (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
1. GATK (https://www.broadinstitute.org/gatk/)

####Hardcoded variables
Some paths and filenames etc are hardcoded, if you want to reproduce the analysis you need to change the paths in the script (or rename any files downloaded from the SRA to match the scripts as well as changing the locations of executables on your machine).

Do not hesitate to contact me if you need assistance in running the scripts.
