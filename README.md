Prophage GPS
========

Prophage GPS: Precisely tracing prophages in prokaryotic genomes using overlapping split-read alignment

# Requirement
========

#### System and software requirements

1. Linux (Tested in CentOS 6.8 and 7.0)
2. [BLAST +](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/)
3. [BWA](http://bio-bwa.sourceforge.net/)
4. [sambamba](http://lomereiter.github.io/sambamba/)
5. [samtools](http://www.htslib.org/)

#### Installation

Install through Conda
```Bash
conda uninstall -c bioconda bwa
conda uninstall -c bioconda sambamba
conda uninstall -c bioconda samtools

# Run Prophage GPS

*Download the `prophage_GPS.sh` and put it in your working path
*Prepare aligned reads in SAM file
bwa index reference_genome_strain1.fasta -p strain1

bwa mem strain1 1.fastq.gz 2.fastq.gz >strain1.sam

samtools view -S -b $sam_file -o $prefix.bam
sambamba markdup -r $prefix.bam $prefix.rmdup.bam
samtools view $prefix.rmdup.bam -o $prefix.rmdup.sam 
