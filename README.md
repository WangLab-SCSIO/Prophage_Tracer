Prophage GPS
========

Prophage GPS: Precisely tracing prophages in prokaryotic genomes using overlapping split-read alignment

Requirement
------

#### System and software requirements

1. Linux (Tested in CentOS 6.8 and 7.0)
2. [BLAST +](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/)
3. [BWA](http://bio-bwa.sourceforge.net/)
4. [sambamba](http://lomereiter.github.io/sambamba/)
5. [samtools](http://www.htslib.org/)

Installation
------

Install through Conda
```Bash
conda uninstall -c bioconda bwa
conda uninstall -c bioconda sambamba
conda uninstall -c bioconda samtools
```
Run Prophage GPS
------

* Download the `prophage_GPS.sh` and put it in your working path
* Prepare aligned reads in SAM file

Assume that you have a sequenced bacterium (strain1) genome in FASTA foramt `reference_genome_strain1.fasta`, paried reads in FASTQ `1.fastq.gz` and `2.fastq.gz`

#### Align reads to reference genome
```Bash
bwa index reference_genome_strain1.fasta -p strain1
bwa mem strain1 1.fastq.gz 2.fastq.gz >strain1.sam
```
#### Remove duplicate reads
```Bash
samtools view -S -b strain1.sam -o strain1.bam
sambamba markdup -r strain1.bam strain1.rmdup.bam
samtools view strain1.rmdup.bam -o strain1.rmdup.sam 
```
#### Run prophage_GPS.sh
```Bash
bash prophage.sh -m strain1.rmdup.sam -r reference_genome_strain1.fasta -p strain1
```
#### Usage

```Bash
usage:   prophage_GPS [options] -m <in.sam> -r <in.fasta> -p <prefix>

options:
     -m  FILE    a full SAM file (required)
     -r  FILE    a reference genome sequence (required)
     -p  STRING  prefix of output files (required; usually a strain name or a sample name)
     -x  INT     maximal siza of prophage (default: 150000)
     -n  INT     minimal size of prophage (default: 5000)
     -a  INT     minimal length of attchment site (default: 10)
     -t  INT     number of threads used for BlastN (default: 1)
     -s  INT     minimal event of split reads required for supporting a prophage candidate
     -d  INT     minimal event of discordant read pairs required for supporting a prophage candidat
```
