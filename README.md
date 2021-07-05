Prophage Tracer
========

Prophage Tracer: Precisely tracing prophages in prokaryotic genomes using overlapping split-read alignment

Requirement
------

#### System and software requirements

1. Linux (Tested in CentOS 6.8 and CentOS Linux release 7.8.2003 (GNU Awk 4.0.2))
2. blastn: 2.6.0+ https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/
3. bwa 0.6 http://bio-bwa.sourceforge.net/
4. sambamba 0.8.0 http://lomereiter.github.io/sambamba/
5. samtools 1.10 http://www.htslib.org/

#### Other softwares may be useful for data pre-processing steps
1. Shovill 1.1.0(https://github.com/tseemann/shovill) for assembling genomes. It is useful for detect prophages assembled into their own seperate contigs in contig-level genomes. In this case, the average sequencing depth of prophage-derived contigs usually but not necessary have significantly higher depth than other contigs. The depth is written into the name of each contig in the output of Shovill.
2. Trimmomatic 0.39(https://github.com/usadellab/Trimmomatic) for remove low-quality regions and adapters in reads.

Installation
------
1. Just download the shell scripts prophage_tracer.sh and prophage_tracer.sh to your working directory of recommended linux
2. Install required softwares through Conda
####first install conda
```Bash
wget -c https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
using Enter or typing "yes", chose "no" when changning envionment variables
```
####Chaing your envionment variables temporarily to use conda
```Bash
export PATH=~/miniconda2/bin:$PATH
```
```Bash
conda install -c bioconda bwa
conda install -c bioconda bwa
conda install -c bioconda sambamba
conda install -c bioconda samtools
```
#### For Windows
If you have gerenated an SAM file using samtools, our script can be equally used on Windows 10 with Git Bash (GNU bash, version 4.4.23(1)-release (x86_64-pc-msys) with GNU Awk 5.0.0, API: 2.0 (GNU MPFR 4.1.0, GNU MP 6.2.0)) and blastn: 2.6.0+ installed . For installing Git Bash, please refer to https://git-scm.com/downloads.

Run Prophage Tracer
------

* Download the `prophage_tracer.sh` and put it in your working path
* Prepare aligned reads in SAM file

Assume that you have a sequenced bacterium (strain1) genome in FASTA foramt `reference_genome_strain1.fasta`, paried reads in FASTQ foramt `1.fastq.gz` and `2.fastq.gz`

#### Align reads to the reference genome
```Bash
bwa index reference_genome_strain1.fasta -p strain1
bwa mem strain1 1.fastq.gz 2.fastq.gz >strain1.sam
```
#### Remove duplicate (this step is not necessary but will improve the output)
```Bash
samtools view -S -b strain1.sam -o strain1.bam
sambamba markdup -r strain1.bam strain1.rmdup.bam
samtools view strain1.rmdup.bam -o strain1.rmdup.sam 
```
#### Run prophage_tracer.sh
```Bash
bash prophage_tracer.sh -m strain1.rmdup.sam -r reference_genome_strain1.fasta -p strain1
```
#### Usage

```Bash
usage:   prophage_tracer [options] -m <in.sam> -r <in.fasta> -p <prefix>

options:
     -m  FILE    a full SAM file (required)
     -r  FILE    a reference genome sequence (required)
     -p  STRING  prefix of output files (required; usually a strain name or a sample name)
     -x  INT     maximal size of a prophage (default: 150000)
     -n  INT     minimal size of a prophage (default: 5000)
     -a  INT     minimal length of attchment site (default: 3)
     -t  INT     number of threads used for BlastN (default: 1)
     -s  INT     minimal event of split reads required for supporting a prophage candidate (default: 1)
     -d  INT     minimal event of discordant read pairs required for supporting a prophage candidat (default: 1)
```

#### Typical output
prophage_candidate|contig|attL_start|attL_start|attR_end|attR_end|prophage_size|SR_evidence_attB|SR_evidence_attP|DRP_evidence_attB|DRP_evidence_attP|
|----------|-------------------------|----|----|----|----|----|----|----|----|----|
candidate_1|contig00007=::=contig00014|209162|209236|2365|2439|16770|0|4|1|2
candidate_2|contig00001|1064123|1064145|1100156|1100178|36033|0|1|0|0
candidate_3|=contig00003::=contig00004|1700|1764|46895|46959|48658|2|28|2|24

#### Explaination
Prophage Phm2 (candidate_2) locates on contig00001 (1064123-1100178) with attL (1064123-1064145) and attR (1100156-1100178). Prophage Phm1 (candidate_1) is seperated in to two or more contigs and consist at least 3'end of contig00007 (att site: 209162-209236 on contig00007) and 5' end of contig00014 (att site: 2365-2439 on contig00014). Predicted prophage size is smaller than true size. Similaryly, Prophage Phm3 (candidate_3)
is seperated in to two or more contigs and consist at least 5'end of contig00003 (att site: 1700-1764 on contig00003) and 5' end of contig00014 (att site: 46895-46959 on contig00014).



#### Notes
prophage_tracer.sh is used for chromosome-level level genoems. prophage_tracer_WGS.sh is usually used for contig-level level genoems. Although prophage_tracer_WGS.sh is suitbale for analysis of chromosome-level level genoems, it will be slow. We recommend using prophage_tracer.sh for analysis of chromosome-level level genoems to reduce running time. 

Using generate_DNA.sh for generating simulated genomes resulting from prophage excision
------

Install `seqkit` first
```Bash
conda install -c bioconda seqkit
```
Run script (default: 20 genomes containing one prophage each)
```Bash
bash generate_DNA.sh
```
