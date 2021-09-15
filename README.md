Prophage Tracer
========

Prophage Tracer: Precisely tracing prophages in prokaryotic genomes using overlapping split-read alignment.
Based on our analysis, in order to detect prophages with low excision rate, 100–1000× sequencing depth for a genome is recommended. At this range of sequencing depth, Prophage Tracer can detect the hidden prophages with excision rates (attB/gyrB) > $10^-3$ and/or replication (attP/gyrB) > $10^-3$ in host genomes.

Requirement
------

#### System and software requirements

1. Linux (Tested in CentOS 6.8 and CentOS Linux release 7.8.2003 (GNU Awk 4.0.2))
2. [blastn: 2.6.0+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/), Download ncbi-blast-2.6.0+-x64-linux.tar.gz for linux system.
3. [bwa 0.6](http://bio-bwa.sourceforge.net/)
4. [sambamba 0.8.0](http://lomereiter.github.io/sambamba/)
5. [samtools 1.10](http://www.htslib.org/)

#### Other softwares may be useful for data pre-processing steps
1. [Shovill 1.1.0](https://github.com/tseemann/shovill) for assembling genomes. It is useful for detecting prophages assembled into their own seperate contigs in contig-level genomes. In this case, the average sequencing depth of prophage-derived contigs is usually but not necessarily significantly higher than other contigs. The depth is written into the name of each contig in the output of Shovill.
2. [Trimmomatic 0.39](https://github.com/usadellab/Trimmomatic) for removeing low-quality regions and adapters in reads.

Installation
------
1. Just download the shell scripts prophage_tracer.sh and prophage_tracer.sh to your working directory
2. Install required softwares through Conda

#### first install conda
```Bash
wget -c https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
```
#### Chaing your envionment variables temporarily to use conda
```Bash
export PATH=~/miniconda2/bin:$PATH
```
```Bash
conda install -c bioconda bwa
conda install -c bioconda sambamba
conda install -c bioconda samtools
```

Run Prophage Tracer
------
#### 1. Download the `prophage_tracer.sh` and put it in your working path

#### 2. Prepare aligned reads in SAM file

Assume that you have a sequenced bacterium (strain1) genome in FASTA format `reference_genome_strain1.fasta`, and paried reads in FASTQ format `1.fastq.gz` and `2.fastq.gz`

* Align reads to the reference genome
```Bash
bwa index reference_genome_strain1.fasta -p strain1
bwa mem strain1 1.fastq.gz 2.fastq.gz >strain1.sam
```
* Remove duplicate (this step is not necessary but will improve the output)
```Bash
samtools view -S -b strain1.sam -o strain1.bam
sambamba markdup -r strain1.bam strain1.rmdup.bam
samtools view strain1.rmdup.bam -o strain1.rmdup.sam 
```
#### 3. Run prophage_tracer.sh
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
     -a  INT     minimal length of attchment site (default: > 2)
     -t  INT     number of threads used for BlastN (default: 1)
     -s  INT     minimal event of split reads required for supporting a prophage candidate (default: 1)
     -d  INT     minimal event of discordant read pairs required for supporting a prophage candidat (default: 1)
```

#### Typical output
*Find result in `strain1.prophage.out`
|prophage_candidate|contig|attL_start|attL_end|attR_start|attR_end|prophage_size|SR_evidence_attB|SR_evidence_attP|DRP_evidence_attB|DRP_evidence_attP|
|----------|-------------------------|----|----|----|----|----|----|----|----|----|
|candidate_1|contig00007=::=contig00014|209162|209236|2365|2439|16770|0|4|1|2
|candidate_2|contig00001|1064123|1064145|1100156|1100178|36033|0|1|0|0
|candidate_3|=contig00003::=contig00004|1700|1764|46895|46959|48658|2|28|2|24
* Find 

#### Explanation
1. If a single contig was given in the `contig` column, it means an intact predicted prophage is in this contig.
2. "::" indicates a predicted prophage is seperated into two or more contigs. "=" indicates that the 5' end or 3' end of a contig contains a part of a prophage. "=contig00003" indicates the 5' end of contig00003 while "contig00003=" indicates the 3' end of contig00003.
3. "SR_evidence_attB/attP" indicate split read counts support *attB* or *attP* of the predicted prophage. "DRP_evidence_attB/attP" indicate discordant read pair counts support *attB* or *attP* of the predicted prophage.
4. Prophage Phm2 (candidate_2) locates on contig00001 (1064123-1100178) with *attL* (1064123-1064145) and *attR* (1100156-1100178). Prophage Phm1 (candidate_1) is seperated into two or more contigs and consists at least 3' end of contig00007 (*att* site: 209162-209236 on contig00007) and 5' end of contig00014 (*att* site: 2365-2439 on contig00014). Prophage Phm3 (candidate_3)
is seperated into two or more contigs and consist at least 5' end of contig00003 (*att* site: 1700-1764 on contig00003) and 5' end of contig00014 (*att* site: 46895-46959 on contig00014).

#### Notes
1. `prophage_tracer.sh` is used for chromosome-level genomes. `prophage_tracer_WGS.sh` can be used for chromosome-level and contig-level genomes. However, using prophage_tracer_WGS.sh for analysis of chromosome-level  genomes would be slow.
2. If a prophage is seperated into two or more contigs, the predicted prophage size might be smaller than the true size.
3. If you have gerenated a SAM file using samtools, our script can be equally used to predict propahges on Windows 10 with Git Bash and blastn (`ncbi-blast-2.6.0+-win64.exe`) installed. For installing Git Bash and blastn in Windows 10, please refer to https://git-scm.com/downloads and https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/.

Using `generate_DNA.sh` for generating simulated genomes resulted from prophage excision
------

Install `seqkit` first
```Bash
conda install -c bioconda seqkit
```
Download ``generate_DNA.sh`` and ``random_DNA.py``

Run script (default: simulating 20 genomes; one prophage in each genome)
```Bash
bash generate_DNA.sh
```


Coming updates
------
1. Addding a parameter of `blastn` to set whether mismatch were allowed in the *att* sites (-penalty).
2. To make the script to be able to analyse single read sequencing data.

Copyright
------
Kaihao Tang, khtang@scsio.ac.cn;
Xiaoxue Wang, xxwang@scsio.ac.cn;
Marine Biofilm Lab;
SCSIO, Chinese Academy of Sciences