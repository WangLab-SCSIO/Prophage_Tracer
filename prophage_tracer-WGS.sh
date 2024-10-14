#!/bin/bash -e 

set -eo pipefail

##########################################################
# Prophage Tracer V1.0.2
#########################################################

## usage
usage() {
    echo "
Prophage Tracer V1.0.2 07/05/2021
usage:   bash prophage_tracer-WGS.sh [options] -m <in.sam> -r <in.fasta> -p <prefix>
requirement：locally installed BLAST+ software

options:
     -m  FILE    a full SAM file (required)
     -r  FILE    a reference genome sequence (required)
     -p  STRING  prefix of output files (required; usually a strain name or a sample name)
     -x  INT     maximal size of prophage (default: 150000)
     -n  INT     minimal size of prophage (default: 5000)
     -a  INT     minimal length of attchment site (default: 10)
     -t  INT     number of threads used for BlastN (default: 1)
	 -s  INT     minimal event of split reads required for supporting a prophage candidate (default: 1)
	 -d  INT     minimal event of discordant read pairs required for supporting a prophage candidate (default: 1)
"
}

while getopts ":hm:r:p:x:n:a:t:s:d:" opt
do
  case $opt in
    h) usage
	   exit 0
	   ;;
    m) sam_file=$OPTARG ;;
    r) fasta=$OPTARG ;;
    p) prefix=$OPTARG ;;
	x) prophage_size_max=$OPTARG ;;
    n) prophage_size_min=$OPTARG ;;
	a) att_length_min=$OPTARG ;;
    t) threads=$OPTARG ;;
    s) SR_COUNT=$OPTARG ;;
    d) DRP_COUNT=$OPTARG ;;
    :)
      echo "Option -$OPTARG requires an argument."
	  usage
      exit 1
      ;;
	?) 
	  echo "Invalid option: -$OPTARG"
	  usage      
      exit 1
	  ;;
  esac
done

#If options of 'x', 'n', 'a' and 't' are not defined, setting to default values
prophage_size_max=${prophage_size_max:-150000}
prophage_size_min=${prophage_size_min:-5000}
att_length_min=${att_length_min:-3}
threads=${threads:-1}
SR_COUNT=${SR_COUNT:-1}
DRP_COUNT=${DRP_COUNT:-1}

#Check locally installed BLAST+ software
BLASTN=`which blastn || true`
if [[ -z "$BLASTN" ]]
then
    usage
	echo "blastn not found"
	exit 1
fi

###########统计contig的长度
cat $fasta | awk '$0 ~ ">" {print c; c=0;split($0,a,"[> ]"); printf a[2] "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > contiglength.file


#Get read length and insert size
#read_length=`head -n 10000 $sam_file | gawk 'BEGIN { max=0 } { if (length($10)>max) max=length($10) } END { print max }'`
read_length=150
#insert_size=`head -n 10000 $sam_file | awk '($2 ~ /163|83|99|147/ )' | cut -f9 | awk '{print sqrt($0^2)}' | awk '$0<10000'| awk '{ sum += $0;i++ } END { print int(sum/i) }'`
insert_size=500

#Step 1: Extracting split reads
echo -e "Step 1: Extracting split reads"

awk '($1 ~ /^@/)||($2 ~ /145|81|99|163/ && $6 ~ /^...?S...?M$/) ||($2 ~ /113|177/ && $6 ~ /^...?S...?M$/ && $9 == 0 ) || ($2 ~ /97|161|147|83/  && $6 ~ /^...?M...?S$/ ) || ($2 ~ /65|129/ && $6 ~ /^...?M...?S$/ && $9 == 0 )' $sam_file > $prefix.sr.temp.1

#Extracting split reads.If paried reads are overlapping, only one kept.

awk '!a[$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8]++' $prefix.sr.temp.1 | awk '($1 !~ /^@/)'| awk '!a[$1]++'| awk '{printf ">"$1"_"$2"\n"$10"\n"}' > $prefix.reads.fasta

#Blastn searching reads against the reference genome.

makeblastdb -in $fasta -dbtype nucl -out $prefix.nuclDB > makeblastdb.log 2>&1 

blastn -query $prefix.reads.fasta -db $prefix.nuclDB -out $prefix.blastn.result -evalue 1e-3 -outfmt 6 -word_size 11 -num_threads $threads > blastn.log 2>&1 

#Manipulate blastn result ######################################################################################下面筛选到2条后再加一个按照reads起始位点排序的命令

awk 'BEGIN{FS="\t";OFS="\t"} {if ($3 > 90) a[$1,++b[$1]]=$0}END{for(i in b) if (b[i] == 2) print a[i,1]"\n"a[i,2]}' $prefix.blastn.result| awk '{if($0!="") print}' |sort -k1,1 -k7n,7 > $prefix.sr.temp.2

#Classify reads containing attB or attP

awk '

BEGIN {FS="\t";OFS="\t"}

function calculation1() {
     attL_start=send1-qend1+qstart2
     attL_end=send1
     attR_start=sstart2
     attR_end=sstart2+qend1-qstart2
     prophage_size=attR_end-attL_end
     print SR_type,SR_DRP_type,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$0
}

function calculation2() {
     attL_start=send2-qend2+qstart1
     attL_end=send2
     attR_start=sstart1
     attR_end=sstart1+qend2-qstart1
     prophage_size=attR_end-attL_end
	 print SR_type,SR_DRP_type,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$0
}



function calculation3() {
     attL_start=sstart2
     attL_end=sstart2+qend1-qstart2
     attR_start=send1-qend1+qstart2
     attR_end=send1
     prophage_size=attR_end-attL_end
	 print SR_type,SR_DRP_type,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$0
}


function calculation4() {
     attL_start=sstart1
     attL_end=sstart1+qend2-qstart1
     attR_start=send2-qend2+qstart1
     attR_end=send2
     prophage_size=attR_end-attL_end
	 print SR_type,SR_DRP_type,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$0
}


{
if (NR==1) 
{
    query = $1; contig=$2; qstart1 = $7; qend1 = $8; sstart1 = $9; send1 = $10; next
}

##make sure to extract correct flag value if "_" was used in read names
Flag_pos=split ($1,flag,"_")

if ($1 == query && $2 == contig)
{   qstart2 = $7; qend2 = $8; sstart2 = $9; send2 = $10 
	
	if (flag[Flag_pos] == 83)
	{	
#attB
		if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_1"
            evidence="attB"
            calculation1()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_1"
            evidence="attB"
			calculation2()
		}
#attP	
	    if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_5"
            evidence="attP"
			calculation3()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_5"
            evidence="attP"
			calculation4()
	    }
	}

	if (flag[Flag_pos] == 147)
	{	
#attB
		if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_2"
            evidence="attB"
            calculation1()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_2"
            evidence="attB"
			calculation2()
		}
#attP	
	    if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_6"
            evidence="attP"
			calculation3()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_6"
            evidence="attP"
			calculation4()
	    }
	}

	if (flag[Flag_pos] == 97)
	{	
#attB
		if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_1"
            evidence="attB"
            calculation1()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2< 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_1"
            evidence="attB"
			calculation2()
		}
#attP	
	    if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_5"
            evidence="attP"
			calculation3()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_5"
            evidence="attP"
			calculation4()
	    }
	}
	
	if (flag[Flag_pos] == 161)
	{	
#attB
		if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_2"
            evidence="attB"
            calculation1()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_2"
            evidence="attB"
			calculation2()
		}
#attP	
	    if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_6"
            evidence="attP"
			calculation3()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_6"
            evidence="attP"
			calculation4()
	    }
	}

	
####Second

	if (flag[Flag_pos] == 99)
	{	
#attB
		if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_7"
            evidence="attB"
            calculation1()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_7"
            evidence="attB"
			calculation2()
		}
#attP	
	    if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_3"
            evidence="attP"
			calculation3()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_3"
            evidence="attP"
			calculation4()
	    }
	}

	if (flag[Flag_pos] == 163)
	{	
#attB
		if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_8"
            evidence="attB"
            calculation1()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_8"
            evidence="attB"
			calculation2()
		}
#attP	
	    if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_4"
            evidence="attP"
			calculation3()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_PR_4"
            evidence="attP"
			calculation4()
	    }
	}

	if (flag[Flag_pos] == 81)
	{	
#attB
		if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_7"
            evidence="attB"
            calculation1()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_7"
            evidence="attB"
			calculation2()
		}
#attP	
	    if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_3"
            evidence="attP"
			calculation3()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_3"
            evidence="attP"
			calculation4()
	    }
	}
	
	if (flag[Flag_pos] == 145)
	{	
#attB
		if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_8"
            evidence="attB"
            calculation1()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_8"
            evidence="attB"
			calculation2()
		}
#attP	
	    if ( (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10) && (send1 > sstart1 && sstart1 > send2 && send2 > sstart2) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_4"
            evidence="attP"
			calculation3()
		}
		else if ( (qend1 > qend2 && qend2 > qstart1 && qstart1 > qstart2 && qstart2 < 10) && (send2 > sstart2 && sstart2 > send1 && send1 > sstart1) )
		{
		    SR_type="SR"
            SR_DRP_type="SR_DRP_4"
            evidence="attP"
			calculation4()
	    }
	}
	

}
else if ($1 != query)
{
	query = $1; contig = $2; qstart1 = $7; qend1 = $8; sstart1 = $9; send1 = $10   
}
	
}' $prefix.sr.temp.2 | awk '

#Cluster split reads
 
BEGIN{FS="\t";OFS="\t"} 

{
    if( $8 > '"$prophage_size_min"' && $8 < '"$prophage_size_max"' ) a[$9]=$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10
}

END{

count=1

for (i in a) 
{
	split(a[i],b,"\t") 
	print i,b[1],b[2],b[3],b[4],b[5],b[6],b[7],"candidate_"count

    for (j in a)
	{
		split(a[j],c,"\t")  
	    if ( i != j && c[7] == b[7])
		{
		    if (b[7]==c[7] && b[2]-c[2] > -5 && b[2]-c[2] < 5 && b[4]-c[4] > -5 && b[4]-c[4] < 5 && b[6]-c[6] > -10 && b[6]-c[6] < 10)
			{
				print j,c[1],c[2],c[3],c[4],c[5],c[6],c[7],"candidate_"count
				delete a[j]				
			}	
        }
	}
    count++		
    delete a[i]
}
}' |awk '{if (!($3 == "")) print $0}' >$prefix.sr.temp.3


awk '

BEGIN{FS="\t";OFS="\t"} 
{
if ( !(a[$NF,1] >0) )
{
    a[$NF,1]=0
}

if ( !(a[$NF,2] >0) )
{
    a[$NF,2]=0
}

if ($2 == "attB")
{
    a[$NF,1]+=1 
	sum[$NF,3]+=$3
	sum[$NF,4]+=$4
	sum[$NF,5]+=$5
	sum[$NF,6]+=$6
	sum[$NF,7]+=$7
	a[$NF,8]=$8
}
else if ($2 == "attP")
{
    a[$NF,2]+=1
	sum[$NF,3]+=$3
	sum[$NF,4]+=$4
	sum[$NF,5]+=$5
	sum[$NF,6]+=$6
	sum[$NF,7]+=$7
	a[$NF,8]=$8
}
}
END{
    for(i in a) 
	{
	    split(i,idx,SUBSEP)
        a[idx[1],3]=int(sum[idx[1],3]/(a[idx[1],1]+a[idx[1],2]))
	    a[idx[1],4]=int(sum[idx[1],4]/(a[idx[1],1]+a[idx[1],2]))
	    a[idx[1],5]=int(sum[idx[1],5]/(a[idx[1],1]+a[idx[1],2]))
	    a[idx[1],6]=int(sum[idx[1],6]/(a[idx[1],1]+a[idx[1],2]))
	    a[idx[1],7]=int(sum[idx[1],7]/(a[idx[1],1]+a[idx[1],2]))
		print idx[1],a[idx[1],1],a[idx[1],2],a[idx[1],1]+a[idx[1],2],a[idx[1],3],a[idx[1],4],a[idx[1],5],a[idx[1],6],a[idx[1],7],a[idx[1],8]
	}
}

'  $prefix.sr.temp.3 | awk '!a[$0]++' > $prefix.sr.temp.out

if [ $? -eq 0 ]
then
	echo -e "Step 2: Extracting discordant read pairs "
else
    echo "!!!Unsuccessfully summarizing informations of split reads!!!" 
	exit
fi


#Step 2: Extracting discordant read pairs

awk --re-interval '$1 ~ /^@/ || ($2 ~ /97|145|81|161/ && ($6 ~ /^[1][0-5][0-9]M/ || $6 ~ /^[6-9][0-9]M/) && $7 == "=" )' $sam_file | awk '!a[$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8]++' | sort -k1,1 -k2n,2 | awk '

BEGIN {FS="\t";OFS="\t"}

{
if (NR==1) 
{
    read = $1; read_flag = $2; len=$9; a[$1] = $0; next
}


if ($1 == read)
{
    if ((read_flag == 97 && len < 0 && $2 ==145 && $9 >0) || (read_flag == 145 && len > 0 && $2 == 97 && $9 <0) || (read_flag == 161 && len < 0 && $2 == 81 && $9 >0) || (read_flag == 81 && len > 0 && $2 == 161 && $9 <0))
	{
	    print "DRP","attP",a[$1]"\n""DRP","attP",$0
	}
	
	else if ((read_flag == 145 && len < 0 && $2 ==97 && $9 >0) || (read_flag == 97 && len > 0 && $2 == 145 && $9 <0) || (read_flag == 81 && len < 0 && $2 == 161 && $9 >0) || (read_flag == 161 && len > 0 && $2 == 81 && $9 <0))
	{
	    print "DRP","attB",a[$1]"\n""DRP","attB",$0
	}
}	

else if ($1 != read)
{
	read = $1; read_flag = $2; len=$9; a[$1] = $0; next   
}
	
}' | awk '($4 ~ /97|161|81|145/ && $11 > 0)' | awk 'BEGIN {FS="\t";OFS="\t"} NR==FNR {split($1,readname,"_"); a[readname[1]]=readname[1]} NR>FNR  { if(!($3 in a)) print $0}' $prefix.sr.temp.2 - >$prefix.drp.temp.1

##First round merging DRP reads accroding to the SR position. Not merged DRP reads collected in the '$prefix'.drp.temp.left file
awk ' BEGIN {FS="\t";OFS="\t"}

FILENAME==ARGV[1] { a[$3]=$0; next} 
FILENAME==ARGV[2] {  

    T_L=int('"$insert_size"'*1.6) #threshold low for clsutering
    T_H=int('"$insert_size"'*3.5) #threshold low for clsutering
 
	attP_count = 0
	attB_count = 0
	for (i in a)
	{
	   split(a[i],b,"\t")
	   if ($6-b[6] > -T_L && $6-b[6] <T_L && $7-b[10] >-T_L && $7-b[10] <T_L && $9-b[11] > -T_H && $9-b[11] < T_H && b[5] == $10)
	   {
	       if (b[2] == "attP" && b[6] >= $5 && $8-b[10]+10 > '"$read_length"' )
		   {
		       attP_count++  
		       delete a[i]
		   }
		   else if (b[2] == "attB" && $6-b[6]+10 > '"$read_length"'  && b[10] >= $7)
		   {
		       attB_count++ 	   
		       delete a[i]	   
		   }
		}
	}
    c[$1]=$0"\t"attB_count"\t"attP_count"\t"attB_count+attP_count
	}
END	{
	for (i in c)
    {
	    print c[i]
	}
	for (i in a)
    {
	    print a[i] >> "'$prefix'.drp.temp.left"
	}
}
' $prefix.drp.temp.1 $prefix.sr.temp.out > $prefix.sr-drp.temp.out

##Cluster the left DRP reads independent of SR reads:  first candidate assign
if test -s $prefix.drp.temp.left; then

awk '

BEGIN{FS="\t";OFS="\t"} 

{
    if($11 > '"$prophage_size_min"' && $11 < '"$prophage_size_max"') a[$3]=$2"\t"$4"\t"$6"\t"$10"\t"$11"\t"$5
}

END{
count=1
T_L=int('"$insert_size"'*2.2) #threshold low for clsutering
T_H=int('"$insert_size"'*3.5) #threshold low for clsutering
for (i in a) 
{
	split(a[i],b,"\t") 
	print i,b[1],b[2],b[3],b[4],b[5],b[6],"DRP_candidate_"count

    for (j in a)
	{	
	    split(a[j],c,"\t")  
	    if ( i != j )
		{
		    if (b[6]==c[6] && b[3]-c[3] > -T_L && b[3]-c[3] < T_L && b[4]-c[4] > -T_L && b[4]-c[4] < T_L && b[5]-c[5] > -T_H && b[5]-c[5] < T_H)
			{
			    if ((b[1] == "attB" && c[1] == "attB")||(b[1] == "attP" && c[1] == "attP")||(b[1] == "attB" && c[1] == "attP" && b[3] < c[3] && b[4] > c[4])||(b[1] == "attP" && c[1] == "attB" && b[3] > c[3] && b[4] < c[4]))
				    { 
					    print j,c[1],c[2],c[3],c[4],c[5],c[6],"DRP_candidate_"count
				        delete a[j]
                	}		
			}	
        }
	}
    count++		
    delete a[i]
}
}' $prefix.drp.temp.left >$prefix.drp.temp.2

else

touch $prefix.drp.temp.2

fi

##Cluster the left DRP reads independent of SR reads:  merging assigned candidate

if test -s $prefix.drp.temp.2; then

awk '{if (!($3 == "")) print $0}' $prefix.drp.temp.2| awk '

BEGIN{FS="\t";OFS="\t"} 
{

#if (NR==1) 
#{
#    a[$NF,3] = $4+0; a[$NF,4] = $5+0
#}

a[$NF,4]=10000000000

if ( !(a[$NF,1] >0) )
{
    a[$NF,1]=0
}

if ( !(a[$NF,2] >0) )
{
    a[$NF,2]=0
}

if ($2 == "attB")
{
    a[$NF,1]+=1 
	a[$NF,3]=a[$NF,3]>=$4?a[$NF,3]:$4
	a[$NF,4]=a[$NF,4]<=$5?a[$NF,4]:$5
	a[$NF,5]=$7
}
else if ($2 == "attP")
{
    a[$NF,2]+=1
	a[$NF,3]=a[$NF,3]>=$4?a[$NF,3]:$4
	a[$NF,4]=a[$NF,4]<=$5?a[$NF,4]:$5
    a[$NF,5]=$7
}
}
END{
    for(i in a) 
	{
	    split(i,idx,SUBSEP) 
		print idx[1],a[idx[1],1],a[idx[1],2],a[idx[1],1]+a[idx[1],2],a[idx[1],3],a[idx[1],4],a[idx[1],4]-a[idx[1],3],a[idx[1],5]
	}
}

' | awk '!a[$0]++' > $prefix.drp.temp.out

else

touch $prefix.drp.temp.out

fi

if [ $? -eq 0 ]
then
	echo -e "Step 3: Combine split reads and discordant read pairs"
else
    echo "!!!Unsuccessfully extracting discordant read pairs!!!" 
	exit
fi


#Classify reads containing attB or attP in WGS (only contigs)

awk '

BEGIN {FS="\t";OFS="\t"}

function calculation1() {
     attL_start=send1-qend1+qstart2
     attL_end=send1
     attR_start=sstart2
     attR_end=sstart2+qend1-qstart2
     prophage_size=CONTIGLEN[contig]-send1+sstart2+qend1-qstart2
	 CONTIGPOSITION=contig"=::="$2
     print SR_type,contig"::"$2,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$1,CONTIGPOSITION
}

function calculation2() {
     attL_start=sstart2
     attL_end=sstart2+qend1-qstart2
     attR_start=send1
     attR_end=send1+qend1-qstart2
     prophage_size=send1+sstart2+qend1-qstart2-1
	 CONTIGPOSITION="="$2"::="contig
	 print SR_type,$2"::"contig,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$1,CONTIGPOSITION
}

function calculation3() {
     attL_start=sstart2
     attL_end=sstart2+qend1-qstart2
     attR_start=send1-qend1+qstart2
     attR_end=send1
     prophage_size=CONTIGLEN[$2]-sstart2-qend1+qstart2+send1
	 CONTIGPOSITION=$2"=::="contig
     print SR_type,$2"::"contig,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$1,CONTIGPOSITION
}

function calculation4() {
     attL_start=send1
     attL_end=send1+qend1-qstart2
     attR_start=sstart2
     attR_end=sstart2+qend1-qstart2
     prophage_size=send1+sstart2+qend1-qstart2-1
	 CONTIGPOSITION="="contig"::="$2
	 print SR_type,contig"::"$2,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$1,CONTIGPOSITION
}
	 
function calculation5() {
     attL_start=sstart2-qend1+qstart2
     attL_end=sstart2
     attR_start=send1-qend1+qstart2
     attR_end=send1
     prophage_size=send1+sstart2-qend1+qstart2-1
	 CONTIGPOSITION="="$2"::="contig
	 print SR_type,$2"::"contig,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$1,CONTIGPOSITION
}

function calculation6() {
     attL_start=send1-qend1+qstart2
     attL_end=send1
     attR_start=sstart2-qend1+qstart2
     attR_end=sstart2
     prophage_size=send1+sstart2-qend1+qstart2-1
	 CONTIGPOSITION="="contig"::="$2
	 print SR_type,contig"::"$2,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$1,CONTIGPOSITION
}

function calculation7() {
     attL_start=sstart2
     attL_end=sstart2+qend1-qstart2
     attR_start=send1
     attR_end=send1+qend1-qstart2
     prophage_size=CONTIGLEN[$2]-sstart2+CONTIGLEN[contig]-send1-qend1+qstart2+1
	 CONTIGPOSITION=$2"=::"contig"="
	 print SR_type,$2"::"contig,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$1,CONTIGPOSITION
}

function calculation8() {
     attL_start=send1
     attL_end=send1+qend1-qstart2
     attR_start=sstart2
     attR_end=sstart2+qend1-qstart2
     prophage_size=CONTIGLEN[$2]-sstart2+CONTIGLEN[contig]-send1-qend1+qstart2+1
	 CONTIGPOSITION=contig"=::"$2"="
	 print SR_type,contig"::"$2,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$1,CONTIGPOSITION
}

function calculation9() {
     attL_start=sstart2-qend1+qstart2
     attL_end=sstart2
     attR_start=send1-qend1+qstart2
     attR_end=send1
     prophage_size=CONTIGLEN[$2]-sstart2+CONTIGLEN[contig]-send1+qend1-qstart2+1
	 CONTIGPOSITION=$2"=::"contig"="
	 print SR_type,$2"::"contig,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$1,CONTIGPOSITION
}

function calculation10() {
     attL_start=send1-qend1+qstart2
     attL_end=send1
     attR_start=sstart2-qend1+qstart2
     attR_end=sstart2
     prophage_size=CONTIGLEN[$2]-sstart2+CONTIGLEN[contig]-send1+qend1-qstart2+1
	 CONTIGPOSITION=contig"=::"$2"="
	 print SR_type,contig"::"$2,evidence,attL_start,attL_end,attR_start,attR_end,prophage_size,$1,CONTIGPOSITION
}

###########contig length
FILENAME==ARGV[1] {CONTIGLEN[$1]=$2}

FILENAME==ARGV[2] {
if (FNR==1) 
{
    query = $1; contig=$2; qstart1 = $7; qend1 = $8; sstart1 = $9; send1 = $10; next
}

##make sure to extract correct flag value if "_" was used in read names
Flag_pos=split ($1,flag,"_")

if ($1 == query && $2 != contig)
{   qstart2 = $7; qend2 = $8; sstart2 = $9; send2 = $10 

#Plus-Plus or Minus-Minus	
	if (flag[Flag_pos] ~ /83|97|147|161|81|99|145|163/ && (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10))
	{	
		if ( (send2 > sstart2 && send1 > sstart1) && CONTIGLEN[contig]-send1+send2 > '"$prophage_size_min"' && CONTIGLEN[contig]-send1+send2 < '"$prophage_size_max"' )
		{
		    SR_type="SR"
            evidence="attB"
            calculation1()
		}
		else if (  (send2 > sstart2 && send1 > sstart1) && CONTIGLEN[$2]-send2+send1 > '"$prophage_size_min"' && CONTIGLEN[$2]-send2+send1 < '"$prophage_size_max"'  )
		{
		    SR_type="SR"
            evidence="attP"
			calculation3()
		}
    }		
	
#Plus-Minus or Minus-Plus	
	if (flag[Flag_pos] ~ /65|83|129|147|99|113|163|177/ && (qend2 > qend1 && qend1 > qstart2 && qstart2 > qstart1 && qstart1 < 10))
	{	
		if ( flag[Flag_pos] ~ /99|113|163|177/ && (send2 > sstart2 && send1 < sstart1) && send1+send2 > '"$prophage_size_min"' && send1+send2 < '"$prophage_size_max"' )
		{
		    SR_type="SR"
            evidence="attB"
			if($2 < contig)
			{ calculation2() }
			else if($2 > contig)
			{ calculation4() }
		}
		
		  else if ( flag[Flag_pos] ~ /99|113|163|177/ && (send2 > sstart2 && send1 < sstart1) && (CONTIGLEN[$2]+CONTIGLEN[contig]-send2-send1 > '"$prophage_size_min"') && (CONTIGLEN[$2]+CONTIGLEN[contig]-send2-send1 < '"$prophage_size_max"') )
		{
		    SR_type="SR"
            evidence="attP"
			if($2 < contig)
			{ calculation7() }
			else if($2 > contig)
			{ calculation8() }
		}
		
		else if ( flag[Flag_pos] ~ /65|83|129|147/ && (send2 < sstart2 && send1 > sstart1) && (send1+send2 > '"$prophage_size_min"') && (send1+send2 < '"$prophage_size_max"') )
		{
		    SR_type="SR"
            evidence="attP"
			if($2 < contig)
			{ calculation5() }
			else if($2 > contig)
			{ calculation6() }
		}
		
		else if ( flag[Flag_pos] ~ /65|83|129|147/ && (send2 < sstart2 && send1 > sstart1) && (CONTIGLEN[$2]+CONTIGLEN[contig]-send2-send1 > '"$prophage_size_min"') && (CONTIGLEN[$2]+CONTIGLEN[contig]-send2-send1 < '"$prophage_size_max"'))
		{
		    SR_type="SR"
            evidence="attB"
			if($2 < contig)
			{ calculation9() }
			else if($2 > contig)
			{ calculation10() }
		}
    }		


}
else if ($1 != query)
{
	query = $1; contig = $2; qstart1 = $7; qend1 = $8; sstart1 = $9; send1 = $10   
}
	
}' contiglength.file $prefix.sr.temp.2 >$prefix.WGS.sr.temp.3

awk '

#Cluster split reads
 
BEGIN{FS="\t";OFS="\t"} 

{
    if( $8 > '"$prophage_size_min"' && $8 < '"$prophage_size_max"' ) a[$9]=$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10
}

END{

count=1
T_L=int('"$read_length"'/3) #threshold low for clsutering
T_H=int('"$read_length"'*0.7) #threshold low for clsutering

for (i in a) 
{
	split(a[i],b,"\t") 
	print i,b[1],b[2],b[3],b[4],b[5],b[6],b[7],"WGS_candidate_"count

    for (j in a)
	{
		split(a[j],c,"\t")  
	    if ( i != j && c[7] == b[7])
		{
		    if (b[7]==c[7] && b[2]-c[2] > -T_L && b[2]-c[2] < T_L && b[4]-c[4] > -T_L && b[4]-c[4] < T_L && b[6]-c[6] > -T_H && b[6]-c[6] < T_H)
			{
				print j,c[1],c[2],c[3],c[4],c[5],c[6],c[7],"WGS_candidate_"count
				delete a[j]				
			}	
        }
	}
    count++		
    delete a[i]
}
}' $prefix.WGS.sr.temp.3|awk '{if (!($3 == "")) print $0}' >$prefix.WGS.sr.temp.4

awk '

BEGIN{FS="\t";OFS="\t"} 
{
if ( !(a[$NF,1] >0) )
{
    a[$NF,1]=0
}

if ( !(a[$NF,2] >0) )
{
    a[$NF,2]=0
}

if ($2 == "attB")
{
    a[$NF,1]+=1 
	sum[$NF,3]+=$3
	sum[$NF,4]+=$4
	sum[$NF,5]+=$5
	sum[$NF,6]+=$6
	sum[$NF,7]+=$7
	a[$NF,8]=$8
}
else if ($2 == "attP")
{
    a[$NF,2]+=1
	sum[$NF,3]+=$3
	sum[$NF,4]+=$4
	sum[$NF,5]+=$5
	sum[$NF,6]+=$6
	sum[$NF,7]+=$7
	a[$NF,8]=$8
}
}
END{
    for(i in a) 
	{
	    split(i,idx,SUBSEP)
        a[idx[1],3]=int(sum[idx[1],3]/(a[idx[1],1]+a[idx[1],2]))
	    a[idx[1],4]=int(sum[idx[1],4]/(a[idx[1],1]+a[idx[1],2]))
	    a[idx[1],5]=int(sum[idx[1],5]/(a[idx[1],1]+a[idx[1],2]))
	    a[idx[1],6]=int(sum[idx[1],6]/(a[idx[1],1]+a[idx[1],2]))
	    a[idx[1],7]=int(sum[idx[1],7]/(a[idx[1],1]+a[idx[1],2]))
		print idx[1],a[idx[1],1],a[idx[1],2],a[idx[1],1]+a[idx[1],2],a[idx[1],3],a[idx[1],4],a[idx[1],5],a[idx[1],6],a[idx[1],7],a[idx[1],8]
	}
}

'  $prefix.WGS.sr.temp.4 | awk '!a[$0]++' > $prefix.WGS.sr.temp.out

#Step 2: Extracting discordant read pairs

awk --re-interval '$1 ~ /^@/ || ($2 ~ /97|145|81|161|177|113|129|65/ && ($6 ~ /^[1][0-5][0-9]M/ || $6 ~ /^[6-9][0-9]M/) )' $sam_file | awk '!a[$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8]++' | sort -k1,1 -k2n,2 | awk '

BEGIN {FS="\t";OFS="\t"}

{
if (NR==1) 
{
    read = $1; read_flag = $2; len=$9; a[$1] = $0; next
}


if ($1 == read)
{
    if (read_flag ~ /97|145|81|161/)
	{
	    print "DRP","attPorattB",a[$1]"\n""DRP","attPorattB",$0
	}
	
	else if (read_flag ~ /177|113|129|65/)
	{
	    print "DRP","attPorattB",a[$1]"\n""DRP","attPorattB",$0
	}
}	

else if ($1 != read)
{
	read = $1; read_flag = $2; len=$9; a[$1] = $0; next   
}
	
}' | awk 'BEGIN {FS="\t";OFS="\t"} ($4 ~ /97|81|65|113/ )' | awk 'BEGIN {FS="\t";OFS="\t"} NR==FNR {split($9,readname,"_"); a[readname[1]]=readname[1]} NR>FNR { if(!($3 in a)) print $0}' $prefix.WGS.sr.temp.3 - >$prefix.WGS.drp.temp.1


awk ' BEGIN {FS="\t";OFS="\t"}

FILENAME==ARGV[1] { a[$3]=$0; next} 
FILENAME==ARGV[2] {  

    T_L=int('"$insert_size"'*1.6) #threshold low for clsutering
    T_H=int('"$insert_size"'*3.5) #threshold low for clsutering
 
	attP_count = 0
	attB_count = 0
	for (i in a)
	{
	   split(a[i],b,"\t")
	   split($10,c,"[::]"); split(c[1],d,"=");split(c[3],e,"=")
	   if(d[2]=="" && e[1]=="")
	   {
	        if (b[5]==d[1] && b[9]==e[2])
	        {
	            if($6-b[6] > -T_L && $6-b[6] <T_L && $7-b[10] >-T_L && $7-b[10] <T_L )
	            {   
	                if ($6-b[6]+10 > '"$read_length"' && b[10] >= $7)
		            {
		                attB_count++
		                delete a[i]
		            }
	                else if (b[6] >= $5 && $8-b[10]+10 > '"$read_length"'  )
		            {
		                attP_count++ 
		                delete a[i]
		            }             			  
		        }
		    }
		    else if (b[5]==e[2] && b[9]==d[1])
		    {
		        if($6-b[10] > -T_L && $6-b[10] <T_L && $7-b[6] >-T_L && $7-b[6] <T_L )
	            {   
	                if ($6-b[10]+10 > '"$read_length"' && b[6] >= $7)
		            {
		                attB_count++
		                delete a[i]
		            }
	                else if ( b[10] >= $5 && $8-b[6]+10 > '"$read_length"'  )
		            {
		               attP_count++  
		               delete a[i]
                    }				   
		        }		
		    }
		}
		else if(d[1]=="" && e[1]=="")
		{
	        if (b[5]==d[2] && b[9]==e[2])
	        {
	            if($6-b[6] > -T_L && $6-b[6] <T_L && $7-b[10] >-T_L && $7-b[10] <T_L )
	            { 
	                if (b[6] >= $5 && b[10] >= $7)
		            {
		                attB_count++
		                delete a[i]
		            }
	                else if ($6-b[6]+10 >= '"$read_length"' && $8-b[10]+10 >= '"$read_length"')
		            {
		                attP_count++ 
		                delete a[i]
		            }
		        }
		    }
	        else if (b[5]==e[2] && b[9]==d[2])
	        {
	            if($7-b[6] > -T_L && $7-b[6] <T_L && $6-b[10] >-T_L && $6-b[10] <T_L )
	            {   
	                if ((b[10] >= $5 && b[6] >= $7))
		            {
		                attB_count++
		                delete a[i]
		            }
	                else if ( $8-b[6]+10 >= '"$read_length"' && $6-b[10]+10 > '"$read_length"')
		            {
		                attP_count++ 
		                delete a[i]
		            }
		        }
		    }
		}		
		else if(d[2]=="" && e[2]=="")
		{
	        if (b[5]==d[1] && b[9]==e[1])
	        {
	            if($6-b[6] > -T_L && $6-b[6] <T_L && $7-b[10] >-T_L && $7-b[10] <T_L )
	            {   
	                if ($6-b[6]+10 >= '"$read_length"' && $8-b[10]+10 > '"$read_length"')
		            {
		                attB_count++
		                delete a[i]
		            }
	                else if (b[6] >= $5 && b[10] >= $7)
		            {
		                attP_count++ 
		                delete a[i]
		            }
		        }
		    }
	        else if (b[5]==e[1] && b[9]==d[1])
	        {
	            if($7-b[6] > -T_L && $7-b[6] <T_L && $6-b[10] >-T_L && $6-b[10] <T_L )
	            {   
	                if ($8-b[6]+10 >= '"$read_length"' && $6-b[10]+10 > '"$read_length"')
		            {
		                attB_count++
		                delete a[i]
		            }
	                else if (b[6] >= $7 && b[10] >= $5)
		            {
		                attP_count++ 
		                delete a[i]
		            }
		        }
		    }
		}		
	}
    f[$1]=$0"\t"attB_count"\t"attP_count"\t"attB_count+attP_count
	}
END	{
	for (i in f)
    {
	    print f[i]
	}
	for (i in a)
    {
	    print a[i] >> "'$prefix'.WGS.drp.temp.left"
	}
}' $prefix.WGS.drp.temp.1 $prefix.WGS.sr.temp.out >$prefix.WGS.sr-drp.temp.out




# Step 3: Combining split reads and discordant read pairs

##Check the input file
if test -s $prefix.drp.temp.out; then

cat $prefix.sr-drp.temp.out $prefix.WGS.sr-drp.temp.out | awk '

BEGIN {FS="\t";OFS="\t"}

FILENAME==ARGV[1] { print $0; next} 
FILENAME==ARGV[2]  { print $1,"0","0","0",$5,"-",$6,"-",$7,$8,$2,$3,$4}

' - $prefix.drp.temp.out >$prefix.temp.out

else

cat $prefix.sr-drp.temp.out $prefix.WGS.sr-drp.temp.out >$prefix.temp.out

fi




#Delet existed SR_evidence.list

if [ -f "$prefix.SR_evidence.list1" ]
then
 
 rm $prefix.SR_evidence.list1
 
fi

awk ' BEGIN {FS="\t";OFS="\t";i=1;print "prophage_candidate","contig","attL_start","attL_end","attR_start","attR_end","prophage_size","SR_evidence_attB","SR_evidence_attP","DRP_evidence_attB","DRP_evidence_attP"}


{		
	if ( $9 > '"$prophage_size_min"' && $9 < '"$prophage_size_max"')  
	{
	    if(($4 >= '"$SR_COUNT"' && $13 >='"$DRP_COUNT"'))
	    {    
		    if ($4 != 0 && $5 > 1 && $7 >1)
		    {
		        if ( $6-$5 >= '"$att_length_min"' )
			    {    
				    print "candidate_"i,$10,$5,$6,$7,$8,$9,$2,$3,$11,$12
					print $1 >> "'$prefix'.SR_evidence.list1"
	                i++
				}
			}
			else if ($5 > 1 && $7 >1)
			{
				print "candidate_"i,$10,$5,$6,$7,$8,$9,$2,$3,$11,$12
				i++
			}
			
        }				    
	}
	
}' $prefix.temp.out > $prefix.prophage.out



#blastn using split reads

awk 'BEGIN {FS="\t";OFS="\t"} NR==FNR {a[$1]=$0} NR>FNR {if ($9 in a) print $1} ' $prefix.SR_evidence.list1 $prefix.sr.temp.3 > $prefix.SR_evidence.list

awk '!a[$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8]++' $prefix.sr.temp.1 | awk '!a[$1]++'| awk 'BEGIN {FS="\t";OFS="\t"} NR==FNR {a[$1]=$0} NR>FNR {if ($1"_"$2 in a) print $0} ' $prefix.SR_evidence.list - | awk '{printf ">"$1"_"$2"\n"$10"\n"}' > $prefix.SR.reads.fasta

blastn -query $prefix.SR.reads.fasta -db $prefix.nuclDB -out $prefix.SR.blastn.result -evalue 1e-3 -outfmt 1 -word_size 11


#Clean up
if [ -f "$prefix.SR_evidence.list1" ]
then 
rm $prefix.SR_evidence.list1
fi

rm contiglength.file $prefix.sr.temp.1 $prefix.reads.fasta blastn.log makeblastdb.log $prefix.sr.temp.2 $prefix.sr.temp.3 $prefix.sr.temp.out $prefix.drp.temp.1 $prefix.drp.temp.left $prefix.sr-drp.temp.out $prefix.drp.temp.2 $prefix.drp.temp.out $prefix.WGS.sr.temp.3 $prefix.WGS.sr.temp.4 $prefix.WGS.sr.temp.out $prefix.WGS.sr-drp.temp.out $prefix.temp.out $prefix.nuclDB.*


if [ $? -eq 0 ]
then
    echo "Successfully running prophage_GPS.sh"
	echo -e "checking predicted prophage in $prefix.prophage.out\n\n"
else
    echo "!!!Unsuccessfully running prophage_GPS.sh!!!" 
	exit
fi

exit 0
 
