#!/bin/bash
echo "Generating DNA" | awk '

BEGIN {FS="\t";OFS="\t";srand();print "genome_name","Att_site_length","Prophage_lenght","attL_start","attL_end","attR_start","attR_end","Genome_size","GC_content" >> "genome_position.tsv"} 

##############################Function_start
function randint(n) {
     return int(n *rand()) + 1
}

function randrang(n,m) {
     return int((m-n) *rand()) + n
}

function mutation1(Att1) {
     Att1=substr(Att1,1,pos1-1)substr("ATCG", randint(4), 1)substr(Att1,pos1+1)
	 return Att1
}
function mutation2(Att2) {
     if (pos1 < pos2)
     {Att2=substr(Att2,1,pos1-1)substr("ATCG", randint(4), 1)substr(Att2,pos1+1,pos2-pos1-1)substr("ATCG", randint(4), 1)substr(Att2,pos2+1)}
     else if (pos1 > pos2)
     {Att2=substr(Att2,1,pos2-1)substr("ATCG", randint(4), 1)substr(Att2,pos2+1,pos1-pos2-1)substr("ATCG", randint(4), 1)substr(Att2,pos1+1)}
	 return Att2
}

function mutation(Att){
    mismatch=int(rand()*100%2+1)
	print "mismatch",mismatch

    if ( mismatch == 1 )
	{    
        len=length(Att)
        pos1=randint(len)
        do
		{
		    pos1=randint(len)
		}while(pos1 == 1 || pos1 == len)  
        return mutation1(Att)
    } 
	else if ( mismatch == 2 )
	{    
        len=length(Att)
        pos1=randint(len)
        pos2=randint(len)
		do
        {   pos1=randint(len)
            pos2=randint(len)
        }while(pos1 == 1 || pos1 == len || pos2 == 1 || pos2 == len || pos1 == pos2)

        return mutation2(Att)
    }
}


##############################Function_end


{

	
for (count = 1; count <= 20; count++)
    {
        do
		{
		    GC_CONTENT=int(rand()*100)/100
	    }while(GC_CONTENT <= 0.20 ||GC_CONTENT >= 0.80)

		GENOME_SIZE=4000000
		    
		
		Genome_left_len=randrang(1000000,3000000)
	    Genome_rigth_len=GENOME_SIZE-Genome_left_len
	    Att_site_len=randrang(2,145)
	    Prophage_len=randrang(5000,150000)

	    "python3 -c \047import random_DNA; print(random_DNA2.biased_DNA_generator("GC_CONTENT","Genome_left_len"))\047"|getline Genome_left
	    "python3 -c \047import random_DNA; print(random_DNA2.biased_DNA_generator("GC_CONTENT","Genome_rigth_len"))\047"|getline Genome_rigth
	    "python3 -c \047import random_DNA; print(random_DNA2.biased_DNA_generator("GC_CONTENT","Att_site_len"))\047"|getline Att_site_L
        
		if (length(Att_site_L) <=4)
		{
		    Att_site_R=Att_site_L
		}
		else
		{
		    Att_site_R=mutation(Att_site_L)
		}	

	    "python3 -c \047import random_DNA; print(random_DNA2.biased_DNA_generator("GC_CONTENT","Prophage_len"))\047"|getline Prophage
		
		Prophage_c=substr(Prophage,int(Prophage/2))substr(Att_site_R,1,int(Att_site_len/2))substr(Att_site_L,int(Att_site_len/2)+1)substr(Prophage,1,int(Prophage/2)-1)
		
		AttB=substr(Att_site_L,1,int(Att_site_len/2))substr(Att_site_R,int(Att_site_len/2)+1)
		AttP=substr(Att_site_R,1,int(Att_site_len/2))substr(Att_site_L,int(Att_site_len/2)+1)
	    
		print ">sg_"count >> "genome_l.fasta"
		print Genome_left >> "genome_l.fasta"
		
		print ">sg_"count >> "attL.fasta"
		print Att_site_L >> "attL.fasta"
		
		print ">sg_"count >> "prophage.fasta"
		print Prophage >> "prophage.fasta"

		print ">sg_"count >> "attR.fasta"
		print Att_site_R >> "attR.fasta"
		
		print ">sg_"count >> "genome_r.fasta"
		print Genome_rigth >> "genome_r.fasta"

		print ">sg_"count >> "attB.fasta"
		print AttB >> "attB.fasta"

		print ">sg_"count >> "attP.fasta"
		print AttP >> "attP.fasta"
		
		print ">sg_"count >> "prophage_c.fasta"
		print Prophage_c >> "prophage_c.fasta"

#################################################################
#reporting positions
#################################################################
		print "sg_"count,Att_site_len,Prophage_len,Genome_left_len+1,Genome_left_len+Att_site_len,Genome_left_len+Att_site_len+Prophage_len+1,Genome_left_len+Att_site_len+Prophage_len+Att_site_len, Genome_left_len+Att_site_len+Prophage_len+Att_site_len+Genome_rigth_len,GC_CONTENT >> "genome_position.tsv"
	}
}'

seqkit concat <(cat genome_l.fasta) <(cat attL.fasta) <(cat prophage.fasta) <(cat attR.fasta) <(cat genome_r.fasta) |seqkit replace -p "sg_" -r 'genome_sg_' > genomes.fasta
seqkit split genomes.fasta -i -O genomes

seqkit concat <(cat genome_l.fasta) <(cat attB.fasta) <(cat genome_r.fasta) |seqkit replace -p "sg_" -r 'attB_sg_' > genomes_attB.fasta
seqkit split genomes_attB.fasta -i -O attB

cat prophage_c.fasta |seqkit replace -p "sg_" -r 'attP_sg_' > prophage_circle.fasta
seqkit split prophage_circle.fasta -i -O attP



### Generating abundance file for metegenomic mode of GemSIM
for ((j=1; j<= 20; j++))
do
    mkdir sg_$j
    cp ./genomes/genomes.id_genome_sg_${j}.fasta ./attB/genomes_attB.id_attB_sg_${j}.fasta ./attP/prophage_circle.id_attP_sg_${j}.fasta ./sg_${j}/
	ls ./sg_${j}/| awk '{FS="\t";OFS="\t"} NR==1{print $0"\t"1}NR==2{print $0"\t"998}NR==3{print $0"\t"1}' > sg_${j}_abundance1.txt
	ls ./sg_${j}/| awk '{FS="\t";OFS="\t"} NR==1{print $0"\t"1}NR==2{print $0"\t"989}NR==3{print $0"\t"10}' > sg_${j}_abundance2.txt
	ls ./sg_${j}/| awk '{FS="\t";OFS="\t"} NR==1{print $0"\t"10}NR==2{print $0"\t"980}NR==3{print $0"\t"10}' > sg_${j}_abundance3.txt
	ls ./sg_${j}/| awk '{FS="\t";OFS="\t"} NR==1{print $0"\t"10}NR==2{print $0"\t"990}NR==3{print $0"\t"990}' > sg_${j}_abundance4.txt

done

