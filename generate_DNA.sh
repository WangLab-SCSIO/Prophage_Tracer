#!/bin/bash
echo "Generating DNA" | awk '

BEGIN {FS="\t";OFS="\t";srand()} 

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

function mutation(Att)
{
    mismatch=int(rand()*100%2+1)

    if ( mismatch == 1 )
	{    
        len=length(Att)
        pos1=randint(len)
        while(pos1 == 1 || pos1 == len)     
        {pos1=randint(len)}
        return mutation1(Att)
    } 
	else if ( mismatch == 2 )
	{    
        len=length(Att)
        pos1=randint(len)
        pos2=randint(len)
        while(pos1 == 1 || pos1 == len || pos2 == 1 || pos2 == len || pos1 == pos2)
        {   pos1=randint(len)
            pos2=randint(len)
        }
        return mutation2(Att)
    }
}


##############################Function_end


function generate_dna(LEN) {
     dna=""
	 do
	 {	 
	     for( i = 1; i <= LEN; i++) 
         {
            dna=dna substr("ATCG", randint(4), 1)
         }
		 dna_check=dna
		 gsub(/[^GC]/,"",dna)
	 }
     while(length(dna)/length(dna_check) < 0.4 || length(dna)/length(dna_check) > 0.7 )
	 return dna_check
}


{

	
	for (count = 1; count <= 20; count++)
	{    
	    Genome_left_len=randrang(5000,3995000)
	    Genome_rigth_len=4000000-Genome_left_len
	    Att_site_len=randrang(10,100)
	    Prophage_len=randrang(5000,150000)

	    Genome_left=generate_dna(Genome_left_len)
	    Genome_rigth=generate_dna(Genome_rigth_len)
	    Att_site_L=generate_dna(Att_site_len)
		Att_site_R=mutation(Att_site_L)
	    Prophage=generate_dna(Prophage_len)
		
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
		print "sg_"count,Att_site_len,Prophage_len,Genome_left_len+1,Genome_left_len+Att_site_len,Genome_left_len+Att_site_len+Prophage_len+1,Genome_left_len+Att_site_len+Prophage_len+Att_site_len >> "genome_position.tsv"

	}
}'

seqkit concat <(cat genome_l.fasta) <(cat attL.fasta) <(cat prophage.fasta) <(cat attR.fasta) <(cat genome_r.fasta) |seqkit replace -p "sg_" -r 'genome_sg_' > genomes.fasta
seqkit split genomes.fasta -i -O genomes

seqkit concat <(cat genome_l.fasta) <(cat attB.fasta) <(cat genome_r.fasta) |seqkit replace -p "sg_" -r 'attB_sg_' > genomes_attB.fasta
seqkit split genomes_attB.fasta -i -O attB

cat prophage_c.fasta |seqkit replace -p "sg_" -r 'attP_sg_' > prophage_circle.fasta
seqkit split prophage_circle.fasta -i -O attP

for ((j=1; j<= 20; j++))
do
    mkdir sg_$j
    cp ./genomes/genomes.id_genome_sg_${j}.fasta ./attB/genomes_attB.id_attB_sg_${j}.fasta ./attP/prophage_circle.id_attP_sg_${j}.fasta ./sg_${j}/
	ls ./sg_${j}/| awk '{FS="\t";OFS="\t"} NR==1{print $0"\t"998}NR==2{print $0"\t"1}NR==3{print $0"\t"1}' > sg_${j}_abundance1.txt
	ls ./sg_${j}/| awk '{FS="\t";OFS="\t"} NR==1{print $0"\t"989}NR==2{print $0"\t"1}NR==3{print $0"\t"10}' > sg_${j}_abundance2.txt
	ls ./sg_${j}/| awk '{FS="\t";OFS="\t"} NR==1{print $0"\t"980}NR==2{print $0"\t"10}NR==3{print $0"\t"10}' > sg_${j}_abundance3.txt
	ls ./sg_${j}/| awk '{FS="\t";OFS="\t"} NR==1{print $0"\t"990}NR==2{print $0"\t"10}NR==3{print $0"\t"990}' > sg_${j}_abundance4.txt

done

echo "Found each group of genomes in each folder containing original bacterial genome, bacterial genome with prophage excised and circular prophage genome"

