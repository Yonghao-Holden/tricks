#!/usr/bin/awk -f

## this script convert gffread generated gff3 file into a modified gff3 file so that misopy sashimi package can index it
## Usage: sed 's/transcript/gene/g' input.gff3 | sed 's/gene_name/Name/g' | sed 's/geneID/Alias/g' | prep_gff3_for_aa_sequence_extraction.awk > output.sashimi.gff3
{
	if ($3 == "gene"){
		exon_counter=0
		split($9,old_attributes,";");
		#print(old_attributes[1],old_attributes[2],old_attributes[3])
		for (i=1;i<=3;i++){
			split(old_attributes[i],temp,"=");
			#print(temp[1], temp[2])
			if (temp[1] == "ID"){
				temp_ID=temp[2];
			}
			if (temp[1] == "Alias"){
				temp_Alias=temp[2];
			}
			if (temp[1] == "Name"){
				temp_Name=temp[2];
			}
		}
		gene_attributes=sprintf("ID=%s_%s;Name=%s;Alias=%s",temp_Name,temp_ID,temp_Name,temp_Alias);
		mRNA_attributes=sprintf("ID=%s;Parent=%s_%s",temp_ID,temp_Name,temp_ID)
		printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,"gene",$4,$5,$6,$7,$8,gene_attributes);
		printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,"mRNA",$4,$5,$6,$7,$8,mRNA_attributes);
	}
	else if($3 == "exon"){
		exon_counter += 1;
		split($9,parent_exon,"=");
		exon_attributes=sprintf("ID=%s:exon_%s;%s",parent_exon[2],exon_counter,$9)
		printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,exon_attributes);
	}
}
