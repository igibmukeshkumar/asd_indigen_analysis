#!/bin/bash


#SBATCH -n 10
#SBATCH -c 1
#SBATCH -p compute
#SBATCH -J ALN_DRAGEN
#SBATCH --output=DR_%j_out.log
#SBATCH --error=DR_%j_err.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sainimuk77@gmail.com

# 1. Download refGene GTF spanning all transcripts
Download Link for MANE Transcripts/Genes: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/
zgrep -v '^#' MANE.GRCh38.v1.4.refseq_genomic.gtf.gz | sed -e 's/"//g' | tr ';' '\t' | awk 'BEGIN{FS=OFS="\t"} $3=="gene" {print $1, $4, $5, $3, $9}' | sed -E 's/gene_id //g' > MANE_gene_bed.txt
awk 'NR==FNR{genes[$1]; next} $NF in genes'  asd_gene_panel.txt MANE_gene_bed.txt | bedtools sort -i - | bedtools merge -i - >  MANE_ASD_GP_STDIN.bed


# 2. Capturing ASD genes based variants from IndiGen data VCF
bcftools view -R pr_gene.bed -Oz -o asd_gene_panel.vcf.gz IndiGen_AF_norm_final.vcf.gz
bcftools view -R pr_gene.bed -Oz -o asd_gene_panel.vcf.gz /home/binukumar/storage300tb/IndiGen_data/AF_IndiGen/annovar_indigen/IndiGen_AF_norm_filal.vcf.gz


# 3. Annotation of variants (ANNOVAR)
perl /lustre/binukumar/tools/annovar_2025/annovar/table_annovar.pl --thread 20 asd_gene_panel.vcf.gz /lustre/binukumar/tools/annovar_2025/annovar/humandb -buildver hg38 \
-out asd_gene_panel \
-remove \
-protocol refGeneWithVer,refGene,avsnp151,dbnsfp47a,clinvar_20250409,indigen1029,1000g2015aug_all,exac03,esp6500siv2_all,gnomad41_exome,gnomad41_genome,gme,dbscsnv11,regsnpintron,revel \
-operation g,g,f,f,f,f,f,f,f,f,f,f,f,f,f \
-otherinfo -nastring . \
-polish -otherinfo -vcfinput

# 4. Variant filtering 
# a. Allele Frequency (<5%) - 1000GP, gnomad41_exome_AF, gnomad41_genome_AF
awk -F'\t' '!($9 == "synonymous SNV" || $9 == "unknown")' asd_gene_panel.hg38_multianno.txt | awk -F '\t' '(NR==1) || ($183 < 0.05) && ($193 < 0.05) && ($211 < 0.05)' > asd_gene_panel.hg38_multianno_0.05.txt

# b. REVEL >0.75; deleterious variantsawk -F'\t' '
awk -F'\t' 'NR==1 || ($211 > 0.75)' asd_gene_panel.hg38_multianno_0.05.txt > asd_gene_panel.hg38_multianno_0.05_revel0.75.txt

# 5 LOFTEE: Annotation of variants with VEP and Selecting High Confidence variants
