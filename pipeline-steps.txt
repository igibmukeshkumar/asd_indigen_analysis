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

# b. ClinVar 
awk -F'\t' 'NR == 1 || tolower($170) ~ /pathogenic/' asd_gene_panel.hg38_multianno_0.05.txt | awk -F'\t' '{print $0 "\t" "ClinVar" }'> clinvar.txt

# c. SIFT_score, PolyPhen_HAVR_score, CADD_phred (17: D<=0.05, 26: P/D>=0.447, 108: D>=20, )
awk -F'\t' 'NR==1 { print } NR>1 && ($17 == "." || $26 == ".") { $17 = 0; $26 = 0; } ($17 <= 0.05 && $26 >= 0.447) { print  }'   asd_gene_panel.hg38_multianno_0.05.txt > 1.txt
awk -F'\t' 'NR==1 { print } NR>1 && ($17 == "." || $108 == ".") { $17 = 0; $108 = 0; } ($17 <= 0.05 && $108 >= 20) { print }'   asd_gene_panel.hg38_multianno_0.05.txt > 2.txt
awk -F'\t' 'NR==1 { print } NR>1 && ($108 == "." || $26 == ".") { $108 = 0; $26 = 0; } ($108 >= 20 && $26 >= 0.447) { print }'   asd_gene_panel.hg38_multianno_0.05.txt > 3.txt
cat 1.txt 2.txt 3.txt | sort -k2n -u | awk -F'\t' '{print $0 "\t" "PP3" }' >  pp3.txt

# d. LoF variants
egrep -w 'Chr|frameshift deletion|frameshift insertion|stopgain|startloss|splicing' asd_gene_panel.hg38_multianno_0.05.txt | awk -F'\t' '{print $0 "\t" "LOF_variants" }' >  lof.txt

# join carefully: after catenate clinvar (benign, conflicting; any match in CLNSIG removed i.e stringent filter)
cat lof.txt clinvar.txt pp3.txt | awk -F'\t' 'NR==1 {for (i=1; i<NF; i++) printf "%s\t", $i; print "Prediction"; next} tolower($170) !~ /benign|conflicting/ {key=""; for(i=1;i<NF;i++) key=key $i OFS; values[key]=(key in values) ? values[key]","$NF:$NF} END {for(k in values) print k values[k]}' OFS='\t' > PCL.txt

# Check in major public database (gnomAD exome, genome, 1000GP_Aug2016, ESP6500, ExAC)
awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "Databases_Status"; next} {status = ($16=="." && $211=="." && $193=="." && $183=="." && $192=="." && $184=="." && $170==".") ? "Novel" : "."; print $0, status}' PCL.txt > asd_panel_final.txt

# 5. For GT data from IndiGen data: Normalize (variants=)
bcftools view -R ASD_STDIN.bed  /home/binukumar/storage300tb/IndiGen_data/indigen_chr_wise_vcf/new_vcfs/IndiGen_all.vcf.gz | vt decompose -s - | vt normalize -n  -r /lustre/binukumar/tools/ref_files/dragmap_ref/hg38_dragen.fa - | bgzip -c > /home/binukumar/storage300tb/Project_ASD/asd_GT/asd_GT_chr_wise_VT.vcf.gz

# AC HOMO & HET: VT NORM Form Chr_wise combined VCF  
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' asd_GT_chr_wise_VT.vcf.gz | \
awk 'BEGIN {
  print "CHROM\tPOS\tREF\tALT\tHET_AC\tHOM_AC\tAC_total"
}
{
  het_ac = 0;
  hom_ac = 0;
  for (i = 5; i <= NF; i++) {
    if ($i == "1/1" || $i == "1|1") {
      hom_ac++;
    } else if ($i == "0/1" || $i == "1/0" ||
               $i == "0|1" || $i == "1|0" ||
               $i == "./1" || $i == "1/." ||
               $i == "1|." || $i == ".|1") {
      het_ac++;
    }
  }
  hom_ac_alleles = hom_ac * 2;
  ac_total = het_ac + hom_ac_alleles;
  print $1 "\t" $2 "\t" $3 "\t" $4 "\t" het_ac "\t" hom_ac_alleles "\t" ac_total
}' > chr_wise_asd_AC_VT.txt


# Calculate Prevalence and Carrier in IndGen For ASD
#code prepared carrier_prevalence_count.sh 
# lof_mis_lp_p.txt=CHROM POS REF ALT 

awk 'FNR==NR {k[$1 FS $2 FS $3 FS $4]; next} 
     /^#/ {print; next} 
     (($1 FS $2 FS $4 FS $5) in k)' lof_mis_lp_p.txt <(zcat asd_GT_chr_wise_VT.vcf.gz) > lp_p_cr_pr.vcf

{ grep '^#CHROM' lp_p_cr_pr.vcf | sed 's/^#//'; bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' lp_p_cr_pr.vcf; } > pr_cr_1.txt


CHROM   POS     REF     ALT
chr12   2593255 G       A


