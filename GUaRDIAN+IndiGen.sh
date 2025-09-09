#!/bin/bash

#-----------------------------------------
# Full pipeline: Mask low-depth genotypes, compute AF/AC/AN, het/hom, annotate VCF
# Requires: GATK 4.3.0.0, bcftools, awk
#-----------------------------------------

# Input VCF (genomeDB-imported)
VCF_IN="case_v2_1_GTGVCFgenDB.vcf.gz"

# Temporary and output files
VCF_MASKED_TMP="case.masked.tmp.vcf.gz"
TABLE_TMP="case.masked.table.tsv"
FINAL_TABLE="per_site_AF_het_hom.tsv"
VCF_FINAL="case.final.AF.vcf.gz"

# Reference genome
REF="/lustre/binukumar/tools/ref_files/dragmap_ref/hg38_dragen.fa"

#-----------------------------------------
echo "Step 1: Mask genotypes with DP < 10"
gatk VariantFiltration \
   -V "$VCF_IN" \
   --genotype-filter-expression "DP < 10" \
   --genotype-filter-name "LowDP" \
   --set-filtered-genotype-to-no-call true \
   -O "$VCF_MASKED_TMP"

echo "Step 2: Convert masked VCF to table with genotypes"
gatk VariantsToTable \
   -V "$VCF_MASKED_TMP" \
   -F CHROM -F POS -F REF -F ALT -F GT \
   -O "$TABLE_TMP"

echo "Step 3: Compute AC, AN, AF, het, hom counts"
awk '{
    AC=0; AN=0; het=0; hom=0;
    for(i=5;i<=NF;i++){
        if($i=="0/1" || $i=="1/0"){ AC+=1; AN+=2; het++ }
        else if($i=="1/1"){ AC+=2; AN+=2; hom++ }
        else if($i=="0/0"){ AN+=2 }
    }
    AF=(AN>0)?AC/AN:0
    print $1,$2,$3,$4,AC,AN,AF,het,hom
}' OFS="\t" "$TABLE_TMP" > "$FINAL_TABLE"

echo "Step 4: Annotate original VCF with new AC, AN, AF, het, hom"
bcftools annotate \
   -a "$FINAL_TABLE" \
   -c CHROM,POS,REF,ALT,AC,AN,AF,het,hom \
   -O z -o "$VCF_FINAL" \
   "$VCF_MASKED_TMP"

echo "Step 5: Index final VCF"
bcftools index "$VCF_FINAL"

echo "Pipeline complete. Final VCF: $VCF_FINAL"
