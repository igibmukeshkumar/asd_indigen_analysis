# asd_indigen_analysis
# Idea Behind Carrier and Prevalence Rate 
Assumptions:
1. variants of a gene invloved in one codition like WD
2. multple genes may invloved with one condition  like ASD
3. one gene may invloved with multiple conditons (pleiotropy) like FMR1 gene (Fragile X syndrome; FXS, Fragile X-associated tremor/ataxia syndrome; FXTAS, and Fragile X-associated primary ovarian insufficiency; FXPOI)
4. 


Step 1: Calculation of Prevalence and Carrier Rate based on the genotype data
--> write python code (script) for the above 

use para for hep and supply from outsie
Usage: ./cr_pr.sh -gf <GenotypeFile> -qf <QueryFile> -o <OutputFile>
  -gf, --GTfile       Path to genotype matrix file
  -qf, --Queryfile    Path to query file (variant-gene-moi-condition)
  -o,  --Output       Output file path
  -h,  --help         Display this help message
  
Step 2: Data to be Processed (input files)
example: genotype file 
variant	s1	s2	s3	s4
v1	0/1	0/0	0/0	0/0
v2	1/1	0/0	./1	0/0
v3	1/1	0/0	0/0	0/0
v4	0/1	0/1	0/0	0/0
v5	0/0	0/1	0/0	1|1
v6	0/1	0/0	0/0	1/0
v7	0/0	0/0	0/0	0/1
v8	1|1	0/0	1/0	0/0
v9	0/1	0/0	0/0	0/0
v10	1/1	0/0	./1	0/0
v11	1/1	0/0	0/0	0/0
v12	0/1	0/1	0/0	0/0
v13	0/0	0/1	0/0	1|1
v14	0/1	0/0	0/0	1/0
v15	0/0	0/0	0/0	0/1
v16	1|1	0/0	1/0	0/0

example: query file
variant	gene	moi	condition
v1	g1	AR	ASD
v2	g1	AR	ASD
v3	g2	AR	ASD
v4	g2	AR	ASD
v5	g3	AR	ASD
v6	g3	AR	ASD
v7	g4	AR	MD
v8	g5	AR	MD
v9	g6	AD	ASD
v10	g6	AD	ASD
v11	g7	AD	ASD
v12	g7	AD	ASD
v13	g8	AD	ASD
v14	g8	AD	ASD
v15	g9	AD	MD
v16	g10	AD	MD

Find matches for variant based on the query file in genotype file and incoporate the gene, moi and condition from query file for the further operations

Step 4:
Count the below het/hom based on the below parameters are explained as follows;
--->
het_raw_per_variant = count het (./1,.|1,1/.,1|.,0/1,0|1) per variant 
hom_raw_per_variant = count hom (1|1, 1/1) per variant
1. count row wise as variant depict row and put in the respective cols as explained
eg
variant	gene	moi	condition	s1	s2	s3	s4	het_raw_per_variant	hom_raw_per_variant
v1	g1	AR	ASD	0/1	0/0	0/0	0/0	1	0
v2	g1	AR	ASD	1/1	0/0	./1	0/0	1	1
v3	g2	AR	ASD	1/1	0/0	0/0	0/0	0	1
v4	g2	AR	ASD	0/1	0/1	0/0	0/0	2	0
v5	g3	AR	ASD	0/0	0/1	0/0	1|1	1	1
v6	g3	AR	ASD	0/1	0/0	0/0	1/0	2	0
v7	g4	AR	MD	0/0	0/0	0/0	0/1	1	0
v8	g5	AR	MD	1|1	0/0	1/0	0/0	1	1
v9	g6	AD	ASD	0/1	0/0	0/0	0/0	1	0
v10	g6	AD	ASD	1/1	0/0	./1	0/0	1	1
v11	g7	AD	ASD	1/1	0/0	0/0	0/0	0	1
v12	g7	AD	ASD	0/1	0/1	0/0	0/0	2	0
v13	g8	AD	ASD	0/0	0/1	0/0	1|1	1	1
v14	g8	AD	ASD	0/1	0/0	0/0	1/0	2	0
v15	g9	AD	MD	0/0	0/0	0/0	0/1	1	0
v16	g10	AD	MD	1|1	0/0	1/0	0/0	1	1
--->
het_raw_per_gene = count het (./1,.|1,1/.,1|.,0/1,0|1) per gene 
hom_raw_per_gene = count hom (1|1, 1/1) per gene 
1. look for variants/variant present in the gene (gene based) and count genotypes in each sample and put the sum of counts in respective cols
eg
variant	gene	moi	condition	s1	s2	s3	s4	het_raw_per_gene	hom_raw_per_gene
v1	g1	AR	ASD	0/1	0/0	0/0	0/0	2	1
v2	g1	AR	ASD	1/1	0/0	./1	0/0	2	1
v3	g2	AR	ASD	1/1	0/0	0/0	0/0	2	1
v4	g2	AR	ASD	0/1	0/1	0/0	0/0	2	1
v5	g3	AR	ASD	0/0	0/1	0/0	1|1	3	1
v6	g3	AR	ASD	0/1	0/0	0/0	1/0	3	1
v7	g4	AR	MD	0/0	0/0	0/0	0/1	1	0
v8	g5	AR	MD	1|1	0/0	1/0	0/0	1	1
v9	g6	AD	ASD	0/1	0/0	0/0	0/0	2	1
v10	g6	AD	ASD	1/1	0/0	./1	0/0	2	1
v11	g7	AD	ASD	1/1	0/0	0/0	0/0	2	1
v12	g7	AD	ASD	0/1	0/1	0/0	0/0	2	1
v13	g8	AD	ASD	0/0	0/1	0/0	1|1	3	1
v14	g8	AD	ASD	0/1	0/0	0/0	1/0	3	1
v15	g9	AD	MD	0/0	0/0	0/0	0/1	1	0
v16	g10	AD	MD	1|1	0/0	1/0	0/0	1	1

--->
uniq_het_for_moi_per_gene = count het based on gene+moi and cosider only where het only (./1,.|1,1/.,1|.,0/1,0|1) or if >=2 het in one sample count once only and not the overlaping het & hom in same sample. 
1. look for gene and its moi (g1+AR or g1+AD) count genotypes in each sample and put the counts in respective cols
eg
variant	gene	moi	condition	s1	s2	s3	s4	uniq_het_for_moi_per_gene
v1	g1	AR	ASD	0/1	0/0	0/0	0/0	1
v2	g1	AR	ASD	1/1	0/0	./1	0/0	1
vx	g1	AR	ASD	0/0	0/0	0/1	0/0	1
uniq_hom_for_moi_per_gene = count hom based on gene+moi and cosider only where hom only (1|1, 1/1) or if >=2 hom in one sample count once only and also count the overlaping hom & het in same sample. 
1. look for gene and its moi (g1+AR or g1+AD) count genotypes in samples and put the counts in respective cols
eg
variant	gene	moi	condition	s1	s2	s3	s4	uniq_hom_for_moi_per_gene
v1	g1	AR	ASD	0/1	0/0	0/0	0/0	1
v2	g1	AR	ASD	1/1	0/0	./1	0/0	1
vx	g1	AR	ASD	0/0	0/0	0/1	0/0	1

--->
uniq_het_for_rec_pcond_CARRIER = count het based on condition+moi and cosider only where het only (./1,.|1,1/.,1|.,0/1,0|1) or if >=2 het in one sample count once only and not the overlaping het & hom in same sample. 
For AR (moi; mode of inheritance) only:
1. look for "condition ASD (autism) + moi (AR means recessive)" count het in each sample if a sample has both hom and het then don't count het
2. if count of het is more then >1 in sample count it once for each sample 
3. sum the counts from each sample and put the counts in respective cols
eg
variant	gene	moi	condition	s1	s2	s3	s4	uniq_het_per_gene_per_moi_per_rec
v1	g1	AR	ASD	0/1	0/0	0/0	0/0	2
v2	g1	AR	ASD	1/1	0/0	./1	0/0	2
v3	g2	AR	ASD	1/1	0/0	0/0	0/0	2
v4	g2	AR	ASD	0/1	0/1	0/0	0/0	2
v5	g3	AR	ASD	0/0	0/1	0/0	1|1	2
v6	g3	AR	ASD	0/1	0/0	0/0	1/0	2
v7	g4	AR	MD	0/0	0/0	0/0	0/1	2
v8	g5	AR	MD	1|1	0/0	1/0	0/0	2
uniq_hom_for_rec_pcond_DISEASED = count hom based on condition+moi and cosider only where hom only (1|1, 1/1) or if >=2 hom in one sample count once only and also count the overlaping hom & het in same sample. 
For AR (moi; mode of inheritance) only:
1. look for "condition ASD (autism) + moi (AR means recessive)" count hom in each sample if a sample has both hom and het then count it hom
2. if count of hom is more then >1 in sample count it once for each sample 
3. sum the counts for all samples and put the counts in respective cols
eg
variant	gene	moi	condition	s1	s2	s3	s4	uniq_hom_per_gene_per_moi_per_rec
v1	g1	AR	ASD	0/1	0/0	0/0	0/0	2
v2	g1	AR	ASD	1/1	0/0	./1	0/0	2
v3	g2	AR	ASD	1/1	0/0	0/0	0/0	2
v4	g2	AR	ASD	0/1	0/1	0/0	0/0	2
v5	g3	AR	ASD	0/0	0/1	0/0	1|1	2
v6	g3	AR	ASD	0/1	0/0	0/0	1/0	2
v7	g4	AR	MD	0/0	0/0	0/0	0/1	1
v8	g5	AR	MD	1|1	0/0	1/0	0/0	1

--->
uniq_hethom_for_dom_pcond_DISEASED = count based on condition+moi if  hom or het any (./1,.|1,1/.,1|.,0/1,0|1, 1|1, 1/1) in sample if multiple hom/het are there then cosider >1 as once only 
For AD (moi; mode of inheritance) only:
1. look for "condition ASD (autism) + moi (AD means dominant)" count hom and het if both are there count hom in the specific sample if hom present >1, count it once per sample and if only het present not the hom in same sampl count het and if het present >1 in same sample count it once 
2. sum the counts for all samples and put the counts in respective cols 
eg
variant	gene	moi	condition	s1	s2	s3	s4	uniq_hethom_for_dom_pmoi_DISEASED
v9	g6	AD	ASD	0/1	0/0	0/0	0/0	4
v10	g6	AD	ASD	1/1	0/0	./1	0/0	4
v11	g7	AD	ASD	1/1	0/0	0/0	0/0	4
v12	g7	AD	ASD	0/1	0/1	0/0	0/0	4
v13	g8	AD	ASD	0/0	0/1	0/0	1|1	4
v14	g8	AD	ASD	0/1	0/0	0/0	1/0	4
v15	g9	AD	MD	0/0	0/0	0/0	0/1	1
v16	g10	AD	MD	1|1	0/0	1/0	0/0	2

--->
Step 5:
Apply Hardy-weinberg equilibrium (HWE) for calculatng Normal i.e PP (P^2) Diseased i.e QQ (Q^2) and Carrier (2PQ): P^2+Q^2+2PQ=1 (hom/hetrozygotes based) or P+Q=1 (allele frequency based i.e allele count of P =((hom*2)+het)/(popsize*2) Q((hom*2)+het))/(popsize*2)
1. Calculate popsize from sample present in genotype file
2. AR (recessive), NROMAL individauls allele count (NORMAL_AC) = 2*(popsize-(uniq_het_for_rec_pcond_CARRIER +uniq_hom_for_rec_pcond_DISEASED)), CARRIER_AC = uniq_het_for_rec_pcond_CARRIER, DISEASED_AC = 2*uniq_hom_for_rec_pcond_DISEASED
3. P=(NORMAL_AC*2)/(popsize*2) (colname = P^2_NORMAL_REC), 2PQ=CARRIER_AC/popsize (colname =2PQ_CARRIER_REC), Q^2=(1-P)^2 (colname=Q^2_DISEASED_REC)
4. AD (dominant), PREVALENCE_DOM=popsize/uniq_hethom_for_dom_pcond_DISEASED
5. put in cols for output: 2PQ_CARRIER_REC	Q^2_DISEASED_REC	P^2_NORMAL_REC	PREVALENCE_DOM
--->
Step 6:
Output file example (hint: repeat entry for based on operation if variant then fill based on then variants, if moi per pene then fill aech row satififying it with same final count)
variant	gene	moi	condition	het_raw_per_variant	hom_raw_per_variant	het_raw_per_gene	hom_raw_per_gene	uniq_het_for_moi_per_gene	uniq_hom_for_moi_per_gene	uniq_het_for_rec_pcond_CARRIER	uniq_hom_for_rec_pcond_DISEASED	uniq_hethom_for_dom_pcond_DISEASED	2PQ_CARRIER_REC		Q^2_DISEASED_REC	P^2_NORMAL_REC	PREVALENCE_DOM
v1	g1	AR	ASD	1	0	2	1	1	1	2	2	.
v2	g1	AR	ASD	1	1	2	1	1	1	2	2	.
v3	g2	AR	ASD	0	1	2	1	1	1	2	2	.
v4	g2	AR	ASD	2	0	2	1	1	1	2	2	.
v5	g3	AR	ASD	1	1	3	1	2	1	2	2	.
v6	g3	AR	ASD	2	0	3	1	2	1	2	2	.
v7	g4	AR	MD	1	0	1	0	1	0	2	1	.
v8	g5	AR	MD	1	1	1	1	1	1	2	1	.
v9	g6	AD	ASD	1	0	2	1	1	1	.	.	4
v10	g6	AD	ASD	1	1	2	1	1	1	.	.	4
v11	g7	AD	ASD	0	1	2	1	1	1	.	.	4
v12	g7	AD	ASD	2	0	2	1	1	1	.	.	4
v13	g8	AD	ASD	1	1	3	1	2	1	.	.	4
v14	g8	AD	ASD	2	0	3	1	2	1	.	.	4
v15	g9	AD	MD	1	0	1	0	1	0	.	.	3
v16	g10	AD	MD	1	1	1	1	1	1	.	.	3

# BEST ARTICLE FOR BASIC OF HWE: https://www.nature.com/scitable/definition/hardy-weinberg-equation-299/
