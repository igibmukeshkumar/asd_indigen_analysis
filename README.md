# ASD Variant Analysis Pipeline

This repository contains a step-wise bioinformatics pipeline for the identification,
annotation, statistical evaluation, and visualization of ASD-associated genetic variants
using IndiGen and gnomAD population datasets.

---

## Pipeline Overview

### Step 1: Targeted Genes Variant Identification and Annotation
**Script:** `pipeline-step1.sh`

- Preparation of targeted genes BED file
- Extraction of variants from the IndiGen dataset
- Annotation of variants
- Variant prioritization based on:
  - **REVEL score** (missense variants)
  - **LOFTEE** (loss-of-function variants)
- Variant classification using **ACMG–AMP guidelines**

---

### Step 2: gnomAD Variant Query and File Preparation
**Scripts (VCF-based):**
- `varID_query_gnomAD.sh`
- `varID_joint_flatten_gnomad.sh`

**Scripts (Gene-based):**
- `query_gnomAD.sh`
- `joint_flatten_gnomad.sh`

- Query gnomAD API using variant IDs or gene-based inputs
- Retrieve allele counts and allele numbers for the available populations
- Generate flattened, analysis-ready tables from gnomAD responses

---

### Step 3: Population-Level MAF (minor allele frquency) Estimation Analysis Across the Globe
**Script:** `maf_fisher_test.R`

- Computation of population-specific MAF
- Construction of contingency tables
- Statistical testing using **Fisher’s Exact Test**
- Multiple testing correction using **Bonferroni adjustment**
- Identification of variants significantly enriched in IndiGen compared to global populations

---

### Step 4: Visualization for Research Publication
**Script:** `data_plot.R`

- Integration of MAF, statistical significance, ACMG classification, and gene annotation
- Generation of publication-quality population comparison plots
- Figures used directly in the research article

  
---

### Step 5: Carrier Frequency and Prevalence Estimation
- Calculation of carrier frequency and disease prevalence
- Focused analysis on ASD-associated variants and genes
- IndiGen cohort–specific burden estimation

####  ***Carrier and Prevalence Rate Estimation***

To analyze genotype data (GT matrix) and compute:

* Per-variant and per-gene raw genotype counts (het/hom)
* Unique per-gene/per-condition/per-MOI statistics
* Carrier and disease prevalence under **autosomal recessive (AR)** and **autosomal dominant (AD)** inheritance models
* Hardy-Weinberg Equilibrium estimates

---

#### ⚙️ **Step-by-Step Plan**

### ✅ **Assumptions**

1. A **condition** (e.g., WD) may involve multiple variants in one gene.
2. A **condition** (e.g., ASD) may involve multiple **genes**.
3. A **gene** (e.g., *FMR1*) can cause multiple **conditions** — **pleiotropy**.

---

##### 🧩 **Step 5.1: Script Interface**

Write a shell wrapper: `cr_pr.sh`, which calls a Python script:

```bash
Usage: ./cr_pr.sh -gf <GenotypeFile> -qf <QueryFile> -o <OutputFile>
  -gf, --GTfile       Path to genotype matrix file
  -qf, --Queryfile    Path to query file (variant-gene-moi-condition)
  -o,  --Output       Output file path
  -h,  --help         Display this help message
```

---

##### 📄 **Step 5.2: Input File Format**

**Genotype File**

```
variant	s1	s2	s3	s4
v1	    0/1	0/0	0/0	0/0
...
```

**Query File**

```
variant	gene	moi	condition
v1	    g1	  AR	  ASD
...
```

---

##### 🔢 **Step 5.3: Data Integration**

* Match `variant` entries between genotype and query file.
* Add `gene`, `moi`, and `condition` as metadata to each variant row.

---

##### 📊 **Step 5.4: Genotype Counting**

###### A. **Per Variant (Row-wise)**

* `het_raw_per_variant`: count of `0/1, 1/0, 0|1, 1|0, ./1, 1/.`
* `hom_raw_per_variant`: count of `1/1, 1|1`

###### B. **Per Gene**

* `het_raw_per_gene`: total hets across all variants in the same gene.
* `hom_raw_per_gene`: total homs across all variants in the same gene.

###### C. **Unique Per Gene+MOI**

* `uniq_het_for_moi_per_gene`: 1 count per sample if at least one **het** only.
* `uniq_hom_for_moi_per_gene`: 1 count per sample if at least one **hom**, even if both het & hom.

---

##### 🧬 **Step 5.5: Disease Model-Based Aggregation**

###### For **AR (Recessive)**

* `uniq_het_for_rec_pcond_CARRIER`: total unique **carrier** samples per condition.
* `uniq_hom_for_rec_pcond_DISEASED`: total unique **disease** samples per condition.

Rules:

* If sample has both het and hom, count it as **hom** only.
* Count only once per sample even if multiple het/hom variants.

###### For **AD (Dominant)**

* `uniq_hethom_for_dom_pcond_DISEASED`: count samples with **any het or hom** in a gene.
* If both are present in a sample → count once as **hom** (dominant disease logic).

---

##### 🧮 **Step 5.6: Hardy-Weinberg Equilibrium (HWE) Calculations**

**Given:**

* popsize = number of samples in genotype file

###### For **AR conditions**:

* `NORMAL_AC` = 2 × (popsize − carriers − diseased)
* `CARRIER_AC` = number of **carriers**
* `DISEASED_AC` = 2 × number of **diseased**

Then:

* **P²\_NORMAL\_REC** = (NORMAL\_AC)² / (2 × popsize)²
* **2PQ\_CARRIER\_REC** = CARRIER\_AC / popsize
* **Q²\_DISEASED\_REC** = (1 - P)²

###### For **AD conditions**:

* **PREVALENCE\_DOM** = (uniq\_hethom\_for\_dom\_pcond\_DISEASED) / popsize

---

##### 🧾 **Step 5.7: Output Format**

All columns summarized for downstream use:

| variant | gene | moi | condition | het\_raw\_per\_variant | hom\_raw\_per\_variant | het\_raw\_per\_gene | hom\_raw\_per\_gene | uniq\_het\_for\_moi\_per\_gene | uniq\_hom\_for\_moi\_per\_gene | uniq\_het\_for\_rec\_pcond\_CARRIER | uniq\_hom\_for\_rec\_pcond\_DISEASED | uniq\_hethom\_for\_dom\_pcond\_DISEASED | 
| ------- | ---- | --- | --------- | ---------------------- | ---------------------- | ------------------- | ------------------- | ------------------------------ | ------------------------------ | ----------------------------------- | ------------------------------------ | --------------------------------------- |
