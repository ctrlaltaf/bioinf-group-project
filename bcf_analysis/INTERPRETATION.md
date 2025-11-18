# BCFtools De Novo Variant Analysis - Interpretation

## Analysis Overview
**Date:** November 17, 2025  
**Chromosome:** chr1  
**Method:** bcftools isec  
**Total de novo candidates:** 70 ClinVar variants unique to HG002 (child)

## Key Findings

### 1. De Novo Pathogenic/Conflicting Variants Found: 3

#### Variant 1: DBT Gene (Maple Syrup Urine Disease)
- **Position:** chr1:100189351
- **Change:** TA → T (deletion)
- **Gene:** DBT (Dihydrolipoamide Branched Chain Transacylase E2)
- **Location:** 3' UTR variant
- **Clinical Significance:** **Conflicting classifications of pathogenicity**
  - Uncertain significance (1)
  - Benign (1)
- **Associated Condition:** Maple syrup urine disease (OMIM:PS248600)
- **Genotype:** Homozygous (1/1)
- **Quality Score:** 12.45
- **Allele Frequency in Child:** 28.89%

**Interpretation:** This variant is in the 3' UTR (non-coding region) of the DBT gene. Maple syrup urine disease is a serious metabolic disorder, but this specific variant has conflicting interpretations. The variant being homozygous and in a non-coding region suggests it may not be directly pathogenic.

---

#### Variant 2: SPTA1 Gene (Hereditary Spherocytosis/Elliptocytosis)
- **Position:** chr1:158668075
- **Change:** GA → G (deletion)
- **Gene:** SPTA1 (Spectrin Alpha, Erythrocytic 1)
- **Location:** Intron variant (c.1834-14delT)
- **Clinical Significance:** **Conflicting classifications of pathogenicity**
  - Uncertain significance (3)
  - Benign (1)
- **Associated Conditions:**
  - Spherocytosis, Recessive (OMIM:266140)
  - Pyropoikilocytosis, hereditary
  - Elliptocytosis
- **Genotype:** Homozygous (1/1)
- **Quality Score:** 24.7
- **Allele Frequency in Child:** 40.38%

**Interpretation:** This intronic variant is 14 bases upstream of an exon. SPTA1 mutations cause red blood cell membrane disorders leading to hemolytic anemia. However, intronic variants far from splice sites (>10bp) are generally less likely to affect splicing. The conflicting classifications and intronic location suggest this may be benign.

---

#### Variant 3: FH Gene (Fumarase Deficiency/Cancer Syndrome)
- **Position:** chr1:241500602
- **Change:** T → TGA (insertion/duplication)
- **Gene:** FH (Fumarate Hydratase)
- **Location:** Intron variant (microsatellite repeat)
- **Clinical Significance:** **Conflicting classifications of pathogenicity**
  - Uncertain significance (4)
  - Benign (5)
  - Likely benign (1)
- **Associated Conditions:**
  - Hereditary leiomyomatosis and renal cell cancer (HLRCC, OMIM:150800)
  - Fumarase deficiency (OMIM:606812)
  - Hereditary cancer-predisposing syndrome
- **Genotype:** Heterozygous (1|0)
- **Quality Score:** 30.77
- **Allele Frequency in Child:** 35.09%

**Interpretation:** This is a microsatellite (GA repeat) variant in an intron of the FH gene. While FH mutations can cause serious conditions including cancer predisposition and metabolic disorders, this variant has been classified as benign or likely benign by the majority of submitters (6 out of 10). Microsatellite variations in introns are commonly polymorphic.

---

## Summary and Conclusions

### Statistical Context
- **Child-specific variants (de novo candidates):** 70
- **Father-specific ClinVar variants:** 476
- **Mother-specific ClinVar variants:** 537

### Clinical Interpretation

1. **No Clearly Pathogenic De Novo Variants Identified**
   - All three variants with "pathogenic" annotations have **conflicting classifications**
   - None are classified as "Pathogenic" or "Likely Pathogenic" by ClinVar consensus

2. **Variant Locations Suggest Low Clinical Impact**
   - All three variants are in **non-coding regions** (UTR or introns)
   - Distance from splice sites and coding regions reduces functional impact likelihood
   - Two variants are homozygous, one is heterozygous

3. **Quality Considerations**
   - Quality scores range from 12.45 to 30.77 (moderate quality)
   - Allele frequencies are relatively low to moderate (28-40%)
   - These may represent sequencing artifacts or mosaic variants

### Speculation & Further Analysis Needed

1. **Possible Explanations for "De Novo" Variants:**
   - **False positives:** Low coverage in parents leading to missed variants
   - **Somatic mosaicism in child:** Post-zygotic mutations
   - **Germline mosaicism in parents:** Low-level variants not detected
   - **Technical artifacts:** Sequencing or alignment errors specific to child sample

2. **Clinical Recommendations:**
   - **Sanger sequencing validation** recommended for all three variants
   - **Parental re-sequencing** at higher coverage to confirm absence
   - **Clinical correlation:** Does the child show any symptoms related to:
     - Metabolic disorders (DBT/MSUD)
     - Hemolytic anemia (SPTA1)
     - Renal tumors or skin leiomyomas (FH/HLRCC)

3. **Functional Studies:**
   - RNA expression analysis to check if intronic variants affect splicing
   - Family segregation studies if extended family available
   - Long-read sequencing to better characterize microsatellite regions

### Conclusion

**No definitively pathogenic de novo variants were identified in chromosome 1.** The three variants with conflicting pathogenicity classifications are all in non-coding regions and are more likely benign based on:
- Conflicting/uncertain ClinVar classifications
- Non-coding locations (UTR/introns)
- Distance from functional elements

The child (HG002) appears to have a normal chromosomal 1 profile without concerning de novo pathogenic variants. However, validation studies are recommended to confirm these findings, especially given the clinical importance of the genes involved.

---

## Files Generated
- `0000.vcf` - All 70 de novo candidate variants in HG002
- `0001.vcf` - Father-specific variants (476)
- `0002.vcf` - Mother-specific variants (537)
- `sites.txt` - Tab-delimited variant intersection summary
- `denovo_analysis_summary.txt` - Technical analysis summary
- `INTERPRETATION.md` - This clinical interpretation document

## Methods
Analysis performed using bcftools isec (v1.19) within the ONT wf-human-variation Singularity/Apptainer container. Variants were filtered to include only those present in the child (HG002) and absent from both parents (HG003, HG004) using the `-n=1` flag.
