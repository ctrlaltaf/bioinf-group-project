# MCDB-5520 Group Project: De Novo Variant Functional Validation

Computational genomics analysis pipeline for identifying and functionally validating de novo variants in trio sequencing data (HG002/HG003/HG004) using Google DeepMind's AlphaGenome multimodal genomic AI.

## Overview

This project combines traditional variant calling approaches with cutting-edge AI-powered functional predictions to:
1. Identify de novo variants from trio VCF data
2. Validate functional impact using AlphaGenome's multimodal predictions
3. Resolve conflicting ClinVar classifications with functional evidence

## Key Results

✅ **55/70 variants successfully validated** with AlphaGenome functional predictions

**Major Findings:**
- 🔴 **5 potentially pathogenic variants** misclassified as benign in ClinVar (TARDBP, CASQ2, LAMC1, LAMC2, TGFB2)
- ⚠️ **3 conflicting ClinVar classifications** partially resolved using functional evidence (DBT, SPTA1, FH)
- 📊 **Comprehensive functional evidence** across RNA-seq, splice junctions, and chromatin accessibility

**Classification Distribution:**
- Uncertain significance: 38 variants (69%)
- Likely benign: 11 variants (20%)
- Likely pathogenic: 6 variants (11%)

## Project Structure

```
.
├── alphagenome_validation/          # AlphaGenome validation pipeline
│   ├── outputs/                     # 55 variant prediction results (JSON + CSV summary)
│   └── scripts/                     # Validation scripts
│       ├── validate_variants.py     # Main validation pipeline
│       ├── all_variants_data.py     # 70 variant definitions
│       └── visualize_predictions.py # Visualization utilities
├── analysis_results/                # Trio variant analysis results
│   ├── de_novo_candidate_variants.csv
│   ├── *_pathogenic_variants.csv    # ClinVar pathogenic variants per sample
│   └── variant_summary_statistics.csv
├── data/                            # VCF data (not in repo - too large)
│   ├── HG002_chr1/                  # Son (proband)
│   ├── HG003_chr1/                  # Father
│   └── HG004_chr1/                  # Mother
├── run_bcftools_denovo_analysis.sh  # De novo variant identification
└── analysis.ipynb                   # Exploratory analysis notebook
```

## AlphaGenome Validation Pipeline

### What is AlphaGenome?

AlphaGenome is Google DeepMind's multimodal genomic AI that predicts functional effects of genetic variants using:
- **RNA-seq**: Gene expression changes
- **Splice Junctions**: Splicing disruption
- **ATAC**: Chromatin accessibility changes

### Methodology

For each variant, AlphaGenome compares:
- **Reference genome** (wild-type sequence)
- **Alternate genome** (with variant inserted)

And predicts functional impacts across a 1 Mbp window centered on the variant.

### Running the Pipeline

```bash
# Setup environment
conda activate alphagenome-env

# Set API key
export ALPHA_GENOME_KEY="your_api_key"

# Validate all 70 variants
cd alphagenome_validation/scripts
python3 validate_variants.py

# Or validate specific variant
python3 validate_variants.py --variant FH

# Results saved to: alphagenome_validation/outputs/
```

### Key Scripts

**`validate_variants.py`** - Main validation pipeline
- Loads 70 de novo variants from `all_variants_data.py`
- Makes AlphaGenome API predictions for RNA-seq, splice junctions, ATAC
- Assesses pathogenicity based on functional evidence
- Outputs individual JSON results + CSV summary

**`visualize_predictions.py`** - Visualization tools
- Generate plots for RNA-seq expression changes
- Visualize splice junction disruption
- Plot chromatin accessibility differences

## De Novo Variant Analysis

### Variant Identification

Used bcftools to identify candidate de novo variants:

```bash
# Run de novo analysis
./run_bcftools_denovo_analysis.sh

# Output: analysis_results/de_novo_candidate_variants.csv
```

**Filtering Criteria:**
- Present in child (HG002) with high quality
- Absent in both parents (HG003/HG004)
- Minimum depth 10x, quality score 30+
- Focus on chromosome 1 for manageable analysis

### Variant Selection for AlphaGenome

From 70 identified de novo variants:
- **Intronic variants** (n=35): Analyzed with SPLICE_JUNCTIONS + ATAC
- **UTR variants** (n=18): Analyzed with RNA_SEQ + ATAC  
- **Coding variants** (n=12): Analyzed with RNA_SEQ + SPLICE_JUNCTIONS
- **Other** (n=5): Context-specific analysis strategies

## Installation & Setup

### Requirements

```bash
# Conda environment
conda create -n alphagenome-env python=3.11
conda activate alphagenome-env

# AlphaGenome (requires /opt/alphagenome installation)
pip install numpy pandas matplotlib seaborn

# bcftools for variant analysis
conda install -c bioconda bcftools
```

### AlphaGenome API Access

1. Request API key from Google DeepMind AlphaGenome team
2. Set environment variable: `export ALPHA_GENOME_KEY="your_key"`
3. Test connection: See `alphagenome_validation/scripts/debug_splice_junctions.py`

## Notable Discoveries

### 1. FH (Fumarate Hydratase) - chr1:241500602

**ClinVar:** Conflicting classifications  
**AlphaGenome:** Likely_pathogenic (Moderate confidence)

**Evidence:**
- Intronic microsatellite deletion
- Local splice junction change: 0.115 (>0.1 threshold)
- Strong evidence for splicing disruption
- **Interpretation:** Supports pathogenic classification

### 2. TARDBP (TDP-43) - chr1:11017202

**ClinVar:** Benign  
**AlphaGenome:** Likely_pathogenic (Moderate confidence)

**Evidence:**
- Intronic deletion affecting splice region
- Splice junction probability change: 0.127
- TARDBP mutations associated with ALS/FTD
- **Interpretation:** May warrant reclassification

### 3. CASQ2 (Calsequestrin 2) - chr1:115725557

**ClinVar:** Benign/Likely_benign  
**AlphaGenome:** Likely_pathogenic (Moderate confidence)

**Evidence:**
- Splice region variant
- Strong splicing disruption signal (0.134)
- CASQ2 mutations cause catecholaminergic polymorphic VT
- **Interpretation:** Functional evidence contradicts benign classification

## Limitations

1. **Multi-allelic variants** - AlphaGenome API doesn't support complex alleles (15/70 failed)
2. **Chromosome 1 only** - Limited to chr1 for computational feasibility
3. **AI predictions** - Require experimental validation to confirm functional effects
4. **API rate limits** - Processing 70 variants takes ~1-2 hours

## Future Directions

- [ ] Extend analysis to all chromosomes
- [ ] Experimental validation of top discordant variants (TARDBP, CASQ2, LAMC1/2)
- [ ] Integrate with AlphaMissense for coding variant predictions
- [ ] Compare AlphaGenome vs traditional in silico tools (CADD, SpliceAI)
- [ ] Submit evidence to ClinVar for reclassification

## References

- **AlphaGenome**: Google DeepMind's genomic AI for functional predictions
- **ClinVar**: NCBI database of variant-disease relationships
- **GIAB Trio**: Genome in a Bottle reference trio (HG002/HG003/HG004)

## Authors

MCDB-5520 Computational Genomics Group Project  
University of Colorado Boulder  
Fall 2024

## License

Academic use only - See course materials for details
