#!/bin/bash

################################################################################
# BCFtools De Novo Variant Analysis for Chr1
# 
# This script uses bcftools isec to identify ClinVar variants present ONLY 
# in the child (HG002) and not in either parent (HG003=father, HG004=mother).
# These are candidate de novo mutations.
#
# Usage: bash run_bcftools_denovo_analysis.sh
# Must be run from: wf-trio-variant-pipeline directory
################################################################################

# Exit on error
set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print formatted messages
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

################################################################################
# Configuration
################################################################################

# Base directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTPUT_DIR="${SCRIPT_DIR}"
BASE_DIR="/mnt/work_1/gest9386/CU_Boulder/MCDB-4520/bioinf-group-project"
DATA_DIR="${BASE_DIR}/data"

# Input ClinVar VCF files
HG002_VCF="${DATA_DIR}/HG002_chr1/output/HG002_chr1.wf_snp_clinvar.vcf.gz"
HG003_VCF="${DATA_DIR}/HG003_chr1/output/HG003_chr1.wf_snp_clinvar.vcf.gz"
HG004_VCF="${DATA_DIR}/HG004_chr1/output/HG004_chr1.wf_snp_clinvar.vcf.gz"

# Container path - using apptainer cache
CONTAINER_IMAGE="${HOME}/.apptainer/cache/ontresearch-wf-human-variation-sha8ecee6d351b0c2609b452f3a368c390587f6662d.img"

# Detect container engine (apptainer or singularity)
if command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
elif command -v singularity &> /dev/null; then
    CONTAINER_CMD="singularity"
else
    log_error "Neither apptainer nor singularity found!"
    exit 1
fi

################################################################################
# Main Analysis
################################################################################

log_info "Starting BCFtools de novo variant analysis for chromosome 1"
log_info "Using container engine: ${CONTAINER_CMD}"
echo ""

# Check if input files exist
log_info "Checking input files..."
if [[ ! -f "$HG002_VCF" ]]; then
    log_error "HG002 ClinVar VCF not found: $HG002_VCF"
    exit 1
fi
if [[ ! -f "$HG003_VCF" ]]; then
    log_error "HG003 ClinVar VCF not found: $HG003_VCF"
    exit 1
fi
if [[ ! -f "$HG004_VCF" ]]; then
    log_error "HG004 ClinVar VCF not found: $HG004_VCF"
    exit 1
fi
log_success "All input VCF files found"
echo ""

# Check if container exists
log_info "Checking container image..."
if [[ ! -f "$CONTAINER_IMAGE" ]]; then
    log_error "Container image not found: $CONTAINER_IMAGE"
    log_info "Please ensure the wf-human-variation container is available"
    exit 1
fi
log_success "Container image found"
echo ""

# Create output directory
log_info "Creating output directory: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"
echo ""

# Run bcftools isec to find de novo variants
log_info "Running bcftools isec to identify de novo variants..."
log_info "This will find variants present ONLY in HG002 (child)"
echo ""

${CONTAINER_CMD} exec \
    --bind /mnt/work_1/gest9386/CU_Boulder/MCDB-4520:/mnt/work_1/gest9386/CU_Boulder/MCDB-4520 \
    "$CONTAINER_IMAGE" \
    bcftools isec -n=1 \
    "$HG002_VCF" \
    "$HG003_VCF" \
    "$HG004_VCF" \
    -p "$OUTPUT_DIR"

log_success "BCFtools isec completed successfully"
echo ""

# Display results
log_info "Output files created:"
ls -lh "$OUTPUT_DIR"
echo ""

log_info "File descriptions:"
echo "  - 0000.vcf      : ClinVar variants ONLY in HG002 (son) - DE NOVO CANDIDATES"
echo "  - 0001.vcf      : ClinVar variants ONLY in HG003 (father)"
echo "  - 0002.vcf      : ClinVar variants ONLY in HG004 (mother)"
echo "  - sites.txt     : Tab-delimited list showing which file(s) contain each variant"
echo "  - README.txt    : Description of the output files"
echo ""

# Count de novo variants
log_info "Counting de novo variants in child (excluding VCF header)..."
DENOVO_COUNT=$(grep -v "^#" "$OUTPUT_DIR/0000.vcf" | wc -l)
log_success "Found $DENOVO_COUNT de novo candidate variants in HG002"
echo ""

# Search for pathogenic variants
log_info "Searching for pathogenic variants in de novo candidates..."
echo ""
echo "========================================================================="
echo "PATHOGENIC DE NOVO VARIANTS IN HG002 (CHILD)"
echo "========================================================================="
if grep -i "pathogenic" "$OUTPUT_DIR/0000.vcf"; then
    echo ""
    log_success "Found pathogenic variants (see above)"
else
    log_warning "No pathogenic variants found in de novo candidates"
fi
echo ""

# Create a summary report
SUMMARY_FILE="${OUTPUT_DIR}/denovo_analysis_summary.txt"
log_info "Creating summary report: $SUMMARY_FILE"

cat > "$SUMMARY_FILE" << EOF
================================================================================
BCFtools De Novo Variant Analysis Summary - Chromosome 1
================================================================================

Analysis Date: $(date)
Script: run_bcftools_denovo_analysis.sh
Container Engine: ${CONTAINER_CMD}

INPUT FILES:
  - Child (HG002)  : $HG002_VCF
  - Father (HG003) : $HG003_VCF
  - Mother (HG004) : $HG004_VCF

OUTPUT DIRECTORY: $OUTPUT_DIR

RESULTS:
  - Total de novo candidate variants in HG002: $DENOVO_COUNT
  - These are ClinVar variants present ONLY in the child
  - Not present in either parent

PATHOGENIC VARIANTS:
EOF

# Add pathogenic variants to summary
if grep -i "pathogenic" "$OUTPUT_DIR/0000.vcf" > /dev/null 2>&1; then
    echo "  Found pathogenic variants (see below):" >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"
    grep -i "pathogenic" "$OUTPUT_DIR/0000.vcf" >> "$SUMMARY_FILE"
else
    echo "  No pathogenic de novo variants found" >> "$SUMMARY_FILE"
fi

log_success "Summary report created: $SUMMARY_FILE"
echo ""

log_success "Analysis complete!"
log_info "You can now examine the de novo variants in: $OUTPUT_DIR/0000.vcf"
echo ""

# Additional statistics
log_info "Additional Statistics:"
echo ""
echo "Father-specific ClinVar variants (HG003):"
grep -v "^#" "$OUTPUT_DIR/0001.vcf" | wc -l
echo ""
echo "Mother-specific ClinVar variants (HG004):"
grep -v "^#" "$OUTPUT_DIR/0002.vcf" | wc -l
echo ""

log_info "For detailed analysis of pathogenic variants, review:"
log_info "  1. $OUTPUT_DIR/0000.vcf (de novo candidates)"
log_info "  2. $SUMMARY_FILE (summary report)"
echo ""

# Display first few pathogenic variants if any exist
if grep -i "pathogenic" "$OUTPUT_DIR/0000.vcf" > /dev/null 2>&1; then
    log_info "Preview of pathogenic variants:"
    echo "========================================================================="
    grep -i "pathogenic" "$OUTPUT_DIR/0000.vcf" | head -10
    echo "========================================================================="
fi
