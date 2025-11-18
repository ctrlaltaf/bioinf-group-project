#!/usr/bin/env python3
"""
Quick test of validate_variants.py on first 5 variants to verify fixes.
"""

import sys
import os
from pathlib import Path

# Add current script directory to path
script_dir = Path(__file__).parent.absolute()
sys.path.insert(0, str(script_dir))

# Import and modify
import all_variants_data
original_variants = all_variants_data.ALL_VARIANTS.copy()
all_variants_data.ALL_VARIANTS = original_variants[:5]

print(f"Testing with {len(all_variants_data.ALL_VARIANTS)} variants (first 5):")
for v in all_variants_data.ALL_VARIANTS:
    print(f"  - {v['name']} ({v['gene']})")
print()

# Now run the main validation script
import validate_variants
validate_variants.main()
