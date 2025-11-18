#!/usr/bin/env python3
"""Debug script to inspect AlphaGenome splice junction output structure."""

import sys
import os
sys.path.insert(0, '/opt/alphagenome')

from alphagenome.data import genome
from alphagenome.models import dna_client
import numpy as np

# Initialize
api_key = os.environ.get('ALPHA_GENOME_KEY')
if not api_key:
    print("ERROR: ALPHA_GENOME_KEY not set")
    sys.exit(1)

print("Initializing AlphaGenome client...")
model = dna_client.create(api_key)

# Test variant (CEP104 - same as in our data)
interval = genome.Interval(chromosome='chr1', start=3312016, end=4360592)
variant = genome.Variant(chromosome='chr1', position=3836304, reference_bases='CAA', alternate_bases='C')

print("Requesting splice junction predictions...")
outputs = model.predict_variant(
    interval=interval,
    variant=variant,
    requested_outputs=[dna_client.OutputType.SPLICE_JUNCTIONS],
    ontology_terms=[]
)

print("\n" + "="*80)
print("INSPECTING OUTPUT STRUCTURE")
print("="*80)

# Reference
ref_sj = outputs.reference.splice_junctions
print("\nReference splice_junctions object:")
print(f"  Type: {type(ref_sj)}")
print(f"  Dir: {[attr for attr in dir(ref_sj) if not attr.startswith('_')]}")

# Check common attribute names
attrs_to_check = [
    'donor_probabilities', 'acceptor_probabilities',
    'donor', 'acceptor',
    'values', 'probabilities',
    'interval', 'chromosome', 'start', 'end'
]

print("\nAttribute availability:")
for attr in attrs_to_check:
    has_it = hasattr(ref_sj, attr)
    print(f"  {attr}: {has_it}")
    if has_it:
        val = getattr(ref_sj, attr)
        print(f"    Type: {type(val)}")
        if isinstance(val, np.ndarray):
            print(f"    Shape: {val.shape}")
            print(f"    Dtype: {val.dtype}")
            print(f"    Min/Max: {val.min():.6f} / {val.max():.6f}")
            print(f"    Mean: {val.mean():.6f}")
            print(f"    First 10 values: {val[:10]}")

# Alternate
print("\n" + "-"*80)
alt_sj = outputs.alternate.splice_junctions
print("\nAlternate splice_junctions object:")
print(f"  Type: {type(alt_sj)}")

# If there are values, compare them
if hasattr(ref_sj, 'values') and hasattr(alt_sj, 'values'):
    print("\n" + "="*80)
    print("COMPARING VALUES")
    print("="*80)
    ref_vals = ref_sj.values
    alt_vals = alt_sj.values
    diff = np.abs(alt_vals - ref_vals)
    print(f"Reference values shape: {ref_vals.shape}")
    print(f"Alternate values shape: {alt_vals.shape}")
    print(f"Max absolute difference: {diff.max():.6f}")
    print(f"Mean absolute difference: {diff.mean():.6f}")
    
    # Find position of variant
    variant_idx = variant.position - ref_sj.interval.start
    print(f"\nVariant position in array: {variant_idx}")
    if 0 <= variant_idx < len(ref_vals):
        window = 10
        start = max(0, variant_idx - window)
        end = min(len(ref_vals), variant_idx + window + 1)
        print(f"Values around variant (±{window} bp):")
        print(f"  Position range: {start} to {end}")
        print(f"  Reference: {ref_vals[start:end]}")
        print(f"  Alternate: {alt_vals[start:end]}")
        print(f"  Difference: {diff[start:end]}")
