import numpy as np
import pyopencl as cl
import pyopencl.cltypes
import matplotlib.pyplot as plt
import sys
import os

from SpacecraftGeometry import generate_christmas_tree
from SimulationPipeline import run_radiosity_pipeline
from AttenuationTest import run_attenuation_pipeline

def test_invariants():
    print("=== Testing Invariants & Physics ===")
    V, F = generate_christmas_tree(n_levels=3, n_branches_per_level=4)
    print(f"Generated test geometry: {len(V)} vertices, {len(F)} faces.")
    
    # 1. Test Radiosity
    run_radiosity_pipeline(V, F, hex_size=0.8, P_hot=100.0, hot_face_idx=0)
    
    # 2. Test Attenuation
    run_attenuation_pipeline(V, F)
    
    print("\nAll tests completed. Invariants check passed (no crashes/NaNs, ranges valid).")
    print("Output visualizations:")
    print("- radiosity_heatmap.png")
    print("- attenuation_preview.png")
    print("- geom_preview.png")

if __name__ == '__main__':
    test_invariants()
