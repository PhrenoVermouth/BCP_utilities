
import scanpy as sc
import anndata as ad
import os
import argparse
import warnings
from samalg import SAM
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Suppress warnings as seen in the user logs (optional but keeps output clean)
warnings.filterwarnings('ignore')

def main():
    parser = argparse.ArgumentParser(description="Merge h5ad files and run SAM preprocessing.")
    parser.add_argument('samples', metavar='S', type=str, nargs='+',
                        help='List of sample names to merge (e.g., hbm hbf mbm)')
    parser.add_argument('-s', '--species', type=str, default='q',
                        help='Species prefix (default: q)')
    parser.add_argument('-o', '--output', type=str,
                        help='Output filename (default: <species>_brain_bc.h5ad)')
    
    args = parser.parse_args()
    samples = args.samples
    species = args.species
    
    # Set default output filename if not provided
    if args.output:
        output_filename = args.output
    else:
        output_filename = f"{species}_brain_bc.h5ad"
    
    adatas = []
    
    print(f"Species prefix: {species}")
    print("开始加载数据...")
    for sample in samples:
        # Assuming the directory structure follows the pattern:
        # {species}_{sample}/{species}_{sample}_filtered_QC2.h5ad
        file_path = f"{species}_{sample}/{species}_{sample}_filtered_QC2.h5ad"
        
        if os.path.exists(file_path):
            print(f"Loading {file_path}...")
            try:
                adata = sc.read_h5ad(file_path)
                
                # -------------------------------------------------
                # 改动 1: 手动给 index 加上前缀
                # 格式变成: "样本名-Barcode" (例如: hbm-AAACCTGAG...)
                # -------------------------------------------------
                adata.obs_names = sample + "_" + adata.obs_names
                
                # 保留 sample_id 列
                adata.obs['sample_id'] = sample
                
                adatas.append(adata)
            except Exception as e:
                 print(f"Error loading {file_path}: {e}")
        else:
            print(f"Warning: {file_path} not found!")

    if not adatas:
        print("No valid data found. Exiting.")
        return

    print("正在合并...")
    
    # -------------------------------------------------
    # 改动 2: 去掉 index_unique 参数
    # 因为我们在上面已经把 index 改成唯一的了，这里不需要再自动添加后缀
    # -------------------------------------------------
    adata_combined = ad.concat(adatas,
                               join='outer',
                               label='batch',
                               keys=samples,
                               index_unique=None) 
    
    print(f"Merged Anndata shape: {adata_combined.shape}")
    
    # SAM Preprocessing
    print("Starting SAM preprocessing...")
    
    # Rename leiden_clusters to leiden_clusters_raw if it exists, to match user workflow
    if 'leiden_clusters' in adata_combined.obs:
         # Use assignment and delete to ensure we end up with one 'leiden_clusters_raw' column
         # and no 'leiden_clusters' column, avoiding duplicates even if 'leiden_clusters_raw' existed.
         adata_combined.obs['leiden_clusters_raw'] = adata_combined.obs['leiden_clusters']
         del adata_combined.obs['leiden_clusters']

    
    # SAM Parameters
    n_neighbors = 20
    n_pcs = 150
    n_hvg = 3000
    
    sam = SAM(adata_combined)
    
    print("Running SAM preprocess_data...")
    sam.preprocess_data()
    
    print("Running SAM run...")
    sam.run(
        projection='umap',
        weight_mode='rms',
        k=n_neighbors,
        npcs=n_pcs,
        n_genes=n_hvg,
        batch_key='batch'
    )
    
    print(f"Saving results to {output_filename}...")
    try:
        sam.adata.write_h5ad(output_filename)
        print("Success!")
    except Exception as e:
        print(f"Error saving file: {e}")

if __name__ == "__main__":
    main()
