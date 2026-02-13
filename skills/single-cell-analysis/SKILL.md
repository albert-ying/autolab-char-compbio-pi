---
name: single-cell-analysis
description: End-to-end single-cell RNA-seq analysis with scanpy and scvi-tools. QC, normalization, batch correction, clustering, differential expression, trajectory inference.
metadata:
    skill-author: Albert Ying
---

# Single-cell analysis

## When to use

- Analyzing scRNA-seq data (.h5ad, 10X, CSV)
- QC filtering, normalization, batch correction
- Clustering, UMAP, marker gene identification
- Cell type annotation and differential expression

## Standard workflow

```python
import scanpy as sc
import scvi

# Load and QC
adata = sc.read_h5ad("data.h5ad")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
adata = adata[adata.obs.pct_counts_mt < 20].copy()

# Normalize and select HVGs
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="batch")

# scVI for batch correction (preferred over Harmony for complex designs)
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
model = scvi.model.SCVI(adata, n_latent=30)
model.train(max_epochs=200, early_stopping=True)
adata.obsm["X_scVI"] = model.get_latent_representation()

# Cluster and visualize
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.8)
```

## Differential expression

```python
# Wilcoxon rank-sum (default, robust)
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)

# scVI DE (accounts for batch, preferred for complex designs)
de_results = model.differential_expression(
    adata, groupby="cell_type", group1="CD4_T", group2="rest"
)
de_results = de_results[de_results["is_de_fdr_0.05"]]
```

## Key decisions

- **Batch correction**: scVI > Harmony > combat for most cases. Use scVI when batch effects are strong or when you need a generative model downstream.
- **Resolution**: start at 0.8, adjust based on known biology. Over-clustering is preferable to under-clustering (you can merge later).
- **HVG selection**: 2000-3000 genes. Use `batch_key` if multiple batches.
- **QC thresholds**: adapt to tissue type. Brain tissue tolerates higher MT% than PBMCs.
