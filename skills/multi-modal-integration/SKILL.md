---
name: multi-modal-integration
description: Multi-modal single-cell data integration (RNA + ATAC, RNA + protein). TOTALVI, MultiVI, mosaic integration approaches.
metadata:
    skill-author: Albert Ying
---

# Multi-modal integration

## When to use

- Integrating scRNA-seq with scATAC-seq
- Analyzing CITE-seq (RNA + surface protein)
- Mosaic integration (not all modalities in all cells)
- Joint embedding across modalities

## CITE-seq with TOTALVI

```python
import scvi

# Setup with protein counts
scvi.model.TOTALVI.setup_anndata(
    adata,
    protein_expression_obsm_key="protein_expression",
    layer="counts",
    batch_key="batch",
)
model = scvi.model.TOTALVI(adata)
model.train(max_epochs=200, early_stopping=True)

# Joint latent space
adata.obsm["X_totalVI"] = model.get_latent_representation()

# Denoised protein expression
_, protein_means = model.get_normalized_expression(return_mean=True)
```

## RNA + ATAC with MultiVI

```python
# Requires paired or mosaic data
scvi.model.MULTIVI.setup_anndata(adata, batch_key="batch")
model = scvi.model.MULTIVI(adata)
model.train()
adata.obsm["X_multiVI"] = model.get_latent_representation()
```

## Key decisions

- **TOTALVI** for CITE-seq. Handles protein background noise properly.
- **MultiVI** for RNA+ATAC. Works with both paired and mosaic data.
- **WNN** (Seurat) is an alternative but does not handle missing modalities.
- Always validate integration with known cell type markers in each modality.
