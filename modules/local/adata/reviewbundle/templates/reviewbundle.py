#!/usr/bin/env python3

import importlib.metadata
import platform
import re
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy.sparse import csc_matrix

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
cell_type_col = "${cell_type_col}"
group_cols = "${group_cols ?: ''}".split()

Path("annotations").mkdir(exist_ok=True)
Path("plots").mkdir(exist_ok=True)


def find_umap_key(adata):
    for key in adata.obsm.keys():
        if key.lower().endswith("umap"):
            return key
    return None


def slugify(value):
    return re.sub(r"[^0-9a-zA-Z]+", "_", value).strip("_").lower()


umap_key = find_umap_key(adata)

# Per-annotator label tables, so reviewers can inspect each tool side by side
annotation_cols = [
    col
    for col in adata.obs.columns
    if col.startswith("annotation:") and col.endswith(":per_cell")
]

for col in annotation_cols:
    tool = slugify(col.split(":")[1]) if len(col.split(":")) > 1 else "annotation"
    table = pd.DataFrame({col: adata.obs[col].astype(str)})
    table.to_csv(f"annotations/{prefix}_{tool}_{slugify(col)}.csv")

# UMAP coloured by the canonical cell type and each requested grouping column
plot_cols = []
if cell_type_col and cell_type_col in adata.obs.columns:
    plot_cols.append(cell_type_col)
plot_cols.extend([col for col in group_cols if col in adata.obs.columns and col not in plot_cols])

if umap_key is not None:
    for col in plot_cols:
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.embedding(
            adata,
            basis=umap_key,
            color=col,
            legend_loc="on data",
            show=False,
            ax=ax,
        )
        fig.tight_layout()
        fig.savefig(f"plots/{prefix}_umap_{slugify(col)}.png", dpi=150)
        plt.close(fig)

# CELLxGENE-compatible object: keep only embeddings, the canonical cell type and grouping columns
cellxgene = adata.copy()
for key in list(cellxgene.obsm.keys()):
    if not key.lower().endswith(("umap", "tsne")):
        del cellxgene.obsm[key]
keep_cols = [col for col in plot_cols]
cellxgene.obs = cellxgene.obs[keep_cols] if keep_cols else cellxgene.obs[[]]
cellxgene.uns = {}
cellxgene.layers = {}

for df in [cellxgene.obs, cellxgene.var]:
    for col in df.columns:
        if df[col].dtype == np.float64:
            df[col] = df[col].astype(np.float32)

cellxgene.X = csc_matrix(cellxgene.X).astype(np.float32)
sc.pp.log1p(cellxgene)
cellxgene.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
        "scanpy": sc.__version__,
        "matplotlib": matplotlib.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
