#!/usr/bin/env python3



import importlib.metadata
import platform

import anndata as ad
import pandas as pd
import yaml

adata = ad.read_h5ad("${h5ad}")
target_col = "${target_col}"

# Priority-ordered candidate columns; the first present in obs wins
source_cols = "${source_cols ?: ''}".split()
present = [col for col in source_cols if col in adata.obs.columns]

if not present:
    raise ValueError(
        f"None of the requested cell type columns {source_cols} were found in obs. "
        f"Available columns: {list(adata.obs.columns)}"
    )

source_col = present[0]
adata.obs[target_col] = adata.obs[source_col].astype(str).values

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
