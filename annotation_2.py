#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 14:01:13 2023

@author: balqees
"""
''' Using CellTypist for cell type classification ''' 
# pip install celltypist
import scanpy as sc
import celltypist
from celltypist import models
import numpy as np
import pandas as pd


anndata = sc.read("processed_adata.h5ad")

## Assign cell type labels using a CellTypist built-in model

# Enabling `force_update = True` will overwrite existing (old) models.
models.download_models(force_update = True)

models.models_path

''' Get an overview of the models and what they represent.'''
models.models_description()

# Indeed, the `model` argument defaults to `Immune_All_Low.pkl`.
model = models.Model.load(model = 'Immune_All_Low.pkl')

# Not run; predict cell identities using this loaded model.
#predictions = celltypist.annotate(adata_2000, model = model, majority_voting = True)
# Alternatively, just specify the model name (recommended as this ensures the model is intact every time it is loaded).
predictions = celltypist.annotate(anndata, model = 'Immune_All_Low.pkl', majority_voting = True, mode = 'best match')

# Get an `AnnData` with predicted labels embedded into the cell metadata columns.
anndata = predictions.to_adata()
anndata.obs

# If the UMAP or any cell embeddings are already available in the `AnnData`, skip this command.
sc.tl.umap(anndata)



sc.pl.umap(anndata, color="HES4")
sc.pl.umap(anndata, color = ['predicted_labels'], size=60, legend_fontsize = "xx-small", save = "labelBM_low_umap.png")
sc.pl.umap(anndata, color = ['majority_voting'], size=60, legend_fontsize = "xx-small", save = "majorityBN_low_umap.png")
