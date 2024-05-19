params.query_input = "/workspace/gitpod/cell/query_integrated.h5ad"
params.reference_file = "/workspace/gitpod/cell/adata_reference.h5ad"
params.outdir = "results"

println "query: $params.query_input"
println "reference: $params.reference_file"

log.info """
    CELLTYPIST - N F   P I P E L I N E
    ===================================
    query: ${params.query_input}
    Reference: ${params.reference_file}
    outdir: ${params.outdir}
    """.stripIndent(true)

process QUALITY {
    // Input file
    input:
    path query

    // Output directory
    output:
    path "quality_results"

    // Script to execute
    script:
    """
export MPLCONFIGDIR="/workspace/gitpod/cell"
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import scanpy as sc
import celltypist
from celltypist import models
import numpy as np
import pandas as pd

anndata = sc.read("$query")


sc.settings.verbosity = 3

# Print the header
sc.logging.print_header()

# Set figure parameters
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Make variable names unique
anndata.var_names_make_unique()

# Plot the highest expressed genes
sc.pl.highest_expr_genes(anndata, n_top=20, save="highest_expr_genes.pdf")

# Filter cells and genes
sc.pp.filter_cells(anndata, min_genes=200)
sc.pp.filter_genes(anndata, min_cells=3)

# Annotate mitochondrial genes
anndata.var['mt'] = anndata.var_names.str.startswith('MT-')

# Calculate QC metrics
sc.pp.calculate_qc_metrics(anndata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Plot violin plots
sc.pl.violin(anndata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
           jitter=0.4, multi_panel=True, save="violin_plots.pdf")

anndata = anndata[anndata.obs.n_genes_by_counts > 2500, :]
anndata = anndata[anndata.obs.pct_counts_mt < 5, :]

sc.pp.normalize_total(anndata, target_sum=1e4)

sc.pp.log1p(anndata)

sc.pp.highly_variable_genes(anndata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(anndata)

anndata.raw = anndata
anndata = anndata[:, anndata.var.highly_variable]

sc.pp.regress_out(anndata, ['total_counts', 'pct_counts_mt'])

sc.pp.scale(anndata, max_value=10)

sc.tl.pca(anndata, svd_solver='arpack')

sc.pl.pca(anndata, color='CST3')
"""
}

workflow {
    // Create a channel from the input file path
    quality_ch = Channel.fromPath(params.query_input)
    
    // Execute the QUALITY process with the input channel
    QUALITY(quality_ch)
}


