## CellyTypist-Nexflow  ![image](https://github.com/balqees-mansour/CellyTypist-Nexflow/assets/87857777/16194b6a-b1a3-4968-b45e-8ef88dd22c16)

celltypist.nf is a Nextflow script that defines a pipeline for  single-cell RNA-sequencing data Annotation using the CellTypist tool. The pipeline consists of several processes, and this script includes the QUALITY process and ANNOTATION process.

*** The script starts by defining the input parameters for the pipeline: *** 
params.query_input: The path to the input query data file in the H5AD format.
params.reference_file: The path to the reference data file in the H5AD format.
params.outdir: The output directory for the pipeline results.



the workflow block sets up the input channel for the QUALITY process by creating a channel from the params.query_input file path. It then executes the QUALITY process with this input channel, Then execute ANNOTATION process with the final OutPut channel which contains (Anndata  cell types annotated with plots).

