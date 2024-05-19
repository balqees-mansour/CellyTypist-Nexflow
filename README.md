
  # CellyTypist-Nexflow 
![image](https://github.com/balqees-mansour/CellyTypist-Nexflow/assets/87857777/16194b6a-b1a3-4968-b45e-8ef88dd22c16)

celltypist.nf is a Nextflow script that defines a pipeline for  single-cell RNA-sequencing data Annotation using the CellTypist tool. The pipeline consists of several processes, and this script includes the QUALITY process and ANNOTATION process.

*** The script starts by defining the input parameters for the pipeline: ***


params.query_input: The path to the input query data file in the H5AD format.
 
params.reference_file: The path to the reference data file in the H5AD format.

params.outdir: The output directory for the pipeline results.



the workflow block sets up the input channel for the QUALITY process by creating a channel from the params.query_input file path. It then executes the QUALITY process with this input channel and then executes the ANNOTATION process with the final output channel which contains (Anndata  cell types annotated with plots).

![image](https://github.com/balqees-mansour/CellyTypist-Nexflow/assets/87857777/399006c2-0ad9-43b7-9700-14c1b01e6c41)

![image](https://github.com/balqees-mansour/CellyTypist-Nexflow/assets/87857777/d9c4000b-d49f-4102-95d3-fd81f0a1b539)


