Here is a short example of how to use the pipeline:

Download the entire directory to the computer you intend to run the pipeline either by clicking on the "Code" menu and downloading the zip file or from the command line with

`git clone https://github.com/Winshipe/omm12_CallConfidence_workflow.git`

Download the mobileOG database from <https://mobileogdb.flsi.cloud.vt.edu/> (beatrix 1.6 All) to `resources/mmseqs/`

and run the following commands

``` bash
conda activate workflow/envs/mmseqs2.yaml

cd resources/mmseqs #assuming you start from the main pipeline directory

mmseqs2 createdb mobileOG-db_beatrix-1.6.All.faa mobileOG

cd ../.. #to return to the main directory
```

I included a copy of the *A. muris* genome from NCBI (Accession CP065321.1) in resources and the supplied `config.yaml` file already has everything else filled out to run for *A. muris* with a mutant allele frequency of 0.5. All that's left is to run the pipeline with the following command

``` bash
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml--use-conda--cores 16
```

an example output is located in `OMM12_results/Acutalibacter_muris_report.html` and the output is also included in `OMM12_results/omm12_report.pdf`

If you're running snakemake on a cluster with a scheduler, you will likely need to add an executor plugin to allow snakemake to talk to the cluster's scheduler. <https://snakemake.github.io/snakemake-plugin-catalog/index.html> If you're unsure, talk to your cluster's support staff
