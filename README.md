# A Nextflow pipeline for scRNA-seq analysis
---

Use this Nextflow pipeline to analyse scRNA-seq data:
- Genome alignment via [Cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct)
- Clustering, marker identification, cell type annotation and visualisation via [Scanpy](https://scanpy.readthedocs.io/en/stable/)<sup>1</sup>  or [Seurat](https://satijalab.org/seurat/)<sup>2</sup>
- Doublet removal with [DoubletDetection](https://github.com/JonathanShor/DoubletDetection?tab=readme-ov-file)<sup>3</sup> (as part of Scanpy analysis) or [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)<sup>4</sup> (as part of Seurat analysis)
- Sample integration via [Harmony](https://github.com/immunogenomics/harmony)<sup>5</sup>
- mRNA velocity analysis via [velocyto](http://velocyto.org/velocyto.py/tutorial/cli.html#run10x-run-on-10x-chromium-samples)<sup>6</sup> and [scvelo](https://scvelo.readthedocs.io/en/stable/)<sup>7</sup>



## Use

1. Ensure Nextflow is installed
As this pipeline uses [Nextflow](https://www.nextflow.io/docs/latest/index.html)<sup>8</sup> to manage analysis workflow, compute resources and software, Nextflow needs to be installed before running the pipeline. 

2.  Clone the Github repository like 
```
    git clone https://github.com/michaelweinberger/scRNA-seq-Nextflow.git
```

3.  Change into the directory and start the pipeline via `nextflow run` like
```
    nextflow run main.nf \
    --input [path/to/input_file] \
    --clustering_mode ["scanpy"|"seurat"] \
    -profile HPC_docker,mouse \
    -resume
```

In this example, the "HPC_no_docker" profile directs Nextflow to run the pipeline in a high-performance computing environment using the slurm scheduler, and the "mouse" profile sets some genome parameters used to generate cellranger genome index files. The profile options are outlined in more detail below.

The `--input` flag gives the file path to a ".txt" or ".csv" input sample sheet, the `--clustering_mode` flag indicates if Scanpy (Python) or Seurat (R) should be used for cell clustering. These parameters can also be set in the nextflow.config file.

When setting the `-resume` flag, the pipeline will resume from a previous run.


### Usage scenarios

- Starting from fastq files\
If the `input` parameter is set in the nextflow.config file or on the command line, the workflow will start from fastq files. See below under "Parameters" for more information about `input` specifications.

The workflow first generates Cellranger genome index files. Input fastq files are then aligned to the reference genome via the `cellranger count` command and the output from individual samples is combined via `cellranger aggr`. This is followed by doublet removal, clustering, data integration and cluster marker identification. The clustering mode as well as specific parameters related to clustering and data integration can be set in the nextflow.config file (see below).

If no cell type annotation file is provided via the `cell_type_anno` parameter, the workflow stops at this point. You can inspect the output marker ".csv" or ".xlsx" files and the UMAP plots, and generate a cell type annotation ".csv" file as detailed below under "Parameters". Then add the file path in the nextflow.config file or add the `--cell_type_anno` flag on the command line. Run the pipeline again with the `-resume` flag.

With the `cell_type_anno` parameter set, the workflow will annotate cell clusters and perform mRNA velocity analysis.
\

- Starting from Cellranger outputs\
To start the workflow from previously generated Cellranger output files instead of from fastq files, set the `input` parameter in the nextflow.config file to an empty string. Instead, specify the `cellranger_out_dir` and `metadata` parameters. These specify a file path to a directory containing cellranger outputs, and a tab delimited file with cell barcodes and metadata, respectively. Please see below under "Parameters" for more details.

As with starting from fastq files, cell type annotation and mRNA velocity analysis will only be run if the `cell_type_anno` parameter is set in the nextflow.config file or if the `--cell_type_anno` flag is used on the command line.

Additionally, the `scRNA_velocity_file` parameter must be set for mRNA velocity analysis to be performed. This specifies the path to a file containing file paths of directories with Cellranger output BAM files, see below under "Parameters" for more details.


### Parameters

All parameters can be set on the command line with `--parameter_name` or in the `params` scope of the nextflow.config file located in the main pipeline directory.

- `project`   The name of the analysis project, defaults to "nf_scRNAseq_analysis".

- `outdir`    The name of the directory to save pipeline outputs to, defaults to "/out" in the main pipeline directory.\
            Within this directory, outputs from individual parts of the pipeline are written to different subdirectories:\
            - "genomes" for Cellranger genome index files\
            - "cellranger" for Cellranger mapping outputs\
            - "scanpy" for Scanpy clustering outputs\
            - "seurat" for Seurat clustering outputs\
            - "velocyto" for velocyto outputs (".loom" files containing counts of spliced and unspliced reads)\
            - "scvelo" for Scvelo mRNA velocity outputs\


Parameters specific to starting from fastq files:

- `input`     The file path to an input ".txt" (tab delimited) or ".csv" file containing sample information.\
            - The first column of the sample sheet must be named "sample_id" and contain sample-specific identifiers that also are prefixes in the corresponding fastq file names.\
            For example: Put "sample_x" as sample ID if your fastq files are named "sample_x_S2_L001_I1_001.fastq.gz", "sample_x_S2_L001_R1_001.fastq.gz", "sample_x_S2_L001_R2_001.fastq.gz" etc.\
            - The second column of the sample sheet must be named "fastq_dir" and contain the file paths to directories with fastq files to be analysed.\
            Fastq files of multiple samples may be located within the same directory.

            Set `cellranger_out_dir` and `metadata` instead of `input` to start the pipeline from previously computed Cellranger mapping results instead of fastq files. Make sure in this case that `input` is empty to prevent genome alignment from being run. Optionally set the `scRNA_velocity_file` parameter to run mRNA velocity analysis when starting from Cellranger outputs.


Parameters specific to starting from Cellranger output files:

- `cellranger_out_dir`    The file path to a directory containing cellranger output "barcodes.tsv.gz", "features.tsv.gz" and "matrix.mtx.gz" files,\ typically a "/outs/count/filtered_feature_bc_matrix" directory.

- `metadata`      The file path to a tab delimited file containing\ 
                - a column named "barcode" of cell barcodes,\
                - a column named "sample_id" of sample identifiers,\
                - optional metadata columns\
                Doublet detection is performed using the "sample_id" column.

- `scRNA_velocity_file`   Optional: The file path to a tab delimited file containing\
                        - a column named "sample_id" of sample identifiers,\
                        - a column named "cellranger_count_dir" of file paths to directories containing  Cellranger count BAM files.\ 
                        Each directory must be the direct parent directory of "/outs/possorted_genome_bam.bam", for example a `cellranger count` output directory.
\
                        Note: If the `input` parameter is not set, mRNA velocity analysis is only run if the `scRNA_velocity_file` parameter is set.


Parameters related to scRNA-seq clustering:

- `clustering_mode`   Clustering analysis mode, can be one of: "scanpy" or "seurat"\
- `min_genes`         Minimum number of genes expressed for a cell to be kept in the dataset, defaults to 200\
- `max_genes`         Maximum number of genes expressed for a cell to be kept in the dataset, defaults to 2500\
- `max_perc_mt`       Maximum percentage of mitochondrial reads for a cell to be kept in the dataset, defaults to 5\
- `min_cells`         Minimum number of cells for a gene to be expressed in to be kept in the dataset, defaults to 3\
- `n_pcs`             Number of principal components to be computed, defaults to 30\
- `harmony_var`       Name of the metadata column to use for Harmony data integration, defaults to ""\
- `leiden_res`        Resolution of cell clustering, defaults to 0.4\


Parameters related to scRNA-seq annotation and mRNA velocity analysis:

- `cell_type_anno`    The file path to a ".csv" file containing\ 
                    - a column named "cluster" of cell cluster numbers\
                    - a column named "cell_type" of cell type annotations\
                    - an optional column named "order" of integers indicating the order in which cell types should appear in UMAP plot legends\
\
                    Note: scRNA-seq annotation and mRNA velocity analysis are only run if the `cell_type_anno` parameter is set.


### Profiles

Multiple parameters can be bundled into profiles. These are defined in the `profiles` scope in the nextflow.config file and can be invoked on the command line via the `-profile` flag.\
Additional executor or genome profiles can be added in the nextflow.config file.

Pre-defined executor profiles are:

- `HPC_no_docker`   Use for pipeline execution on a high performance cluster without using Docker.\
\
                Pre-defined process options:\
                `executor = "slurm"`
                `queue    = "long"`
                `time     = "3days"`
                `memory   = "100 GB"`
                `cpus     = 30`
\
                Pre-defined parameters:\
                `docker_enabled    = false`                 -> Indicates whether software dependencies should be run out of Docker containers\
                `cellranger_module = "cellranger/7.2.0"`    -> Name of the Cellranger software module to load if Docker is disabled\
                `samtools_module   = "samtools/1.17"`       -> Name of the Samtools software module to load if Docker is disabled\
                `python_module     = "python-cbrg"`         -> Name of the Python software module to load if Docker is disabled\
                `r_module          = "R-cbrg"`              -> Name of the R software module to load if Docker is disabled\
\
                Pre-defined Docker options:\
                `docker.enabled = false`

- `HPC_docker`   Use for pipeline execution on a high performance cluster that allows the use of Docker.\
\
                Pre-defined process options:\
                `executor = "slurm"`
                `queue    = "long"`
                `time     = "3days"`
                `memory   = "100 GB"`
                `cpus     = 30`
\
                Pre-defined parameters:\
                `docker_enabled    = true`\   
                `cellranger_module = "none"`\
                `samtools_module   = "none"`\
                `python_module     = "none"`\
                `r_module          = "none"`\
\                
                Pre-defined Docker options:\
                `docker.enabled = true`\
                `docker.runOptions = '-u $(id -u):$(id -g)'`\

- `standard`     Use for local pipeline execution.\
\
                Pre-defined parameters:\
                `docker_enabled    = true`\   
                `cellranger_module = "none"`\
                `samtools_module   = "none"`\
                `python_module     = "none"`\
                `r_module          = "none"`\ 
\
                Pre-defined Docker options:\
                `docker.enabled = true`\
                `docker.runOptions = '-u $(id -u):$(id -g)'`\


Pre-defined genome profiles are:

- `human`\
\
                Pre-defined parameters:\
                `species         = "human"`\
                `species_latin   = "homo_sapiens"`  -> Species name used in genome browser files\
                `genome          = "GRCh38"`        -> Genome name used in ensembl genome browser\
                `genome_ucsc     = "hg38"`          -> Genome name used in UCSC genome browser\
                `ensembl_version = "110"`           -> Ensembl genome browser release to use for file downloads\-

- `mouse`\
\
                Pre-defined parameters:\
                `species         = "mouse"`\
                `species_latin   = "mus_musculus"`\
                `genome          = "GRCm39"`\
                `genome_ucsc     = "mm39"`\
                `ensembl_version = "110"`\

- `zebrafish`\
\
                Pre-defined parameters:\
                `species         = "zebrafish"`\
                `species_latin   = "danio_rerio"`\
                `genome          = "GRCz11"`\
                `genome_ucsc     = "danRer11"`\
                `ensembl_version = "110"`\



## References
1.	Wolf, F.A., Angerer, P., and Theis, F.J. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biol 19, 15. 10.1186/s13059-017-1382-0.
2.	Hao, Y., Hao, S., Andersen-Nissen, E., Mauck, W.M., 3rd, Zheng, S., Butler, A., Lee, M.J., Wilk, A.J., Darby, C., Zager, M., et al. (2021). Integrated analysis of multimodal single-cell data. Cell 184, 3573-3587 e3529. 10.1016/j.cell.2021.04.048.
3.	Gayoso, Adam, Shor, Jonathan, Carr, Ambrose J., Sharma, Roshan, Pe'er, Dana (2020, December 18). DoubletDetection (Version v3.0). Zenodo. http://doi.org/10.5281/zenodo.2678041
3.	McGinnis, C.S., Murrow, L.M., and Gartner, Z.J. (2019). DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. Cell Syst 8, 329-337 e324. 10.1016/j.cels.2019.03.003.
4.	Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., Baglaenko, Y., Brenner, M., Loh, P.R., and Raychaudhuri, S. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16, 1289-1296. 10.1038/s41592-019-0619-0.
5.	La Manno, G., Soldatov, R., Zeisel, A., Braun, E., Hochgerner, H., Petukhov, V., Lidschreiber, K., Kastriti, M.E., Lonnerberg, P., Furlan, A., et al. (2018). RNA velocity of single cells. Nature 560, 494-498. 10.1038/s41586-018-0414-6.
6.	Bergen, V., Lange, M., Peidli, S., Wolf, F.A., and Theis, F.J. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling. Nat Biotechnol 38, 1408-1414. 10.1038/s41587-020-0591-3.
7.  P. Di Tommaso, et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319



