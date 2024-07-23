#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_INIT_WF         } from "./subworkflows/pipeline_init.nf"
include { CELLRANGER_REF_WF        } from "./subworkflows/cellranger_ref"
include { CELLRANGER_COUNT_AGGR_WF } from "./subworkflows/cellranger_count_aggr"
include { CLUSTER_SCANPY_WF        } from "./subworkflows/cluster_scanpy"
include { CLUSTER_SEURAT_WF        } from "./subworkflows/cluster_seurat"
include { ANNOTATE_SCANPY_WF       } from "./subworkflows/annotate_scanpy"
include { ANNOTATE_SEURAT_WF       } from "./subworkflows/annotate_seurat"
include { SCVELO_WF                } from "./subworkflows/scvelo"



// Workflow

workflow {

    main:

        // Create output directory if it does not exist
        new File ( params.outdir ).mkdirs()

        // Create channel to collect software versions
        ch_versions = Channel.empty()


        if ( params.input || params.cell_type_anno ) {
            //
            // SUBWORKFLOW: Generate cellranger genome index
            //
            CELLRANGER_REF_WF ( 
                params.species,
                params.genome,
                params.genome_ucsc,
                params.species_latin,
                params.ensembl_version,
                params.outdir,
                params.docker_enabled,
                params.cellranger_module,
            )
            ch_versions = ch_versions.mix(CELLRANGER_REF_WF.out.versions)
        }


        if ( params.input ) {

            println "Path to input sample sheet: $params.input (see nextflow.config file)"

            //
            // SUBWORKFLOW: Process input sample sheet
            //
            PIPELINE_INIT_WF (
                params.input,
            )

            //
            // SUBWORKFLOW: Align to genome with cellranger
            //
            CELLRANGER_COUNT_AGGR_WF (
                CELLRANGER_REF_WF.out.cellranger_index,
                params.input,
                PIPELINE_INIT_WF.out.input,
                params.project,
                params.outdir,
                params.docker_enabled,
                params.cellranger_module,
            )
            ch_versions = ch_versions.mix(CELLRANGER_COUNT_AGGR_WF.out.versions)
        } else {
            println "File specified by 'input' parameter not found: Skipping cellranger alignment."
        }


        if ( !params.input && !params.cellranger_out_dir ) {
            println "Error: File not found - parameter 'cellranger_out_dir' needs to be set if parameter 'input' is not set."
        } else if ( !params.input && params.cellranger_out_dir ) {
            cellranger_out_dir = params.cellranger_out_dir

            if ( params.metadata ) {
                cellranger_metadata = params.metadata
            } else {
                println "Error: File not found - parameter 'metadata' needs to be set if parameter 'input' is not set."
            }

        } else {
            cellranger_out_dir  = CELLRANGER_COUNT_AGGR_WF.out.cellranger_aggr_bc_matrix
            cellranger_metadata = CELLRANGER_COUNT_AGGR_WF.out.cellranger_metadata
        }


        //
        // SUBWORKFLOW: Cluster cells with Scanpy or Seurat
        //
        if ( params.scRNA_analysis == "scanpy" ) {
            CLUSTER_SCANPY_WF (
                cellranger_out_dir,
                cellranger_metadata,
                params.project,
                params.min_genes,
                params.max_genes,
                params.max_perc_mt,
                params.min_cells,
                params.n_pcs,
                params.harmony_var,
                params.leiden_res,
                params.outdir,
                params.docker_enabled,
                params.python_module,
            )
            ch_versions = ch_versions.mix(CLUSTER_SCANPY_WF.out.versions)
        } else if ( params.scRNA_analysis == "seurat" ) {
            CLUSTER_SEURAT_WF (
                cellranger_out_dir,
                cellranger_metadata,
                params.project,
                params.min_genes,
                params.max_genes,
                params.max_perc_mt,
                params.min_cells,
                params.n_pcs,
                params.harmony_var,
                params.leiden_res,
                params.outdir,
                params.docker_enabled,
                params.r_module,
            )
            ch_versions = ch_versions.mix(CLUSTER_SEURAT_WF.out.versions)
        } else {
            println "'scRNA_analysis' parameter needs to be: 'scanpy' or 'seurat'"
        }


        //
        // SUBWORKFLOW: Annotate scRNA-seq object
        //
        if ( params.cell_type_anno ) {
            if ( params.scRNA_analysis == "scanpy" ) {
                ANNOTATE_SCANPY_WF (
                    CLUSTER_SCANPY_WF.out.scrnaseq_object,
                    params.cell_type_anno,
                    params.project,
                    params.leiden_res,
                    params.outdir,
                    params.docker_enabled,
                    params.python_module,
                )
                ch_versions = ch_versions.mix(ANNOTATE_SCANPY_WF.out.versions)
                scrnaseq_object = ANNOTATE_SCANPY_WF.out.scrnaseq_object
            } else if ( params.scRNA_analysis == "seurat" ) {
                ANNOTATE_SEURAT_WF (
                    CLUSTER_SEURAT_WF.out.scrnaseq_object,
                    params.cell_type_anno,
                    params.project,
                    params.leiden_res,
                    params.outdir,
                    params.docker_enabled,
                    params.r_module,
                )
                ch_versions = ch_versions.mix(ANNOTATE_SEURAT_WF.out.versions)
                scrnaseq_object = ANNOTATE_SEURAT_WF.out.scrnaseq_object
            } else {
                println "'scRNA_analysis' parameter needs to be: 'scanpy' or 'seurat'"
            }
        } else {
            println "File specified with 'cell_type_anno' parameter not found: Skipping scRNA-seq cell type annotation"
        }


        //
        // SUBWORKFLOW: Perform mRNA velocity analysis
        //
        if ( params.cell_type_anno ) {
            if ( params.input ) {
                velocyto_input = CELLRANGER_COUNT_AGGR_WF.out.cellranger_count_out
            } else {
                if ( params.scRNA_velocity_file ) {
                    Channel
                        .fromPath( params.scRNA_velocity_file )
                        .splitCsv( header: true, sep: "\t" )
                        .map { line ->
                            tuple( line.sample_id, file( line.cellranger_count_dir, type: "dir", checkIfExists: true ) )
                        }
                        .set { velocyto_input }
                } else {
                    println "Both 'input' and 'scRNA_velocity_file' parameters not set: Skipping mRNA velocity analysis"
                }
            }
            SCVELO_WF (
                CELLRANGER_REF_WF.out.cellranger_index,
                velocyto_input,
                cellranger_metadata,
                scrnaseq_object,
                params.project,
                params.outdir,
                params.docker_enabled,
                params.samtools_module,
                params.python_module,
                params.r_module,
            )
            ch_versions = ch_versions.mix(SCVELO_WF.out.versions)
        } else {
            println "File specified with 'cell_type_anno' parameter not found: Skipping mRNA velocity analysis"
        }
}
