#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * run Seurat clustering
 */
process CLUSTER_SEURAT_PR {
    debug true
    publishDir "${outdir}/seurat", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "michaelweinberger/r-base-4.4.2-seurat:v1" : "" }

    input:
    path (cellranger_out_dir)
    path (metadata)
    val  (out_name)
    val  (min_genes)
    val  (max_genes)
    val  (max_perc_mt)
    val  (min_cells)
    val  (n_pcs)
    val  (harmony_var)
    val  (leiden_res)
    path (outdir)
    val  (docker_enabled)
    val  (r_module)
    
    output:
    path ( "*${out_name}*.pdf"                             ), emit: plots_pdf
    path ( "*${out_name}*.png"                             ), emit: plots_png
    path ( "${out_name}_doublets.csv"                      ), emit: doublets_csv
    path ( "${out_name}_markers.xlsx"                      ), emit: markers_csv
    path ( "${out_name}_metadata.csv"                      ), emit: cell_metadata_csv
    path ( "${out_name}_scRNAseq_analysed_no_doublets.rds" ), emit: scrnaseq_object
    path ( "versions.yml"                                  ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$r_module"
    fi

    [ ! -d "${outdir}/seurat" ] && mkdir -p "${outdir}/seurat"

    4_seurat.R --args data_dir="$cellranger_out_dir" \
                      metadata_file="$metadata" \
                      out_dir="\$PWD" \
                      project_name="$out_name" \
                      min_genes="$min_genes" \
                      max_genes="$max_genes" \
                      max_perc_mt="$max_perc_mt" \
                      min_cells="$min_cells" \
                      n_pcs="$n_pcs" \
                      harmony_var="$harmony_var" \
                      leiden_res="$leiden_res"

    echo "${task.process}:" > versions.yml
    R --version | sed 's/^/r-base,/' >> versions.txt
    """
}



// Workflow

workflow CLUSTER_SEURAT_WF {

    take:
        cellranger_out_dir
        metadata
        out_name
        min_genes
        max_genes
        max_perc_mt
        min_cells
        n_pcs
        harmony_var
        leiden_res
        outdir
        docker_enabled
        r_module

    main:
        ch_versions = Channel.empty()

        CLUSTER_SEURAT_PR ( 
            cellranger_out_dir,
            metadata,
            out_name,
            min_genes,
            max_genes,
            max_perc_mt,
            min_cells,
            n_pcs,
            harmony_var,
            leiden_res,
            outdir,
            docker_enabled,
            r_module,
        )
        ch_versions = ch_versions.mix(CLUSTER_SEURAT_PR.out.versions)

    emit:
        versions          = ch_versions
        scrnaseq_object   = CLUSTER_SEURAT_PR.out.scrnaseq_object
}
