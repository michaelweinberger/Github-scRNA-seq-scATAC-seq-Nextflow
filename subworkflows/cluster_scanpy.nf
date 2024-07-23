#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * run Scanpy clustering
 */
process CLUSTER_SCANPY_PR {
    debug true
    publishDir "${outdir}/scanpy", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "michaelweinberger/python-3.11.9-scanpy:v1" : "" }

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
    val  (python_module)
    
    output:
    path ( "*${out_name}*.pdf"                              ), emit: plots_pdf
    path ( "*${out_name}*.png"                              ), emit: plots_png
    path ( "${out_name}_doublets.csv"                       ), emit: doublets_csv
    path ( "${out_name}_markers.csv"                        ), emit: markers_csv
    path ( "${out_name}_obs.csv"                            ), emit: cell_metadata_csv
    path ( "${out_name}_scRNAseq_analysed_no_doublets.h5ad" ), emit: scrnaseq_object
    path ( "versions.yml"                                   ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$python_module"
    fi

    [ ! -d "${outdir}/scanpy" ] && mkdir -p "${outdir}/scanpy"

    4_scanpy.py -i "$cellranger_out_dir" \
                -m "$metadata" \
                -o "\$PWD" \
                -n "$out_name" \
                -mig "$min_genes" \
                -mag "$max_genes" \
                -mam "$max_perc_mt" \
                -mic "$min_cells" \
                -npc "$n_pcs" \
                -hv "$harmony_var" \
                -r "$leiden_res"

    echo "${task.process}:" > versions.yml
    python --version | sed 's/^/python,/' >> versions.txt
    """
}



// Workflow

workflow CLUSTER_SCANPY_WF {

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
        python_module

    main:
        ch_versions = Channel.empty()

        CLUSTER_SCANPY_PR ( 
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
            python_module,
        )
        ch_versions = ch_versions.mix(CLUSTER_SCANPY_PR.out.versions)

    emit:
        versions          = ch_versions
        scrnaseq_object   = CLUSTER_SCANPY_PR.out.scrnaseq_object
}
