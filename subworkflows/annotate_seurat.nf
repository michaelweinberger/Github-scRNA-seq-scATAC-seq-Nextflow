#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * annotate cell types in Seurat object
 */
process ANNOTATE_SEURAT_PR {
    debug true
    publishDir "${outdir}/seurat", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "michaelweinberger/r-base-4.4.2-seurat:v1" : "" }

    input:
    path (scrnaseq_object)
    path (cell_type_annotation)
    val  (out_name)
    val  (leiden_res)
    path (outdir)
    val  (docker_enabled)
    val  (r_module)
    
    output:
    path ( "*${out_name}*.pdf"                              ), emit: plots_pdf, optional: true
    path ( "*${out_name}*.png"                              ), emit: plots_png, optional: true
    path ( "${out_name}_metadata.csv"                       ), emit: cell_metadata_csv
    path ( "${out_name}_scRNAseq_no_doublets_annotated.rds" ), emit: scrnaseq_object
    path ( "versions.yml"                                   ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$r_module"
    fi

    [ ! -d "${outdir}/seurat" ] && mkdir -p "${outdir}/seurat"

    5_seurat_plotting.R --args in_file="$scrnaseq_object" \
                               annotation="$cell_type_annotation" \
                               out_dir="\$PWD" \
                               project_name="$out_name" \
                               leiden_res="$leiden_res"

    echo "${task.process}:" > versions.yml
    R --version | sed 's/^/r-base,/' >> versions.txt
    """
}



// Workflow

workflow ANNOTATE_SEURAT_WF {

    take:
        scrnaseq_object
        cell_type_annotation
        out_name
        leiden_res
        outdir
        docker_enabled
        r_module

    main:
        ch_versions = Channel.empty()

        ANNOTATE_SEURAT_PR ( 
            scrnaseq_object,
            cell_type_annotation,
            out_name,
            leiden_res,
            outdir,
            docker_enabled,
            r_module,
        )
        ch_versions = ch_versions.mix(ANNOTATE_SEURAT_PR.out.versions)

    emit:
        versions             = ch_versions
        scrnaseq_object      = ANNOTATE_SEURAT_PR.out.scrnaseq_object
}
