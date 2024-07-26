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

    container { ( "$docker_enabled" ) ? "michaelweinberger/seurat-doubletfinder-harmony:v1" : "" }

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
    path ( "versions.txt"                                   ), emit: versions

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

    echo "${task.process}:" > versions.txt
        R --version | head -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(dplyr); print(paste('dplyr', packageVersion('dplyr')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(ggplot2); print(paste('ggplot2', packageVersion('ggplot2')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(patchwork); print(paste('patchwork', packageVersion('patchwork')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(pheatmap); print(paste('pheatmap', packageVersion('pheatmap')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(Seurat); print(paste('Seurat', packageVersion('Seurat')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(viridis); print(paste('viridis', packageVersion('viridis')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
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
