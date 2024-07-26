#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * annotate cell types in Scanpy object
 */
process ANNOTATE_SCANPY_PR {
    debug true
    publishDir "${outdir}/scanpy", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "michaelweinberger/python-3.11.9-scanpy:v1" : "" }

    input:
    path ( scrnaseq_object )
    path ( cell_type_annotation )
    val  ( out_name )
    val  ( leiden_res )
    path ( outdir )
    val  ( docker_enabled )
    val  ( python_module )
    
    output:
    path ( "*${out_name}*.pdf"                                    ), emit: plots_pdf
    path ( "*${out_name}*.png"                                    ), emit: plots_png
    path ( "${out_name}_obs.csv"                                  ), emit: cell_metadata_csv
    path ( "${out_name}_scRNAseq_no_doublets_annotated.h5ad"      ), emit: scrnaseq_object
    path ( "${out_name}_scRNAseq_no_doublets_annotated_scVI.h5ad" ), emit: scrnaseq_scvi_object
    path ( "versions.txt"                                         ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$python_module"
    fi

    [ ! -d "${outdir}/scanpy" ] && mkdir -p "${outdir}/scanpy"

    5_scanpy_plotting.py -i "$scrnaseq_object" \
                         -a "$cell_type_annotation" \
                         -o "\$PWD" \
                         -n "$out_name" \
                         -r "$leiden_res"

    echo "${task.process}:" > versions.txt
        python --version | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import scanpy; print(f'scanpy,{scanpy.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import anndata; print(f'anndata,{anndata.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import numpy; print(f'numpy,{numpy.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import pandas; print(f'pandas,{pandas.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import matplotlib; print(f'matplotlib,{matplotlib.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import seaborn; print(f'seaborn,{seaborn.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
    """
}



// Workflow

workflow ANNOTATE_SCANPY_WF {

    take:
        scrnaseq_object
        cell_type_annotation
        out_name
        leiden_res
        outdir
        docker_enabled
        python_module

    main:
        ch_versions = Channel.empty()

        ANNOTATE_SCANPY_PR ( 
            scrnaseq_object,
            cell_type_annotation,
            out_name,
            leiden_res,
            outdir,
            docker_enabled,
            python_module,
        )
        ch_versions = ch_versions.mix(ANNOTATE_SCANPY_PR.out.versions)

    emit:
        versions             = ch_versions
        scrnaseq_object      = ANNOTATE_SCANPY_PR.out.scrnaseq_object
        scrnaseq_scvi_object = ANNOTATE_SCANPY_PR.out.scrnaseq_scvi_object
}
