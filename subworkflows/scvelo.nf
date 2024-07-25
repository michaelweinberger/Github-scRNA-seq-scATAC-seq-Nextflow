#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * count spliced and unspliced reads in scRNA-seq BAM files with velocyto
 */
process VELOCYTO_PR {
    debug false
    publishDir "${outdir}/velocyto", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "michaelweinberger/python-3.11-velocyto:v1" : "" }
    //container "mparikhbroad/velocyto@sha256:52610aa2fb04ebe577613e149ba54421e537de3742204e9a52905432941b64f1"

    input:
    path  ( transcriptome_ref )
    tuple val( sample_id ), path( cellranger_count_outdir )
    path  ( outdir )
    val   ( docker_enabled )
    val   ( samtools_module )
    val   ( python_module )
    
    output:
    path ( "*.loom"       ), emit: loom_file // emitted for loom merge
    path ( "versions.txt" ), emit: versions

    script:
    """
    if [ "${docker_enabled}" == "false" ] ; then
        module load "${samtools_module}"
        module load "${python_module}"
    fi

    [ ! -d "${outdir}/velocyto" ] && mkdir -p "${outdir}/velocyto"

    velocyto run10x "${cellranger_count_outdir}" \
                    "${transcriptome_ref}/genes/genes.gtf"

    tmp_var="${cellranger_count_outdir}"
    loom_prefix="\$(echo "\${tmp_var##*/}")"
    mv "${cellranger_count_outdir}/velocyto/\${loom_prefix}.loom" "./${sample_id}.loom"

    echo "${task.process}:" > versions.txt
        samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
        python --version >> versions.txt
    """
}



/*
 * merge loom files of individual samples
 */
process MERGE_LOOM_PR {
    debug false
    publishDir "${outdir}/velocyto", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "michaelweinberger/python-3.11.9-scvelo:v1" : "" }

    input:
    path ( loom_files )
    path ( metadata )
    path ( outdir )
    val  ( docker_enabled )
    val  ( python_module )
    
    output:
    path ( "merged.loom"  ), emit: merged_loom // emitted for scvelo
    path ( "versions.txt" ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$python_module"
    fi

    [ ! -d "${outdir}/velocyto" ] && mkdir -p "${outdir}/velocyto"

    6_velocity_loom_merge.py -l "$loom_files" \
                             -m "$metadata" \
                             -o "\$PWD"

    echo "${task.process}:" > versions.txt
        python --version >> versions.txt
        python -c "import scvelo; print(f'scvelo,{scvelo.__version__}')" >> versions.txt
        python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    """
}



/*
 * export count_matrix, cell and gene metadata, and PCA results from Seurat object
 */
process EXPORT_SEURAT_DATA_PR {
    debug false
    publishDir "${outdir}/seurat", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "michaelweinberger/seurat-doubletfinder-harmony:v1" : "" }

    input:
    path ( seurat_object )
    val  ( out_name )
    path ( outdir )
    val  ( docker_enabled )
    val  ( r_module )
    
    output:
    path ( "velocity_counts_${out_name}.mtx"   ), emit: count_matrix
    path ( "velocity_metadata_${out_name}.csv" ), emit: cell_meta_data
    path ( "velocity_genes_${out_name}.csv"    ), emit: gene_meta_data
    path ( "velocity_pca_${out_name}.csv"      ), emit: pca_data
    path ( "versions.txt"                      ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$r_module"
    fi

    [ ! -d "${outdir}/seurat" ] && mkdir -p "${outdir}/seurat"

    7_export_seurat_data_for_anndata.R --args in_file="$seurat_object" \
                                              out_dir="\$PWD" \
                                              project_name="$out_name"

    echo "${task.process}:" > versions.txt
        R --version | head -n 1 >> versions.txt
        R -e "library(Matrix); print(paste('Matrix', packageVersion('Matrix')))" | grep '[1]' | tail -n 1 >> versions.txt
        R -e "library(patchwork); print(paste('patchwork', packageVersion('patchwork')))" | grep '[1]' | tail -n 1 >> versions.txt
        R -e "library(Seurat); print(paste('Seurat', packageVersion('Seurat')))" | grep '[1]' | tail -n 1 >> versions.txt
    """
}



/*
 * construct Anndata object from scRNA-seq count_matrix, cell and gene metadata, and PCA results
 */
process CONSTRUCT_ANNDATA_PR {
    debug false
    publishDir "${outdir}/seurat", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "michaelweinberger/python-3.11.9-scvelo:v1" : "" }

    input:
    path ( count_matrix )
    path ( cell_meta_data )
    path ( gene_meta_data )
    path ( pca_data )
    val  ( out_name )
    path ( outdir )
    val  ( docker_enabled )
    val  ( python_module )
    
    output:
    path ( "${out_name}_from_matrix_metadata_pca.h5ad" ), emit: anndata_object // emitted for scvelo
    path ( "versions.txt"                              ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$python_module"
    fi

    [ ! -d "${outdir}/seurat" ] && mkdir -p "${outdir}/seurat"

    8_construct_scanpy_object.py -mtx "$count_matrix" \
                                 -md "$cell_meta_data" \
                                 -g "$gene_meta_data" \
                                 -p "$pca_data" \
                                 -o "\$PWD" \
                                 -n "$out_name"

    echo "${task.process}:" > versions.txt
        python --version >> versions.txt
        python -c "import scvelo; print(f'scvelo,{scvelo.__version__}')" >> versions.txt
        python -c "import scanpy; print(f'scanpy,{scanpy.__version__}')" >> versions.txt
        python -c "import anndata; print(f'anndata,{anndata.__version__}')" >> versions.txt
        python -c "import scipy; print(f'scipy,{scipy.__version__}')" >> versions.txt
        python -c "import numpy; print(f'numpy,{numpy.__version__}')" >> versions.txt
        python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    """
}



/*
 * perform scvelo analysis
 */
process SCVELO_PR {
    debug false
    publishDir "${outdir}/scvelo", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "michaelweinberger/python-3.11.9-scvelo:v1" : "" }

    input:
    path ( anndata_object )
    path ( loom_file )
    val  ( out_name )
    path ( outdir )
    val  ( docker_enabled )
    val  ( python_module )
    
    output:
    path ( "*${out_name}*.pdf"                ), emit: plots_pdf
    path ( "*${out_name}*.png"                ), emit: plots_png
    path ( "${out_name}_dynamical_genes.csv"  ), emit: dynamical_genes_csv
    path ( "${out_name}_scvelo.h5ad"          ), emit: scvelo_object
    path ( "versions.txt"                     ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$python_module"
    fi

    [ ! -d "${outdir}/scvelo" ] && mkdir -p "${outdir}/scvelo"

    9_scvelo.py -i "$anndata_object" \
                -l "$loom_file" \
                -o "\$PWD" \
                -n "$out_name"

    echo "${task.process}:" > versions.txt
        python --version >> versions.txt
        python -c "import scvelo; print(f'scvelo,{scvelo.__version__}')" >> versions.txt
        python -c "import scanpy; print(f'scanpy,{scanpy.__version__}')" >> versions.txt
        python -c "import numpy; print(f'numpy,{numpy.__version__}')" >> versions.txt
        python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    """
}



// Workflow

workflow SCVELO_WF {

    take:
        cellranger_index
        velocyto_input
        cellranger_metadata
        scrnaseq_object
        out_name
        outdir
        docker_enabled
        samtools_module
        python_module
        r_module

    main:
        ch_versions = Channel.empty()

        // Generate loom file across samples containing count data of spliced and unspliced reads
        VELOCYTO_PR ( 
            cellranger_index,
            velocyto_input,
            outdir,
            docker_enabled,
            samtools_module,
            python_module,
        )
        ch_versions = ch_versions.mix(VELOCYTO_PR.out.versions)

        // Pass generated loom files simultaneously into merging process
        ch_loom_files = VELOCYTO_PR.out.loom_file.collect()

        MERGE_LOOM_PR (
            ch_loom_files,
            cellranger_metadata,
            outdir,
            docker_enabled,
            python_module,
        )
        ch_versions = ch_versions.mix(MERGE_LOOM_PR.out.versions)

        // Construct or declare input anndata object
        if ( params.clustering_mode == "seurat" ) {
            EXPORT_SEURAT_DATA_PR (
                scrnaseq_object,
                out_name,
                outdir,
                docker_enabled,
                r_module,
            )
            ch_versions = ch_versions.mix(EXPORT_SEURAT_DATA_PR.out.versions)

            CONSTRUCT_ANNDATA_PR (
                EXPORT_SEURAT_DATA_PR.out.count_matrix,
                EXPORT_SEURAT_DATA_PR.out.cell_meta_data,
                EXPORT_SEURAT_DATA_PR.out.gene_meta_data,
                EXPORT_SEURAT_DATA_PR.out.pca_data,
                out_name,
                outdir,
                docker_enabled,
                python_module,
            )
            ch_versions    = ch_versions.mix(CONSTRUCT_ANNDATA_PR.out.versions)
            anndata_object = CONSTRUCT_ANNDATA_PR.out.anndata_object
        } else if ( params.scRNA_analysis == "scanpy" ) {
            anndata_object = scrnaseq_object
        } else {
            println "'scRNA_analysis' parameter needs to be: 'scanpy' or 'seurat'"
        }

        // Run Scvelo analysis
        SCVELO_PR (
            anndata_object,
            MERGE_LOOM_PR.out.merged_loom,
            out_name,
            outdir,
            docker_enabled,
            python_module,
        )
        ch_versions = ch_versions.mix(SCVELO_PR.out.versions)

    emit:
        versions      = ch_versions
        scvelo_object = SCVELO_PR.out.scvelo_object
}
