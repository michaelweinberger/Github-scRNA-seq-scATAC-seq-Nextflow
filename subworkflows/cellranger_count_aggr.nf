#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * run cellranger count
 */
process CELLRANGER_COUNT_PR {
    debug       false
    tag         "$sample_id"
    publishDir  "${outdir}/cellranger", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }
    publishDir  "${outdir}/cellranger", pattern: "cellranger_count_dir", mode: "copy", saveAs: "${sample_id}_count"

    container { ( "$docker_enabled" ) ? "litd/docker-cellranger:v7.2.0" : "" }

    input:
    path  ( transcriptome_ref )
    tuple val( sample_id ), path( fastq_dir )
    path  ( outdir )
    val   ( docker_enabled )
    val   ( cellranger_module )
    
    output:
    tuple val( sample_id ), path( cellranger_dir ), emit: cellranger_count_out // emitted for velocyto
    path ( "cellranger_count_dir"                ), emit: cellranger_count_dir // emitted for publishing
    path ( "versions.yml"                        ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$cellranger_module"
    fi

    [ ! -d "${outdir}/cellranger" ] && mkdir -p "${outdir}/cellranger"

    # get full output directory path
    cellranger_dir="\$(echo \${PWD}/cellranger_count_dir)"

    # run cellranger count
    cellranger count --id="cellranger_count_dir" \
		     --fastqs="$fastq_dir" \
		     --transcriptome="$transcriptome_ref" \
		     --sample="$sample_id"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}



/*
 * generate cellranger aggr input .csv file
 */
process CELLRANGER_AGGR_INPUT_PR {
    debug       false
    tag         "$sample_id"
    publishDir  "${outdir}/cellranger", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    input:
    path  ( sample_sheet )
    tuple val ( sample_id ), path ( cellranger_count_outdir )
    path  ( outdir )
    
    output:
    path  ( "cellranger_aggr_input.csv" ), emit: cellranger_aggr_input // emitted for cellranger aggr

    shell:
    '''
    [ ! -d "!{outdir}/cellranger" ] && mkdir -p "!{outdir}/cellranger"

    # create header for aggregator input file
    sed 1q "!{sample_sheet}" | awk -F\\t -vOFS=, '{$1=$2=""; print $0}' | sed 's/^,//' > cellranger_aggr_input.csv
    sed -i "s/^/sample_id,molecule_h5/" cellranger_aggr_input.csv

    # adjust sample id
    sample="$(grep "!{sample_id}" "!{sample_sheet}" | awk -F\\t '{print $1}' | tr , _)"

    # extract additional metadata columns and paste together with comma as separator
    metadata="$(grep "!{sample_id}" "!{sample_sheet}" | awk -F\\t -vOFS=, '{$1=$2=""; print $0}' | sed 's/^,//')"

    # create molecule_info .h5 file path
    h5_path="!{cellranger_count_outdir}/outs/molecule_info.h5"

    # append everything to .csv file
    echo "${sample},${h5_path}${metadata}" >> cellranger_aggr_input.csv
    '''
}



/*
 * generate barcode metadata file
 */
process CELLRANGER_METADATA_PR {
    debug       false
    publishDir  "${outdir}/cellranger", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    input:
    path ( cellranger_aggr_input_csv )
    path ( outdir )
    
    output:
    path ( "cellranger_aggr_cell_metadata.tsv" ), emit: metadata // emitted for scanpy,seurat and loom merge

    script:
    """
    [ ! -d "!{outdir}/cellranger" ] && mkdir -p "!{outdir}/cellranger"

    3_cellranger_metadata.sh -i "${cellranger_aggr_input_csv}"
    """
}



/*
 * run cellranger aggr
 */
process CELLRANGER_AGGR_PR {
    debug       false
    publishDir  "${outdir}/cellranger", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "litd/docker-cellranger:v7.2.0" : "" }

    input:
    path ( cellranger_aggr_input_csv )
    val  ( out_name )
    path ( outdir )
    val  ( docker_enabled )
    val  ( cellranger_module )
    
    output:
    path ( "${out_name}/outs/count/filtered_feature_bc_matrix" ), emit: cellranger_aggr_bc_matrix // emitted for scanpy,seurat
    path ( out_name                                            ), emit: cellranger_aggr_outdir    // emitted for publishing
    path ( "versions.yml"                                      ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$cellranger_module"
    fi

    [ ! -d "${outdir}/cellranger" ] && mkdir -p "${outdir}/cellranger"

    # Run cellranger aggr
    cellranger aggr --id="$out_name" \
	--csv="$cellranger_aggr_input_csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}



// Workflow

workflow CELLRANGER_COUNT_AGGR_WF {

    take:
        cellranger_index
        sample_sheet
        input                     // tuple sample ID + fastq dir
        cellranger_aggr_out_name
        outdir
        docker_enabled
        cellranger_module

    main:
        ch_versions = Channel.empty()

        // run cellranger count
        CELLRANGER_COUNT_PR ( 
            cellranger_index,
            input,
            outdir,
            docker_enabled,
            cellranger_module,
        )
        ch_versions          = ch_versions.mix(CELLRANGER_COUNT_PR.out.versions)
        cellranger_count_out = CELLRANGER_COUNT_PR.out.cellranger_count_out.flatten()

        // generate cellranger aggr input csv file
        CELLRANGER_AGGR_INPUT_PR (
            sample_sheet, 
            cellranger_count_out,
            outdir,
        )
        cellranger_aggr_input = CELLRANGER_AGGR_INPUT_PR.out.cellranger_aggr_input.collectFile (
            name: "cellranger_aggr_input.csv", 
            storeDir: "${outdir}/cellranger", 
            keepHeader: true
        )

        // generate cell barcode level metadata file
        CELLRANGER_METADATA_PR (
            cellranger_aggr_input,
            outdir,
        )

        // run cellranger aggr
        CELLRANGER_AGGR_PR ( 
            cellranger_aggr_input,
            cellranger_aggr_out_name,
            outdir,
            docker_enabled,
            cellranger_module,
        )
        ch_versions = ch_versions.mix(CELLRANGER_AGGR_PR.out.versions)

    emit:
        versions                   = ch_versions
        cellranger_count_out       = cellranger_count_out
        cellranger_metadata        = CELLRANGER_METADATA_PR.out.metadata
        cellranger_aggr_bc_matrix  = CELLRANGER_AGGR_PR.out.cellranger_aggr_bc_matrix
}
