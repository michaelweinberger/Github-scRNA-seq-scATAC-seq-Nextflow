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
    publishDir (
        "${outdir}/cellranger", 
        pattern: "versions.txt", 
        mode: "copy", 
        saveAs: { filename -> "${filename}" }
    )
    publishDir (
        "${outdir}/cellranger", 
        pattern: "cellranger_count_dir", 
        mode: "copy", 
        saveAs: { fn -> "${sample_id}_count" }
    )

    container { ( "$docker_enabled" ) ? "litd/docker-cellranger:v7.2.0" : "" }

    input:
    path  ( transcriptome_ref )
    tuple val( sample_id ), path( fastq_dir )
    path  ( outdir )
    val   ( docker_enabled )
    val   ( cellranger_module )
    
    output:
    path ( "cellranger_count_info.tsv" ), emit: cellranger_count_info // emitted for velocyto
    path ( "cellranger_count_dir"      ), emit: cellranger_count_dir  // emitted for publishing
    path ( "versions.txt"              ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$cellranger_module"
    fi

    [ ! -d "${outdir}/cellranger" ] && mkdir -p "${outdir}/cellranger"

    # run cellranger count
    cellranger count --id="cellranger_count_dir" \
		     --fastqs="$fastq_dir" \
		     --transcriptome="$transcriptome_ref" \
		     --sample="$sample_id"

    # save sample ID and full output path to file
    printf "sample_id\tcellranger_dir\n" > cellranger_count_info.tsv

    cellranger_dir="\$(echo \${PWD}/cellranger_count_dir)"

    printf "${sample_id}\t\${cellranger_dir}\n" >> cellranger_count_info.tsv

    echo "${task.process}:" > versions.txt
        echo cellranger: "\$(cellranger --version 2>&1 | awk '{print \$(NF)}' )" | sed -e \$'s/^/\t/' >> versions.txt
    """
}



/*
 * generate cellranger aggr input .csv file
 */
process CELLRANGER_AGGR_INPUT_PR {
    debug false
    publishDir (
        "${outdir}/cellranger", 
        pattern: "", 
        mode: "copy", 
        saveAs: { filename -> "${filename}" }
    )

    input:
    path  ( sample_sheet )          // sample IDs & metadata
    path  ( cellranger_count_info ) // sample IDs & full filepaths cellranger count output dirs
    path  ( outdir )
    
    output:
    path  ( "cellranger_aggr_input.csv" ), emit: cellranger_aggr_input // emitted for cellranger aggr

    script:
    """
    [ ! -d "${outdir}/cellranger" ] && mkdir -p "${outdir}/cellranger"

    2_cellranger_aggr_input.sh -s "${sample_sheet}" -i "${cellranger_count_info}"
    """
}



/*
 * generate barcode metadata file
 */
process CELLRANGER_METADATA_PR {
    debug false
    publishDir (
        "${outdir}/cellranger", 
        pattern: "", 
        mode: "copy", 
        saveAs: { filename -> "${filename}" }
    )

    input:
    path ( cellranger_aggr_input_csv )
    path ( outdir )
    
    output:
    path ( "cellranger_aggr_cell_metadata.tsv" ), emit: metadata // emitted for scanpy,seurat and loom merge

    script:
    """
    [ ! -d "${outdir}/cellranger" ] && mkdir -p "${outdir}/cellranger"

    3_cellranger_metadata.sh -i "${cellranger_aggr_input_csv}"
    """
}



/*
 * run cellranger aggr
 */
process CELLRANGER_AGGR_PR {
    debug false
    publishDir (
        "${outdir}/cellranger", 
        pattern: "", 
        mode: "copy", 
        saveAs: { filename -> "${filename}" }
    )

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
    path ( "versions.txt"                                      ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$cellranger_module"
    fi

    [ ! -d "${outdir}/cellranger" ] && mkdir -p "${outdir}/cellranger"

    # Run cellranger aggr
    cellranger aggr --id="$out_name" \
	--csv="$cellranger_aggr_input_csv"

    echo "${task.process}:" > versions.txt
        echo cellranger: "\$(cellranger --version 2>&1 | awk '{print \$(NF)}' )" | sed -e \$'s/^/\t/' >> versions.txt
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
        ch_versions = ch_versions.mix(CELLRANGER_COUNT_PR.out.versions)

        // collect sample IDs and full paths of cellranger count output dirs
        CELLRANGER_COUNT_PR.out.cellranger_count_info
        .collectFile (
            name: "cellranger_count_info.tsv", 
            storeDir: "${outdir}/cellranger", 
            keepHeader: true
        )
        .set { cellranger_count_info }

        // generate tuple of sample IDs and full paths of cellranger count output dirs for velocyto
        cellranger_count_info
        .splitCsv( header: true , sep: "\t" )
        .map { line ->
            tuple( line.sample_id, line.cellranger_dir )
        }
        .set { cellranger_count_out }

        // generate cellranger aggr input csv file
        CELLRANGER_AGGR_INPUT_PR (
            sample_sheet, 
            cellranger_count_info,
            outdir,
        )

        // generate cell barcode level metadata file
        CELLRANGER_METADATA_PR (
            CELLRANGER_AGGR_INPUT_PR.out.cellranger_aggr_input,
            outdir,
        )

        // run cellranger aggr
        CELLRANGER_AGGR_PR ( 
            CELLRANGER_AGGR_INPUT_PR.out.cellranger_aggr_input,
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
