#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * Download genome files
 */
process GENOME_FILES_PR {
    debug false
    publishDir "${outdir}/genomes", pattern: "", mode: "copy", overwrite: false, saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "michaelweinberger/ubuntu-22.04:v1" : "" }
    
    input:
    val  ( species )
    val  ( genome )
    val  ( genome_ucsc )
    val  ( species_latin )
    val  ( ensembl_version )
    path ( outdir )
    val  ( docker_enabled )
    
    output:
    path ( "refdata-*"          ), emit: cellranger_index, optional: true
    path ( "${genome}.fa"       ), emit: genome_fasta, optional: true
    path ( "${genome}.gtf"      ), emit: genome_gtf, optional: true
    path ( "${genome}_rmsk.txt" ), emit: genome_repeat_mask, optional: true

    script:
    """
    [ ! -d "${outdir}/genomes" ] && mkdir -p "${outdir}/genomes"

    1_genome_files.sh -s "${species}" -g "${genome}" -u "${genome_ucsc}" -l "${species_latin}" -e "${ensembl_version}" -o "${outdir}"
    """
}



/*
 * Prepare cellranger genome index
 */
process CELLRANGER_REF_PR {
    debug false
    publishDir "${outdir}/genomes", pattern: "", mode: "copy", overwrite: false, saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "litd/docker-cellranger:v7.2.0" : "" }
    
    input:
    val  ( genome )
    path ( genome_fasta )
    path ( genome_gtf )
    path ( outdir )
    val  ( docker_enabled )
    val  ( cellranger_module )
    
    output:
    path ( "refdata-*"    ), emit: cellranger_index
    path ( "versions.txt" ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$cellranger_module"
    fi

    [ ! -d "${outdir}/genomes" ] && mkdir -p "${outdir}/genomes"

    # prepare cellranger reference data
    cellranger mkgtf \
        "${genome_gtf}" \
        "${genome}.filtered.gtf" \
        --attribute=gene_biotype:protein_coding

    cellranger mkref \
        --genome="refdata-cellranger-${genome}" \
        --fasta="${genome_fasta}" \
        --genes="${genome}.filtered.gtf" \
        --output-dir="\$PWD"

    cat <<-END_VERSIONS > versions.txt
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}



// Workflow

workflow CELLRANGER_REF_WF {

    take:
        species
        genome
        genome_ucsc
        species_latin
        ensembl_version
        outdir
        docker_enabled
        cellranger_module
	
    main:
        ch_versions = Channel.empty()

        GENOME_FILES_PR (
            species,
            genome,
            genome_ucsc,
            species_latin,
            ensembl_version,
            outdir,
            docker_enabled,
        )

        if ( species == "human" || species == "mouse" ) {
            cellranger_index = GENOME_FILES_PR.out.cellranger_index
        } else {
            CELLRANGER_REF_PR ( 
                genome,
                GENOME_FILES_PR.out.genome_fasta,
                GENOME_FILES_PR.out.genome_gtf,
                outdir,
                docker_enabled,
                cellranger_module,
            )
            ch_versions      = ch_versions.mix(CELLRANGER_REF_PR.out.versions)
            cellranger_index = CELLRANGER_REF_PR.out.cellranger_index
        }

    emit:
        versions           = ch_versions
        cellranger_index   = cellranger_index
        //genome_repeat_mask = GENOME_FILES_PR.out.genome_repeat_mask
}
