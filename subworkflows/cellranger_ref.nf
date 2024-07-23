#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * Prepare cellranger genome index
 */
process CELLRANGER_REF_PR {
    debug false
    publishDir "${outdir}/genomes", pattern: "", mode: "copy", overwrite: false, saveAs: { filename -> "${filename}" }

    container { ( "$docker_enabled" ) ? "litd/docker-cellranger:v7.2.0" : "" }
    
    input:
    val  (species)
    val  (genome)
    val  (genome_ucsc)
    val  (species_latin)
    val  (ensembl_version)
    path (outdir)
    val  (docker_enabled)
    val  (cellranger_module)
    
    output:
    path ("refdata-*"          ), emit: cellranger_index
    path ("${genome}_rmsk.txt" ), emit: genome_repeat_mask, optional: true
    path ("versions.yml"       ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$cellranger_module"
    fi

    [ ! -d "${outdir}/genomes" ] && mkdir -p "${outdir}/genomes"

    2_genome_files_SH.sh -s "${species}" -g "${genome}" -u "${genome_ucsc}" -l "${species_latin}" -e "${ensembl_version}" -o "${outdir}"

    cat <<-END_VERSIONS > versions.yml
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

        CELLRANGER_REF_PR ( 
            species,
            genome,
            genome_ucsc,
            species_latin,
            ensembl_version,
            outdir,
            docker_enabled,
            cellranger_module,
        )
        ch_versions = ch_versions.mix(CELLRANGER_REF_PR.out.versions)

    emit:
        versions           = ch_versions
        cellranger_index   = CELLRANGER_REF_PR.out.cellranger_index
        //genome_repeat_mask = CELLRANGER_REF_PR.out.genome_repeat_mask
}
