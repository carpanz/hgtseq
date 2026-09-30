process PARSE_INTEGRATION_SITES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_parsed_integration_sites.txt"), emit: integration_sites
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // For each unmapped read whose mate is mapped, report the read name (QNAME),
    // the chromosome (RNAME) and the position of the mapped mate (PNEXT):
    // these are the candidate integration sites in the host genome
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        view \\
        ${args} \\
        --threads ${task.cpus - 1} \\
        ${bam} \\
        | cut -f 1,3,8 > ${prefix}_parsed_integration_sites.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_parsed_integration_sites.txt
    """
}
