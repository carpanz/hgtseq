process RANALYSIS {
    tag "analysis_report"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8309ae0799cf8625166500cc634dae1212101f93583ebc0a49898b9fcfc6c256/data'
        : 'community.wave.seqera.io/library/bioconductor-biovizbase_bioconductor-genomeinfodb_bioconductor-genomicranges_bioconductor-ggbio_pruned:baa766c046284ca8'}"

    input:
    path classified_reads_single, stageAs: 'classified_single/*'
    path classified_reads_both, stageAs: 'classified_both/*'
    path integration_sites
    val sampleids
    path markdownfile
    val istest
    val taxonomy_id

    output:
    path "analysis_report.html" , emit: report
    path "analysis_report.RData", emit: rdata
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('bioconductor-ggbio'), eval("Rscript -e \"cat(as.character(packageVersion('ggbio')))\""), emit: versions_ggbio, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samplestring = sampleids.join(',')
    """
    Rscript -e "rmarkdown::render(
        '${markdownfile}',
        params = list(
            sampleids = '${samplestring}',
            istest = '${istest}',
            taxonomy_id = '${taxonomy_id}'
        ),
        knit_root_dir = getwd(),
        output_dir = getwd()
    )"
    """

    stub:
    """
    touch analysis_report.html
    touch analysis_report.RData
    """
}
