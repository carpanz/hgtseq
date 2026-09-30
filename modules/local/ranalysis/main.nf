process RANALYSIS {
    tag "analysis_report"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'library://lescailab/hgtseq/r-ggbio-reporting:sha256.eb829b05cf12e8d827813a6afb6e38592aac6568f685a6519f5ed7dd20125cb3'
        : 'ghcr.io/lescailab/r-ggbio-reporting:1.0.0'}"

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
