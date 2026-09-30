/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQ_QC_TRIM_ALIGN        } from '../subworkflows/local/fastq_qc_trim_align'
include { BAM_QC                     } from '../subworkflows/local/bam_qc'
include { CLASSIFY_UNMAPPED          } from '../subworkflows/local/classify_unmapped'
include { REPORTING                  } from '../subworkflows/local/reporting'
include { UNTAR as UNTAR_KRAKEN2_DB  } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_KRONA_DB    } from '../modules/nf-core/untar/main'
include { SAMTOOLS_SORT              } from '../modules/nf-core/samtools/sort/main'
include { MULTIQC                    } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap           } from 'plugin/nf-schema'
include { paramsSummaryMultiqc       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText     } from '../subworkflows/local/utils_nfcore_hgtseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HGTSEQ {

    take:
    ch_samplesheet              // channel: [ val(meta), [ fastq_1, fastq_2 ] ] or [ val(meta), [ bam ] ]
    fasta                       // string:  path to the host reference genome FASTA
    gff                         // string:  path to the host annotation GFF/GTF, used by Qualimap (optional)
    aligner_index               // string:  path to a pre-built BWA-MEM or BWA-MEM2 index (optional)
    aligner                     // string:  'bwa-mem' or 'bwa-mem2'
    krakendb                    // string:  path to the Kraken2 database, as a directory or a .tar.gz archive
    kronadb                     // string:  path to the Krona taxonomy .tab file, or a .tar.gz archive containing it
    taxonomy_id                 // string:  NCBI taxonomy ID of the host organism
    istest                      // boolean: whether the pipeline runs on test data
    multiqc_runkraken           // boolean: whether to include the Kraken2 reports in MultiQC
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Prepare the reference genome and the classification databases
    //
    def ch_fasta = fasta
        ? channel.value([ [ id: file(fasta).baseName ], file(fasta, checkIfExists: true) ])
        : channel.empty()
    def gff_file = gff ? file(gff, checkIfExists: true) : []

    def ch_krakendb = channel.empty()
    if (krakendb.toLowerCase().endsWith('.tar.gz')) {
        UNTAR_KRAKEN2_DB ( channel.value([ [ id: 'kraken2_db' ], file(krakendb, checkIfExists: true) ]) )
        ch_krakendb = UNTAR_KRAKEN2_DB.out.untar.map { _meta, db -> db }
    } else {
        ch_krakendb = channel.value(file(krakendb, checkIfExists: true))
    }

    def ch_kronadb = channel.empty()
    if (kronadb.toLowerCase().endsWith('.tar.gz')) {
        UNTAR_KRONA_DB ( channel.value([ [ id: 'krona_db' ], file(kronadb, checkIfExists: true) ]) )
        ch_kronadb = UNTAR_KRONA_DB.out.untar.map { _meta, db -> db }
    } else {
        ch_kronadb = channel.value(file(kronadb, checkIfExists: true))
    }

    //
    // Samples can be provided either as paired-end FastQ files or as aligned BAM files
    //
    def ch_input = ch_samplesheet
        .branch { meta, _files ->
            fastq: meta.data_type == 'fastq'
            bam:   meta.data_type == 'bam'
        }

    //
    // SUBWORKFLOW: QC, trimming and alignment of FastQ files to the host genome
    //
    FASTQ_QC_TRIM_ALIGN (
        ch_input.fastq,
        ch_fasta,
        aligner_index,
        aligner
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC_TRIM_ALIGN.out.fastqc_raw.map { _meta, zip -> zip })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC_TRIM_ALIGN.out.fastqc_trimmed.map { _meta, zip -> zip })

    //
    // MODULE: Sort and index both the newly aligned reads and the input BAM files
    //
    SAMTOOLS_SORT (
        FASTQ_QC_TRIM_ALIGN.out.bam.mix(ch_input.bam),
        [ [], [], [] ],
        'bai'
    )
    def ch_bam_bai = SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_SORT.out.index, failOnDuplicate: true, failOnMismatch: true)

    //
    // SUBWORKFLOW: Alignment QC
    //
    BAM_QC (
        ch_bam_bai,
        gff_file
    )
    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.stats.map { _meta, stats -> stats })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.flagstat.map { _meta, flagstat -> flagstat })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.idxstats.map { _meta, idxstats -> idxstats })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.qualimap.map { _meta, qualimap -> qualimap })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.bamstats.map { _meta, bamstats -> bamstats })

    //
    // SUBWORKFLOW: Extraction and taxonomic classification of the unmapped reads
    //
    CLASSIFY_UNMAPPED (
        ch_bam_bai,
        ch_krakendb
    )
    // With small test databases too few reads are classified to generate a meaningful Kraken2 report
    if (multiqc_runkraken) {
        ch_multiqc_files = ch_multiqc_files.mix(CLASSIFY_UNMAPPED.out.report_single.map { _meta, report -> report })
        ch_multiqc_files = ch_multiqc_files.mix(CLASSIFY_UNMAPPED.out.report_both.map { _meta, report -> report })
    }

    //
    // SUBWORKFLOW: Krona plots and RMarkdown analysis report
    //
    REPORTING (
        CLASSIFY_UNMAPPED.out.classified_single,
        CLASSIFY_UNMAPPED.out.classified_both,
        CLASSIFY_UNMAPPED.out.candidate_integrations,
        ch_kronadb,
        file("${projectDir}/assets/analysis_report.Rmd", checkIfExists: true),
        taxonomy_id,
        istest
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'hgtseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'hgtseq'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                                                     // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
