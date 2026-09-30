//
// Extract the unmapped reads from the alignments and classify them taxonomically.
// Two categories of read pairs are analysed separately:
//   - single: the read is unmapped but its mate is mapped (SAM flag 5, excluding 8 and 256);
//             the position of the mapped mate is used to infer a candidate integration site
//   - both:   the read and its mate are both unmapped (SAM flag 13, excluding 256)
//

include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_SINGLE             } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_BOTH               } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_BAMTOFASTQ as SAMTOOLS_BAMTOFASTQ_SINGLE } from '../../../modules/local/samtools/bamtofastq/main'
include { SAMTOOLS_BAMTOFASTQ as SAMTOOLS_BAMTOFASTQ_BOTH   } from '../../../modules/local/samtools/bamtofastq/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_SINGLE                 } from '../../../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_BOTH                   } from '../../../modules/nf-core/kraken2/kraken2/main'
include { PARSE_INTEGRATION_SITES                           } from '../../../modules/local/parse_integration_sites/main'

workflow CLASSIFY_UNMAPPED {

    take:
    ch_bam_bai   // channel: [ val(meta), path(bam), path(bai) ]
    ch_krakendb  // channel: path(kraken2_db)

    main:
    SAMTOOLS_VIEW_SINGLE ( ch_bam_bai, [ [], [], [] ], [ [], [] ], [ [], [] ], '' )

    SAMTOOLS_VIEW_BOTH ( ch_bam_bai, [ [], [], [] ], [ [], [] ], [ [], [] ], '' )

    PARSE_INTEGRATION_SITES ( SAMTOOLS_VIEW_SINGLE.out.bam )

    SAMTOOLS_BAMTOFASTQ_SINGLE ( SAMTOOLS_VIEW_SINGLE.out.bam )

    SAMTOOLS_BAMTOFASTQ_BOTH ( SAMTOOLS_VIEW_BOTH.out.bam )

    // The extracted reads are written to a single FastQ file and are classified as single-end reads
    KRAKEN2_SINGLE (
        SAMTOOLS_BAMTOFASTQ_SINGLE.out.fastq.map { meta, fastq -> [ meta + [ single_end: true ], fastq ] },
        ch_krakendb,
        false,
        true
    )

    KRAKEN2_BOTH (
        SAMTOOLS_BAMTOFASTQ_BOTH.out.fastq.map { meta, fastq -> [ meta + [ single_end: true ], fastq ] },
        ch_krakendb,
        false,
        true
    )

    emit:
    classified_single      = KRAKEN2_SINGLE.out.classified_reads_assignment  // channel: [ val(meta), path(txt) ]
    classified_both        = KRAKEN2_BOTH.out.classified_reads_assignment    // channel: [ val(meta), path(txt) ]
    report_single          = KRAKEN2_SINGLE.out.report                       // channel: [ val(meta), path(txt) ]
    report_both            = KRAKEN2_BOTH.out.report                         // channel: [ val(meta), path(txt) ]
    candidate_integrations = PARSE_INTEGRATION_SITES.out.integration_sites   // channel: [ val(meta), path(txt) ]
}
