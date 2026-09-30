//
// Quality control, adapter trimming and alignment of paired-end FastQ files to the host genome
//

include { FASTQC as FASTQC_RAW     } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../../../modules/nf-core/fastqc/main'
include { TRIMGALORE               } from '../../../modules/nf-core/trimgalore/main'
include { BWA_INDEX                } from '../../../modules/nf-core/bwa/index/main'
include { BWA_MEM                  } from '../../../modules/nf-core/bwa/mem/main'
include { BWAMEM2_INDEX            } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM              } from '../../../modules/nf-core/bwamem2/mem/main'

workflow FASTQ_QC_TRIM_ALIGN {

    take:
    ch_reads  // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    ch_fasta  // channel: [ val(meta), path(fasta) ]
    index     // string:  path to a pre-built aligner index (optional, built from the FASTA when missing)
    aligner   // string:  'bwa-mem' or 'bwa-mem2'

    main:
    FASTQC_RAW ( ch_reads )

    TRIMGALORE ( ch_reads )

    FASTQC_TRIMMED ( TRIMGALORE.out.reads )

    //
    // Use the pre-built index when available, otherwise index the reference.
    // The reference is only indexed when there is at least one FastQ sample to align.
    //
    def ch_fasta_to_index = ch_reads
        .first()
        .combine(ch_fasta)
        .map { _meta, _reads, meta_fasta, fasta -> [ meta_fasta, fasta ] }

    def ch_index = channel.empty()
    if (index) {
        ch_index = channel.value([ [ id: 'index' ], file(index, checkIfExists: true) ])
    } else if (aligner == 'bwa-mem2') {
        BWAMEM2_INDEX ( ch_fasta_to_index )
        ch_index = BWAMEM2_INDEX.out.index.first()
    } else {
        BWA_INDEX ( ch_fasta_to_index )
        ch_index = BWA_INDEX.out.index.first()
    }

    def ch_bam = channel.empty()
    if (aligner == 'bwa-mem2') {
        BWAMEM2_MEM ( TRIMGALORE.out.reads, ch_index, ch_fasta, false )
        ch_bam = BWAMEM2_MEM.out.bam
    } else {
        BWA_MEM ( TRIMGALORE.out.reads, ch_index, ch_fasta, false )
        ch_bam = BWA_MEM.out.bam
    }

    emit:
    fastqc_raw     = FASTQC_RAW.out.zip      // channel: [ val(meta), [ zip ] ]
    fastqc_trimmed = FASTQC_TRIMMED.out.zip  // channel: [ val(meta), [ zip ] ]
    trimmed_reads  = TRIMGALORE.out.reads    // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    bam            = ch_bam                  // channel: [ val(meta), path(bam) ]
}
