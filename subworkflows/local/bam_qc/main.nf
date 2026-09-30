//
// Quality control of sorted and indexed BAM files
//

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { QUALIMAP_BAMQC    } from '../../../modules/nf-core/qualimap/bamqc/main'
include { BAMTOOLS_STATS    } from '../../../modules/nf-core/bamtools/stats/main'

workflow BAM_QC {

    take:
    ch_bam_bai  // channel: [ val(meta), path(bam), path(bai) ]
    gff         // path:    annotation used by Qualimap to restrict the analysis (optional, [] when missing)

    main:
    def ch_bam = ch_bam_bai.map { meta, bam, _bai -> [ meta, bam ] }

    SAMTOOLS_STATS ( ch_bam_bai, [ [], [], [] ] )

    SAMTOOLS_FLAGSTAT ( ch_bam_bai )

    SAMTOOLS_IDXSTATS ( ch_bam_bai )

    QUALIMAP_BAMQC ( ch_bam, gff )

    BAMTOOLS_STATS ( ch_bam )

    emit:
    stats    = SAMTOOLS_STATS.out.stats        // channel: [ val(meta), path(stats) ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat  // channel: [ val(meta), path(flagstat) ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats  // channel: [ val(meta), path(idxstats) ]
    qualimap = QUALIMAP_BAMQC.out.results      // channel: [ val(meta), path(results) ]
    bamstats = BAMTOOLS_STATS.out.stats        // channel: [ val(meta), path(stats) ]
}
