//
// Collate the Kraken2 classifications of all samples into Krona plots and an RMarkdown analysis report
//

include { GAWK as GAWK_SINGLE                                   } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_BOTH                                     } from '../../../modules/nf-core/gawk/main'
include { KRONA_KTIMPORTTAXONOMY as KRONA_KTIMPORTTAXONOMY_SINGLE } from '../../../modules/nf-core/krona/ktimporttaxonomy/main'
include { KRONA_KTIMPORTTAXONOMY as KRONA_KTIMPORTTAXONOMY_BOTH   } from '../../../modules/nf-core/krona/ktimporttaxonomy/main'
include { RANALYSIS                                             } from '../../../modules/local/ranalysis/main'

workflow REPORTING {

    take:
    ch_classified_single    // channel: [ val(meta), path(txt) ]
    ch_classified_both      // channel: [ val(meta), path(txt) ]
    ch_integration_sites    // channel: [ val(meta), path(txt) ]
    ch_kronadb              // channel: path(taxonomy.tab)
    rmarkdown_template      // path:    RMarkdown template of the analysis report
    taxonomy_id             // string:  NCBI taxonomy ID of the host
    istest                  // boolean: whether the pipeline runs on test data

    main:
    // Collate the classified reads of all samples, keeping only reads assigned to a taxon other than root
    GAWK_SINGLE (
        ch_classified_single.map { _meta, txt -> txt }.collect().map { txts -> [ [ id: 'group' ], txts, 'txt' ] },
        [],
        false
    )

    GAWK_BOTH (
        ch_classified_both.map { _meta, txt -> txt }.collect().map { txts -> [ [ id: 'group' ], txts, 'txt' ] },
        [],
        false
    )

    KRONA_KTIMPORTTAXONOMY_SINGLE ( GAWK_SINGLE.out.output, ch_kronadb )

    KRONA_KTIMPORTTAXONOMY_BOTH ( GAWK_BOTH.out.output, ch_kronadb )

    RANALYSIS (
        ch_classified_single.map { _meta, txt -> txt }.collect(),
        ch_classified_both.map { _meta, txt -> txt }.collect(),
        ch_integration_sites.map { _meta, txt -> txt }.collect(),
        ch_classified_single.map { meta, _txt -> meta.id }.collect(),
        rmarkdown_template,
        istest,
        taxonomy_id
    )

    emit:
    krona_single = KRONA_KTIMPORTTAXONOMY_SINGLE.out.html  // channel: [ val(meta), path(html) ]
    krona_both   = KRONA_KTIMPORTTAXONOMY_BOTH.out.html    // channel: [ val(meta), path(html) ]
    report       = RANALYSIS.out.report                    // channel: path(html)
}
