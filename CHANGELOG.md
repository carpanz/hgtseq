# nf-core/hgtseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0dev - [unreleased]

### `Added`

- Samples provided as FastQ files and as BAM files can be mixed in the same samplesheet: the type of each sample is detected from the samplesheet columns
- The samplesheet is validated with [nf-schema](https://nextflow-io.github.io/nf-schema/), including the uniqueness of the sample names
- Tool citations and bibliography in the MultiQC methods description
- `meta.yml` documentation for the local subworkflows
- nf-test pipeline tests for the `test` and `test_bam` profiles

### `Changed`

- Sync `TEMPLATE` with nf-core/tools 4.1.0
- The pipeline code is written in the Nextflow [strict syntax](https://www.nextflow.io/docs/latest/strict-syntax.html), and requires Nextflow `>=25.10.4`
- The software versions are collected through the `versions` topic channel
- Parameters are only accessed in `main.nf`, and passed explicitly to the workflows and subworkflows
- The local subworkflows have been reorganised: `PREPARE_READS` and `READS_QC` are merged in `FASTQ_QC_TRIM_ALIGN`, while `SORTBAM` is replaced by the `samtools/sort` module, which now also indexes the alignments
- The local modules `SAMTOOLS_VIEW_SINGLE` and `GAWK` are replaced by the nf-core modules `samtools/view` (with the equivalent flag filter `-f 5 -F 264`) and `gawk`
- The local modules `SAMTOOLS_FASTQ` and `PARSEOUTPUTS` are renamed `SAMTOOLS_BAMTOFASTQ` and `PARSE_INTEGRATION_SITES`, and use the Seqera Containers samtools image
- The RMarkdown analysis report uses a [Seqera Containers](https://seqera.io/containers/) image built from its conda environment, instead of the custom `lescailab/r-ggbio-reporting` images, so that all the container engines and conda use the same software versions
- The RMarkdown analysis report has a conda environment, and is now generated also with the `conda`/`mamba` profiles (previously the whole reporting, including the Krona plots, was skipped)
- The analysis report loads `GenomeInfoDb` explicitly, as required by recent Bioconductor versions
- The Qualimap resource requirements are defined in `conf/base.config`
- The samplesheet columns are `sample,fastq_1,fastq_2` for FastQ files and `sample,bam` for BAM files

### `Fixed`

- The output directory `preprocess/exctracted_reads` is renamed `preprocess/extracted_reads`
- `SORTBAM` was called with a different number of inputs for FastQ and BAM samples
- The FastQ samplesheet parsing was broken for BAM input
- The local modules do not modify the `meta` map in place anymore
- `test_full` uses `resourceLimits` instead of the removed `max_cpus`/`max_memory`/`max_time` parameters

### `Dependencies`

| Dependency    | Old version | New version |
| ------------- | ----------- | ----------- |
| `bamtools`    | 2.5.1       | 2.5.2       |
| `bwa`         | 0.7.17      | 0.7.19      |
| `bwa-mem2`    | 2.2.1       | 2.3         |
| `fastqc`      | 0.11.9      | 0.12.1      |
| `gawk`        | 5.1.0       | 5.3.1       |
| `kraken2`     | 2.1.2       | 2.1.6       |
| `krona`       | 2.8         | 2.8.1       |
| `multiqc`     | 1.14        | 1.35        |
| `qualimap`    | 2.2.2d      | 2.3         |
| `r-base`      | 4.2.0       | 4.5.3       |
| `ggbio`       | 1.44.0      | 1.58.0      |
| `samtools`    | 1.17        | 1.24        |
| `trim-galore` | 0.6.7       | 2.3.0       |

### `Deprecated`

- The `--isbam` parameter is removed: the input type is detected from the samplesheet
- The `--hook_url` parameter is removed, following the nf-core template

## [1.1.0](https://github.com/nf-core/hgtseq/releases/tag/1.1.0) - Beary Rose

### `Added`

- [#31](https://github.com/nf-core/hgtseq/pull/31) - Fixed issue where _single_unmapped_ reads also include _both_unmapped_ reads, by creating a local module with two steps samtools flag filtering

### `Fixed`

- [#31](https://github.com/nf-core/hgtseq/pull/31) - Sync `TEMPLATE` with `tools 2.8` and all nf-core/modules updated

### `Dependencies`

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `samtools` | 1.15.1      | 1.17        |
| `multiqc`  | 1.13        | 1.14        |

### `Deprecated`

## [1.0.0](https://github.com/nf-core/hgtseq/releases/tag/1.0.0) - Dalmatian Daffodil

Initial release of nf-core/hgtseq, created with the [nf-core](https://nf-co.re/) template.
