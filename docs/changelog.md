# Changelog

## Unreleased

- Added `--control-presence` filter with `any`, `present`, and `none` modes to retain cases based on whether ENCODE metadata records a control relationship for the experiment.
- Added pure post-processing function `filter_control_presence` in `encodefetch.postprocess`.
- Automatically retained matched control rows for all surviving cases while dropping orphaned controls.

## 0.5.0

- Added accession-file parsing for `--accessions`.
- Added metadata-only mode as the preferred way to skip downloads.
- Added exporter control strategies: `all`, `pool`, deterministic `first`, and metadata-ranked `best`.
- Improved file-level control mapping through `controlled_by_files`.
- Added nf-core and Snakemake exporters for ChIP-seq, ATAC-seq, and RNA-seq assay families.

## 0.1.0

- Initial release.
- Experiment-first search with control expansion.
- Multi-threaded metadata fetching and downloads.
- Paired FASTQ collapsing.
- nf-core and Snakemake samplesheet exporters.
- Rich progress bars for experiments and downloads.
