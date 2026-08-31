# Changelog

## Unreleased

- Added `--summary / --no-summary` option (default: enabled) that prints a compact KPI panel to the terminal after a run and writes two report files alongside the manifest:
  - `summary.html`: a self-contained, dependency-free report (no external JS or fonts, so it renders the same on an offline cluster) covering the selection criteria used, case/control experiment and file counts, download size, release date range, replication status for cases and controls, and a metadata breakdown split into "Per experiment" (Organism, Biosample, Target, Lab, Award, Platform, Donor Sex/Age/Life Stage/Ethnicity, Perturbed, Classification, Control Type - deduplicated by experiment, not inflated by replicate file counts) and "Per file" (File Format, Run Type, Output Type, Assembly, File Status - genuinely file-level properties).
  - `summary.json`: the same statistics in machine-readable form, including the CLI filters that produced the report, for automated pipelines and reproducibility reporting.
- Added unique experiment/file deduplication (including paired FASTQ R1+R2), estimated download size, and per-experiment biological/technical replication statistics (case and control tracked separately, since ENCODE only expects 2+ biological replicates for the case).
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
