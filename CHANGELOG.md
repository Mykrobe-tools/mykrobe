# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


## [0.13.0]

### Added

- When using `--debug` flag, log species probe coverage/depth info and also each
  time a probe is rejected due to low coverage and/or depth.

- Added flag `--ncbi_names`. This is for soon to be updated tb probes for
  species calling, where the JSON file will also report alternative NCBI taxon
  names, as well as the default GTBD names.

- Added option `--dump_species_covgs`, to dump a JSON file of all the species
  probe coverage information. This includes the raw coverage info from mccortex,
  before it is aggregated into a single call that you see in the usual
  mykrobe output.

### Changed

- The tb species probes are going to be updated after the next release of
  mykrobe. Code changed to handle these new probes, specifically where a node in
  the taxon tree has no probes. If a child node is called as present from the
  reads, then push that call up to the parent node. Code still works as normal
  on the existing (soon to be old) panels.

- When reporting panel version in JSON file and stderr logging, only had
  the "main" version of the collection of panels, eg "20220705" for tb.
  Changed to also report the specific panel used, eg "20220705/202206".

## [0.12.2]

### Fixed

- Use all input files when `-1/-i/--seq` is passed multiple times [[#168][168]]
- When running on Staph, was throwing TypeError [[#170][170]]
- Windows build [[#172][172]]

## 0.12.1

### Fixed

-  output from `--version` flag (0.12.0 was reporting 0.11.0)

## 0.12.0

### Added

- This version is required to use the latest panel which combines the 2021 WHO
  catalogue ([doi:10/h298](https://doi.org/10/h298)) with the previous default mykrobe
  catalogue, as constructed in [doi:10/h299](https://doi.org/10/h299).

### Changed

- Loading JSON files will now handled (gzip) compressed files.
- When an X mutation resolves to an existing variant in the panel, instead of
  raising an error, do not output the X mutation. This means can still use
  Walker-2015 panel.

## 0.11.0

### Added

- Support for including stop codons in the panel with the `*` character.
- Panel version is now included in JSON output [[#119][119]]
- Shorthand CLI versions for common/long options:
    - `-D`: `--min_proportion_expected_depth`
    - `-P`: `--custom_probe_set_path`
    - `-R`: `--custom_variant_to_resistance_json`
    - `-L`: `--custom_lineage_json`
    - `-O`: `--format`
- Singularity container automatically built and added to each new release.

### Changed

- Default kmer size (21) made consistent across all subcommands [[#141][141]]
- Nanopore preset `--ont` defaults to an expected error rate of 0.08 and ploidy
  'haploid' [[#134][134]]
- VCF interface from pyvcf to pyvcf3
- Automated tests/builds all (ie linux, windows, mac) moved to github actions.
  Turned off travis and appveyor.

### Fixed

- Improved error messaging when an X amino acid resolves to a mutation already present
  in the panel.
- Do not assume reference fasta file name is the same as the chromosome [[#140][140]]

## 0.10.0

### Added

- `-o` as another option for `--output`
- `-i` as another option for `-1,--seq`
- `-A` as another option for `--report_all_calls`
- `-e` as another option for `--expected_error_rate`
- new species `typhi`

### Changed

- If `--tmp` is not specified now, it is set to [`tempfile.mkdtemp()`][mkdtemp]. By
  default, this is a directory that is created within the `$TMPDIR`. The reason for this
  is that currently, when mykrobe cleans up temporary files, it doesn't remove the
  temporary directory ([#113][113]). I wasn't really comfortable with removing the
  temporary directory in case the user had set it to something like the current
  directory. At least with this default now the temporary directory won't clog up the
  current directory after mykrobe is finished and the OS should clean it up
  periodically.
- Ensure `--memory` and `--threads` are passed on to mccortex [[#114][114]].
- Positional arguments `species` and `sample` are now options `-S,species` and
  `-s,--sample` respectively.

### Fixed

- Bug where the `--ont` flag resulted in mykrobe predict crashing if the panel
  contained any presence/absence genes [[#123][123]].
- Bug reporting variants that are "amino acid to any other amino acid" changes,
  eg `katG_S315X`. If gene is on reverse strand, then the wrong new amino
  acid change was being reported [[#127][127]].

[113]: https://github.com/Mykrobe-tools/mykrobe/issues/113

[114]: https://github.com/Mykrobe-tools/mykrobe/issues/114

[119]: https://github.com/Mykrobe-tools/mykrobe/issues/119

[123]: https://github.com/Mykrobe-tools/mykrobe/issues/123

[127]: https://github.com/Mykrobe-tools/mykrobe/issues/127

[134]: https://github.com/Mykrobe-tools/mykrobe/issues/134

[140]: https://github.com/Mykrobe-tools/mykrobe/issues/140

[141]: https://github.com/Mykrobe-tools/mykrobe/issues/141

[142]: https://github.com/Mykrobe-tools/mykrobe/issues/142

[168]: https://github.com/Mykrobe-tools/mykrobe/issues/168

[170]: https://github.com/Mykrobe-tools/mykrobe/issues/170

[172]: https://github.com/Mykrobe-tools/mykrobe/issues/172

[Unreleased]: https://github.com/Mykrobe-tools/mykrobe/compare/v0.13.0...HEAD

[0.12.2]: https://github.com/Mykrobe-tools/mykrobe/compare/v0.12.1...v0.12.2

[0.13.0]: https://github.com/Mykrobe-tools/mykrobe/compare/v0.12.2...v0.13.0

[mkdtemp]: https://docs.python.org/3.6/library/tempfile.html#tempfile.mkdtemp

