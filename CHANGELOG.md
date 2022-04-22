# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
[Unreleased]: https://github.com/Mykrobe-tools/mykrobe/compare/v0.10.0...HEAD
[mkdtemp]: https://docs.python.org/3.6/library/tempfile.html#tempfile.mkdtemp

