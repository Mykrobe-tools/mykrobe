<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Mykrobe](#mykrobe)
  - [Requirements](#requirements)
  - [Installation](#installation)
    - [pipenv (virtual environment)](#pipenv-virtual-environment)
    - [Bioconda](#bioconda)
    - [Containers](#containers)
    - [pip](#pip)
  - [Usage](#usage)
    - [AMR prediction (Mykrobe predictor)](#amr-prediction-mykrobe-predictor)
    - [Examples](#examples)
    - [Output](#output)
  - [Citing](#citing)
  - [Genotyping a pre-built probe set](#genotyping-a-pre-built-probe-set)
    - [Examples](#examples-1)
    - [Make a custom probe set (for use with `mykrobe genotype`)](#make-a-custom-probe-set-for-use-with-mykrobe-genotype)
    - [Add variants to the database (for background/context)](#add-variants-to-the-database-for-backgroundcontext)
    - [Make probes and dump-probes](#make-probes-and-dump-probes)
    - [Examples](#examples-2)
      - [1 Simple case - building a probe without using backgrounds](#1-simple-case---building-a-probe-without-using-backgrounds)
      - [2. 'Dumping' the Variant database](#2-dumping-the-variant-database)
      - [3. Building a custom probe set](#3-building-a-custom-probe-set)
        - [Build a variant probe set defined based on reference co-ordinates (1-based)](#build-a-variant-probe-set-defined-based-on-reference-co-ordinates-1-based)
        - [Build a variant probe set defined based on gene co-ordinates (1-based)](#build-a-variant-probe-set-defined-based-on-gene-co-ordinates-1-based)
- [Citation](#citation)
- [Tests](#tests)
- [Release](#release)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# Mykrobe

<http://www.mykrobe.com>


## Installation

### From source

Requirements:

-   C++ compiler (to compile mccortex during install)
-   Python > 3.4
-   [mongodb](https://www.mongodb.com/) > 3.0 (optional, not needed to run `mykrobe predict`)

Either git clone:
```sh
git clone https://github.com/Mykrobe-tools/mykrobe.git mykrobe
```
or download the [latest release](https://github.com/Mykrobe-tools/mykrobe/releases/latest)
source code archive and extract it.

Then install with:
```sh
cd mykrobe
pip3 install .
```

If you get installation relating to compiling `mccortex`, then it may be helpful
to debug by trying to compile `mccortex` first. The `mykrobe` install looks for
the compiled binary file `./mccortex/bin/mccortex31`, and if it finds it then
it does not try to recompile, and simply copies the file.
This means you can make this binary before running `pip3`, like this:
```sh
git clone --recursive -b geno_kmer_count https://github.com/Mykrobe-tools/mccortex mccortex
cd mccortex
make
```
and once the `make` runs successfully, `pip3 install .` can be run from the `mykrobe/`
directory.

### Bioconda
[![conda badge](https://anaconda.org/bioconda/mykrobe/badges/installer/conda.svg)][bioconda]

Before attempting to install with bioconda, please ensure you have your channels set up
as specified in the [documentation][bioconda-channels]. If you don't, you may run into
issues with an older version of `mykrobe` being installed.

To install [this package][bioconda] with `conda` in a dedicated environment:

```sh
conda install -c bioconda mykrobe
```

[bioconda]: https://anaconda.org/bioconda/mykrobe
[bioconda-channels]: https://bioconda.github.io/user/install.html#set-up-channels

### Containers

[biocontainers]: https://biocontainers.pro/#/
[container]: https://quay.io/repository/biocontainers/mykrobe?tab=tags

[Biocontainers][biocontainers] maintain images for all bioconda recipes. The container
and all tags for `mykrobe` can be found [here][container]. To use a specific version,
just select your required version/tag and use the URI as follows.

```sh
tag="0.7.0--py37h2666aa9_0"
uri="quay.io/biocontainers/mykrobe:${tag}"

# using Singularity
singularity exec docker://"$uri" mykrobe --help
# using docker
docker pull "$uri"
```

## Test on example data

You can test that `mykrobe predict` runs as expected by downloading and running
on a small toy set of reads, as follows.

```sh
wget -O test_reads.fq.gz https://ndownloader.figshare.com/files/21059229
mykrobe predict SAMPLE tb --output out.json --format json --seq test_reads.fq.gz
```

The test reads are simulated, and perfectly match the reference, except for
the isoniazid resistant associated variant inhA I21T. You should see a section
like this in the output file `out.json`:

```json
"Isoniazid": {
    "predict": "R",
    "called_by": {
        "inhA_I21T-ATC1674262ACT": {
    ... etc
```


## Usage

```
    mykrobe --help
    usage: mykrobe [-h] [--version]
                         {predict,panels,variants,vars,genotype} ...

    optional arguments:
      -h, --help            show this help message and exit
      --version             mykrobe version

    [sub-commands]:
      {predict,panels,variants,vars,genotype}
        predict             predict the sample's drug susceptibility
        panels              A description of the AMR panels available within
                            Mykrobe predict
        variants (vars)     build variant probes
        genotype            genotype a sample using a probe set
```


### AMR prediction (Mykrobe predictor)

```
mykrobe predict --help
usage: mykrobe predict [-h] [-k kmer] [--tmp TMP] [--keep_tmp]
                         [--skeleton_dir SKELETON_DIR]
                         [--mccortex31_path MCCORTEX31_PATH] [-t THREADS]
                         [-m MEMORY] [--expected_depth EXPECTED_DEPTH]
                         [-1 seq [seq ...]] [-c ctx] [-f] [--ont]
                         [--guess_sequence_method] [--ignore_minor_calls]
                         [--ignore_filtered IGNORE_FILTERED]
                         [--model model] [--ploidy ploidy]
                         [--filters FILTERS [FILTERS ...]]
                         [--report_all_calls]
                         [--expected_error_rate EXPECTED_ERROR_RATE]
                         [--min_variant_conf MIN_VARIANT_CONF]
                         [--min_gene_conf MIN_GENE_CONF]
                         [--min_proportion_expected_depth MIN_PROPORTION_EXPECTED_DEPTH]
                         [--min_gene_percent_covg_threshold MIN_GENE_PERCENT_COVG_THRESHOLD]
                         [--output OUTPUT] [-q] [--panel panel]
                         [--custom_probe_set_path custom_probe_set_path]
                         [--custom_variant_to_resistance_json custom_variant_to_resistance_json]
                         [--min_depth min_depth]
                         [--conf_percent_cutoff conf_percent_cutoff]
                         [--format {json,csv}]
                         sample species

positional arguments:
  sample                sample id
  species               species

optional arguments:
  -h, --help            show this help message and exit
  -k kmer, --kmer kmer  kmer length (default:21)
  --tmp TMP             tmp directory (default: tmp/)
  --keep_tmp            Dont remove tmp files
  --skeleton_dir SKELETON_DIR
                        directory for skeleton binaries
  --mccortex31_path MCCORTEX31_PATH
                        Path to mccortex31. Default mccortex31
  -t THREADS, --threads THREADS
                        threads
  -m MEMORY, --memory MEMORY
                        memory for graph constuction
  --expected_depth EXPECTED_DEPTH
                        expected depth
  -1 seq [seq ...], --seq seq [seq ...]
                        sequence files (fasta,fastq,bam)
  -c ctx, --ctx ctx     cortex graph binary
  -f, --force           force
  --ont                 Set default for ONT data. Sets expected_error_rate to
                        0.15 and to haploid
  --guess_sequence_method
                        Guess if ONT or Illumia based on error rate. If error
                        rate is > 10%, ploidy is set to haploid and a
                        confidence threshold is used
  --ignore_minor_calls  Ignore minor calls when running resistance prediction
  --ignore_filtered IGNORE_FILTERED
                        don't include filtered genotypes
  --model model         Genotype model used, default kmer_count. Options
                        kmer_count, median_depth
  --ploidy ploidy       Use a diploid (includes 0/1 calls) or haploid
                        genotyping model
  --filters FILTERS [FILTERS ...]
                        don't include filtered genotypes
  --report_all_calls    report all calls
  --expected_error_rate EXPECTED_ERROR_RATE
                        Expected sequencing error rate. Set to 0.15 for ONT
                        genotyping.
  --min_variant_conf MIN_VARIANT_CONF
                        minimum genotype confidence for variant genotyping
  --min_gene_conf MIN_GENE_CONF
                        minimum genotype confidence for gene genotyping
  --min_proportion_expected_depth MIN_PROPORTION_EXPECTED_DEPTH
                        minimum depth required on the sum of both alleles.
                        Default 0.3 (30%)
  --min_gene_percent_covg_threshold MIN_GENE_PERCENT_COVG_THRESHOLD
                        all genes alleles found above this percent coverage
                        will be reported (default 100 (only best alleles
                        reported))
  --output OUTPUT       File path to save output json file as. Default is to
                        stdout.
  -q, --quiet           do not output warnings to stderr
  --panel panel         variant panel (default:201901). custom requires
                        custom_probe_set_path and
                        custom_variant_to_resistance_json to be set
  --custom_probe_set_path custom_probe_set_path
                        For use with `--panel custom`. File path to fasta file
                        from `mykrobe make-probes`.
  --custom_variant_to_resistance_json custom_variant_to_resistance_json
                        For use with `--panel custom`. File path to JSON with
                        key,value pairs of variant names and induced drug
                        resistance.
  --min_depth min_depth
                        min_depth
  --conf_percent_cutoff conf_percent_cutoff
                        Number between 0 and 100. Determines
                        --min_variant_conf, by simulating variants and
                        choosing the cutoff that would keep x% of the
                        variants. Default is 90 if --ont, otherwise
                        --min_variant_conf is used as the cutoff
  --format {json,csv}   Choose output format. Default: csv.
```

### Examples

```sh
mykrobe predict tb_sample_id tb -1 tb_sequence.bam/fq --format json --output results.json
# send output to stdout instead
mykrobe predict staph_sample_id staph -1 staph_sequence.bam/fq > result.csv
```

### Output

Output is in CSV by default. For a more detailed output use the JSON format with `--format json`.
```
{
    "sample_id": {
        "susceptibility": {
            "Rifampicin": {
                "predict": "S"
            },
            ...
            "Streptomycin": {
                "predict": "S"
            }
        "phylogenetics": {
            "lineage": {
                "Unknown": {
                    "percent_coverage": -1,
                    "median_depth": -1
                }
            },
            ...
            "species": {
                "Mycobacterium_tuberculosis": {
                    "percent_coverage": 98.0,
                    "median_depth": 53
                }
            }
        },
        "typed_variants": {
            "rpoB_N438S-AAC761118AGT": {
                "info": {
                    "contamination_depths": [],
                    "coverage": {
                        "alternate": {
                            "percent_coverage": 47.62,
                            "median_depth": 0.0,
                            "min_depth": 0
                        },
                        "reference": {
                            "percent_coverage": 100.0,
                            "median_depth": 49.0,
                            "min_depth": 44.0
                        }
                    },
                    "expected_depths": [
                        56.0
                    ]
                },
                "_cls": "Call.VariantCall",
                "genotype": [
                    0,
                    0
                ],
                "genotype_likelihoods": [
                    -4.25684443365591,
                    -99999999.0,
                    -99999999.0
                ]
            },   ...
        },
```
## Citing

If you use one of the following panels please cite the relevant publications:
```sh
mykrobe predict tb_sample_id  tb --panel walker-2015 -1 tb_sequence.bam
```

> Walker, Timothy M., et al. "Whole-genome sequencing for prediction of Mycobacterium tuberculosis drug susceptibility and resistance: a retrospective cohort study." The Lancet Infectious Diseases 15.10 (2015): 1193-1202.

```sh
mykrobe predict tb_sample_id  tb --panel bradley-2015 -1 tb_sequence.bam
```

> Bradley, Phelim, et al. "Rapid antibiotic-resistance predictions from genome sequence data for Staphylococcus aureus and Mycobacterium tuberculosis." Nature communications 6 (2015).

## Genotyping a pre-built probe set
```
mykrobe genotype --help
usage: mykrobe genotype [-h] [-k kmer] [--tmp TMP] [--keep_tmp]
                              [--skeleton_dir SKELETON_DIR]
                              [--mccortex31_path MCCORTEX31_PATH] [-t THREADS]
                              [-m MEMORY] [--expected_depth EXPECTED_DEPTH]
                              [-1 seq [seq ...]] [-c ctx] [-f] [--ont]
                              [--guess_sequence_method] [--ignore_minor_calls]
                              [--ignore_filtered IGNORE_FILTERED]
                              [--model model] [--ploidy ploidy]
                              [--filters FILTERS [FILTERS ...]]
                              [--report_all_calls]
                              [--expected_error_rate EXPECTED_ERROR_RATE]
                              [--min_variant_conf MIN_VARIANT_CONF]
                              [--min_gene_conf MIN_GENE_CONF]
                              [--min_proportion_expected_depth MIN_PROPORTION_EXPECTED_DEPTH]
                              [--min_gene_percent_covg_threshold MIN_GENE_PERCENT_COVG_THRESHOLD]
                              [--output OUTPUT] [-q]
                              sample probe_set

positional arguments:
  sample                sample id
  probe_set             probe_set

optional arguments:
  -h, --help            show this help message and exit
  -k kmer, --kmer kmer  kmer length (default:21)
  --tmp TMP             tmp directory (default: tmp/)
  --keep_tmp            Dont remove tmp files
  --skeleton_dir SKELETON_DIR
                        directory for skeleton binaries
  --mccortex31_path MCCORTEX31_PATH
                        Path to mccortex31. Default mccortex31
  -t THREADS, --threads THREADS
                        threads
  -m MEMORY, --memory MEMORY
                        memory for graph constuction
  --expected_depth EXPECTED_DEPTH
                        expected depth
  -1 seq [seq ...], --seq seq [seq ...]
                        sequence files (fasta,fastq,bam)
  -c ctx, --ctx ctx     cortex graph binary
  -f, --force           force
  --ont                 Set demykrobe genotype --help
usage: mykrobe genotype [-h] [-k kmer] [--tmp TMP] [--keep_tmp]
                              [--skeleton_dir SKELETON_DIR]
                              [--mccortex31_path MCCORTEX31_PATH] [-t THREADS]
                              [-m MEMORY] [--expected_depth EXPECTED_DEPTH]
                              [-1 seq [seq ...]] [-c ctx] [-f] [--ont]
                              [--guess_sequence_method] [--ignore_minor_calls]
                              [--ignore_filtered IGNORE_FILTERED]
                              [--model model] [--ploidy ploidy]
                              [--filters FILTERS [FILTERS ...]]
                              [--report_all_calls]
                              [--expected_error_rate EXPECTED_ERROR_RATE]
                              [--min_variant_conf MIN_VARIANT_CONF]
                              [--min_gene_conf MIN_GENE_CONF]
                              [--min_proportion_expected_depth MIN_PROPORTION_EXPECTED_DEPTH]
                              [--min_gene_percent_covg_threshold MIN_GENE_PERCENT_COVG_THRESHOLD]
                              [--output OUTPUT] [-q]
                              sample probe_set

positional arguments:
  sample                sample id
  probe_set             probe_set

optional arguments:
  -h, --help            show this help message and exit
  -k kmer, --kmer kmer  kmer length (default:21)
  --tmp TMP             tmp directory (default: tmp/)
  --keep_tmp            Dont remove tmp files
  --skeleton_dir SKELETON_DIR
                        directory for skeleton binaries
  --mccortex31_path MCCORTEX31_PATH
                        Path to mccortex31. Default mccortex31
  -t THREADS, --threads THREADS
                        threads
  -m MEMORY, --memory MEMORY
                        memory for graph constuction
  --expected_depth EXPECTED_DEPTH
                        expected depth
  -1 seq [seq ...], --seq seq [seq ...]
                        sequence files (fasta,fastq,bam)
  -c ctx, --ctx ctx     cortex graph binary
  -f, --force           force
  --ont                 Set default for ONT data. Sets expected_error_rate to
                        0.15 and to haploid
  --guess_sequence_method
                        Guess if ONT or Illumia based on error rate. If error
                        rate is > 10%, ploidy is set to haploid and a
                        confidence threshold is used
  --ignore_minor_calls  Ignore minor calls when running resistance prediction
  --ignore_filtered IGNORE_FILTERED
                        don't include filtered genotypes
  --model model         Genotype model used, default kmer_count. Options
                        kmer_count, median_depth
  --ploidy ploidy       Use a diploid (includes 0/1 calls) or haploid
                        genotyping model
  --filters FILTERS [FILTERS ...]
                        don't include filtered genotypes
  --report_all_calls    report all calls
  --expected_error_rate EXPECTED_ERROR_RATE
                        Expected sequencing error rate. Set to 0.15 for ONT
                        genotyping.
  --min_variant_conf MIN_VARIANT_CONF
                        minimum genotype confidence for variant genotyping
  --min_gene_conf MIN_GENE_CONF
                        minimum genotype confidence for gene genotyping
  --min_proportion_expected_depth MIN_PROPORTION_EXPECTED_DEPTH
                        minimum depth required on the sum of both alleles.
                        Default 0.3 (30%)
  --min_gene_percent_covg_threshold MIN_GENE_PERCENT_COVG_THRESHOLD
                        all genes alleles found above this percent coverage
                        will be reported (default 100 (only best alleles
                        reported))
  --output OUTPUT       File path to save output json file as. Default is to
                        stdout.
  -q, --quiet           do not output warnings to stderrfault for ONT data. Sets expected_error_rate to
                        0.15 and to haploid
  --guess_sequence_method
                        Guess if ONT or Illumia based on error rate. If error
                        rate is > 10%, ploidy is set to haploid and a
                        confidence threshold is used
  --ignore_minor_calls  Ignore minor calls when running resistance prediction
  --ignore_filtered IGNORE_FILTERED
                        don't include filtered genotypes
  --model model         Genotype model used, default kmer_count. Options
                        kmer_count, median_depth
  --ploidy ploidy       Use a diploid (includes 0/1 calls) or haploid
                        genotyping model
  --filters FILTERS [FILTERS ...]
                        don't include filtered genotypes
  --report_all_calls    report all calls
  --expected_error_rate EXPECTED_ERROR_RATE
                        Expected sequencing error rate. Set to 0.15 for ONT
                        genotyping.
  --min_variant_conf MIN_VARIANT_CONF
                        minimum genotype confidence for variant genotyping
  --min_gene_conf MIN_GENE_CONF
                        minimum genotype confidence for gene genotyping
  --min_proportion_expected_depth MIN_PROPORTION_EXPECTED_DEPTH
                        minimum depth required on the sum of both alleles.
                        Default 0.3 (30%)
  --min_gene_percent_covg_threshold MIN_GENE_PERCENT_COVG_THRESHOLD
                        all genes alleles found above this percent coverage
                        will be reported (default 100 (only best alleles
                        reported))
  --output OUTPUT       File path to save output json file as. Default is to
                        stdout.
  -q, --quiet           do not output warnings to stderr
```

### Examples
```
mykrobe genotype sample_id example-data/staph-amr-bradley_2015.fasta -1 seq.fq
{
    "sample_id": {
    "files": [
        "seq.fq "
    ],
    "kmer": 21,
    "sequence_calls": {
        "mecA": {
            "info": {
                "copy_number": 0.0,
                "contamination_depths": [],
                "coverage": {
                    "percent_coverage": 0.0,
                    "median_depth": 0.0,
                    "min_non_zero_depth": 0.0
                },
                "expected_depths": [
                    1
                ]
            },
            "_cls": "Call.SequenceCall",
            "genotype": [
                0,
                0
            ],
            "genotype_likelihoods": [
                -0.001,
                -99999999.0,
                -99999999.0
            ]
        },
        "fusA": {
            "info": {
                "copy_number": 1.0276923076923077,
                "contamination_depths": [],
                "version": "10",
                "coverage": {
                    "percent_coverage": 100.0,
                    "median_depth": 167.0,
                    "min_non_zero_depth": 116.0
                },
                "expected_depths": [
                    162.5
                ]
            },
            "_cls": "Call.SequenceCall",
            "genotype": [
                1,
                1
            ],
            "genotype_likelihoods": [
                -994.7978064088725,
                -349.45246450237215,
                -10.95808091830304
            ]
        },
    ....
}
}
```
### Make a custom probe set (for use with `mykrobe genotype`)

### Add variants to the database (for background/context)

This is optional but will make any probe sets built more robust to variation in within k-1 bases of the key variants. This will require [mongoDB](https://www.mongodb.com/) > 3.0 running in the background.
```
usage: mykrobe variants add [-h] [--db_name db_name] [-f] [-q]
                                  [-m METHOD]
                                  vcf reference_set

positional arguments:
  vcf                   a vcf file
  reference_set         reference set

optional arguments:
  -h, --help            show this help message and exit
  --db_name db_name     db_name
  -f, --force           force
  -q, --quiet           do not output warnings to stderr
  -m METHOD, --method METHOD
                        variant caller method (e.g. CORTEX)
```
To add a VCF to the database db_name run
```sh
mykrobe variants add --db_name :db_name sample.vcf :reference
```

Use the `--method` argument to specify the variant caller or pipeline used (if you'll have multiple Call Sets per sample)

```sh
mykrobe variants add --db_name :db_name --method CORTEX sample_cortex.vcf :reference
```

### Make probes and dump-probes
```
    mykrobe variants make-probes --help
    usage: mykrobe variants make-probes [-h] [--db_name db_name] [-q]
                                              [-f VCF] [-v VARIANT] [-t TEXT_FILE]
                                              [-g GENBANK] [-k KMER]
                                              [--no-backgrounds]
                                              reference_filepath

    positional arguments:
      reference_filepath    reference_filepath

    optional arguments:
      -h, --help            show this help message and exit
      --db_name db_name     db_name
      -q, --quiet           do not output warnings to stderr
      -f VCF, --vcf VCF     Use variants defined in a VCF file
      -v VARIANT, --variant VARIANT
                            Variant in DNA positions e.g. A1234T
      -t TEXT_FILE, --text_file TEXT_FILE
                            Text file containing variants as rows A1234T
      -g GENBANK, --genbank GENBANK
                            Genbank file containing genes as features
      -k KMER, --kmer KMER  kmer length
      --no-backgrounds      Build probe set against reference only ignoring nearby
                            variants
```
### Examples

#### 1 Simple case - building a probe without using backgrounds

    mykrobe variants make-probes -v A1234T example-data/NC_000962.3.fasta

#### 2. 'Dumping' the Variant database

To build a ProbeSet of all non-singleton variants in the database run:
```
`mykrobe variants dump-probes`

usage: mykrobe dump-probes [-h] [--db_name db_name] [-q] [--kmer kmer] [--force]
                     [-v]
                     reference_filepath

positional arguments:
  reference_filepath  reference_filepath

optional arguments:
  -h, --help          show this help message and exit
  --db_name db_name   db_name
  -q, --quiet         do not output warnings to stderr
  --kmer kmer         kmer length
  --force
  -v, --verbose
```

```sh
mykrobe variants dump-probes reference_set.fasta > variant_probe_set.fasta
```

This will generate a probe set for each variant in the database. The resulting fasta file will look like the following:
```
     >ref-37d2eea6a23d526cbee4e00b901dc97885a88e7aa8721432b080dcc342b459ce?num_alts=10&ref=56cf2e4ca9fefcd2b15de4d6
    TCGCCGCAGCGGTTGGCAACGATGTGGTGCGATCGCTAAAGATCACCGGGCCGGCGGCACCAT
    ...
    TCGCCGCAGCGGTTGGCAACGATGTGGTGCAATCGCTAAAGATCACCGGGCCGGCGGCATCAT
    >alt-37d2eea6a23d526cbee4e00b901dc97885a88e7aa8721432b080dcc342b459ce
    TCGCCGCAGCGGTTGGCAACGATGTGGTGCAATCGCTAAAGATCACCGGGCCGGCGGCACGAT
    >ref-2dab6387a677ac17f6bc181f47235a4196885723b34ceff3a05ffcbfd6834347?num_alts=10&ref=56cf2e4ca9fefcd2b15de4d6
    CTGTCGCTGGGAAGAGCGAATACGTCTGGACCAGGACGGGCTACCCGAACACGATATCTTTCG
    >alt-2dab6387a677ac17f6bc181f47235a4196885723b34ceff3a05ffcbfd6834347
    ...
```
Where you have a series of variants represented as a set of alleles. The reference allele followed by multiple alternate alleles. You will end up with multiple alternate alleles if there are other variants that fall within k of the target variant.

Each variant is referenced by a `var_hash` with is the hash of ":ref:pos:alt" which is indexed in the database and can be used to query for Variant object.

See `mykrobe genotype` to use these probes to genotype a new sample.

#### 3. Building a custom probe set

`mykrobe variants  make-probes` allows you to build a probe set using Variants that are not already in the database but using the population variation to produce multiple alleles per variant.
```
     usage: mykrobe variants  make-probes [-h] [--db_name db_name] [-q] [-v VARIANT] [-f FILE]
                             [-g GENBANK] [-k KMER] [--no-backgrounds]
                             reference_filepath

    positional arguments:
      reference_filepath    reference_filepath

    optional arguments:
      -h, --help            show this help message and exit
      --db_name db_name     db_name
      -q, --quiet           do not output warnings to stderr
      -v VARIANT, --variant VARIANT
                            Variant in DNA positions e.g. A1234T
      -f FILE, --file FILE  File containing variants as rows A1234T
      -g GENBANK, --genbank GENBANK
                            Genbank file containing genes as features
      -k KMER, --kmer KMER  kmer length
      --no-backgrounds      Build probe set against reference only ignoring nearby
                            variants
```
##### Build a variant probe set defined based on reference co-ordinates (1-based)

First, define your variants for which you want to build probes. Columns are
```
ref/gene pos ref alt alphabet

    ref     2522798 G       T       DNA
    ref     3785555 A       G       DNA
    ref     839793  C       A       DNA
    ref     2734398 C       G       DNA
    ref     3230861 T       A       DNA
    ref     1018694 A       T       DNA
```

```
mykrobe variants make-probes --db_name :db_name -f variants.txt ref.fa > variant_probe_set.fa
```
##### Build a variant probe set defined based on gene co-ordinates (1-based)

You can also define your variants in terms of gene coordinates in amino acid or DNA space.
```
    rpoB    S431X   PROT
    rpoB    F425X   PROT
    embB    M306X   PROT
    rrs     C513X   DNA
    gyrA    D94X    PROT
    gid     P75L    PROT
    gid     V88A    PROT
    katG    S315X   PROT
```
To do this you must provide a genbank file defining the position of the variants in the reference (-g (GENBANK) )
```
mykrobe variants  make-probes --db_name :db_name -f aa_variants.txt -g ref.gb  ref.fa> gene_variant_probe_set.fa
```
# Citation

> [Bradley, Phelim, et al. "Rapid antibiotic-resistance predictions from genome sequence data for Staphylococcus aureus and Mycobacterium tuberculosis."Nature communications 6 (2015).](http://www.nature.com/ncomms/2015/151221/ncomms10063/full/ncomms10063.html)

Please cite us if you use Mykrobe predictor in a publication

# Tests

To run tests:

Requires `mongod`. [Download](https://www.mongodb.com/download-center/community).

    mongod

and in another window.

      pip install tox
      tox

To run tests for a particular python version run, e.g. python 3.6:

    tox -e py36

# Release

See dist/README.md
