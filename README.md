Master: [![Build Status](https://travis-ci.org/Phelimb/mykrobe-atlas-cli.svg?branch=master)](https://travis-ci.org/Phelimb/mykrobe-atlas-cli)

Dev: [![Build Status](https://travis-ci.org/Phelimb/mykrobe-atlas-cli.svg?branch=dev)](https://travis-ci.org/Phelimb/mykrobe-atlas-cli)

Tested on python 2.7, 3.4, 3.5, and 3.6.

## Requirements

* python 2.7, python > 3.4
* [mongodb](https://www.mongodb.com/) > 3.0 (optional)


## Installation

	git clone https://github.com/Mykrobe-tools/mykrobe-atlas-cli.git mykrobe
	cd mykrobe
	
	## Download pre-built probesets
	wget -O mykrobe-data.tar.gz https://goo.gl/DXb9hN && tar -zxvf mykrobe-data.tar.gz && rm -fr src/mykrobe/data && mv mykrobe-data src/mykrobe/data
	
	pip install .
	

This will install two executables: mykrobe and mccortex31 (a fork of [mccortex](https://github.com/mcveanlab/mccortex)).

## Usage

	mykrobe --help
	usage: mykrobe [-h] [--version] {predict,variants,vars,genotype} ...

	optional arguments:
	  -h, --help            show this help message and exit
	  --version             mykrobe-atlas version

	[sub-commands]:
	  {predict,variants,vars,genotype}
	    predict             predict the sample's drug susceptibility
	    variants (vars)     build variant probes
	    genotype            genotype a sample using a probe set

### AMR prediction (Mykrobe predictor)

	mykrobe predict --help
	usage: mykrobe-atlas predict [-h] [-k kmer] [--tmp TMP] [--keep_tmp]
                             [--skeleton_dir SKELETON_DIR]
                             [--mccortex31_path MCCORTEX31_PATH] [-t THREADS]
                             [-m MEMORY] [--expected_depth EXPECTED_DEPTH]
                             [-1 seq [seq ...]] [-c ctx] [-f] [--ont]
                             [--ignore_filtered IGNORE_FILTERED]
                             [--model model] [--filters FILTERS [FILTERS ...]]
                             [--report_all_calls]
                             [--expected_error_rate EXPECTED_ERROR_RATE]
                             [--min_variant_conf MIN_VARIANT_CONF]
                             [--min_gene_conf MIN_GENE_CONF]
                             [--min_gene_percent_covg_threshold MIN_GENE_PERCENT_COVG_THRESHOLD]
                             [-q] [--panel panel] [--min_depth min_depth]
                             [--output OUTPUT]
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
      --ont                 Set default for ONT data
      --ignore_filtered IGNORE_FILTERED
                            don't include filtered genotypes
      --model model         Genotype model used, default median_depth. Options
                            kmer_count, median_depth
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
      --min_gene_percent_covg_threshold MIN_GENE_PERCENT_COVG_THRESHOLD
                            all genes alleles found above this percent coverage
                            will be reported (default 100 (only best alleles
                            reported))
      -q, --quiet           do not output warnings to stderr
      --panel panel         variant panel (default:walker-2015)
      --min_depth min_depth
                            min_depth
      --output OUTPUT       File path to save output json file as. Default is to
                            stdout.

#### Examples

```bash
mykrobe predict tb_sample_id tb -1 tb_sequence.bam/fq --output results.json
# send output to stdout instead
mykrobe predict staph_sample_id staph -1 staph_sequence.bam/fq
```
	


e.g.

	mykrobe predict ERR117639 /download/ena/ERR117639*.gz tb
	
### Output

Output is in JSON format. To convert to a less verbose tabular format use [json_to_tsv](scripts/json_to_tsv.py).

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
	                            "min_depth": 47.0
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

#### Citations

If you use one of the following panels please cite the relevant publications:

	mykrobe predict tb_sample_id  tb --panel walker-2015 -1 tb_sequence.bam

> Walker, Timothy M., et al. "Whole-genome sequencing for prediction of Mycobacterium tuberculosis drug susceptibility and resistance: a retrospective cohort study." The Lancet Infectious Diseases 15.10 (2015): 1193-1202.

	mykrobe predict tb_sample_id  tb --panel bradley-2015 -1 tb_sequence.bam

> Bradley, Phelim, et al. "Rapid antibiotic-resistance predictions from genome sequence data for Staphylococcus aureus and Mycobacterium tuberculosis." Nature communications 6 (2015).	


### Genotyping a pre-built probe set

	mykrobe genotype --help
	usage: mykrobe genotype [-h] [-k kmer] [--tmp TMP] [--keep_tmp]
	                              [--skeleton_dir SKELETON_DIR]
	                              [--mccortex31_path MCCORTEX31_PATH] [-t THREADS]
	                              [-m MEMORY] [--expected_depth EXPECTED_DEPTH]
	                              [-1 seq [seq ...]] [-c ctx] [-f] [--ont]
	                              [--ignore_filtered IGNORE_FILTERED]
	                              [--model model]
	                              [--filters FILTERS [FILTERS ...]]
	                              [--report_all_calls]
	                              [--expected_error_rate EXPECTED_ERROR_RATE]
	                              [--min_variant_conf MIN_VARIANT_CONF]
	                              [--min_gene_conf MIN_GENE_CONF]
	                              [--min_gene_percent_covg_threshold MIN_GENE_PERCENT_COVG_THRESHOLD]
	                              [-q]
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
	  --ont                 Set default for ONT data
	  --ignore_filtered IGNORE_FILTERED
	                        don't include filtered genotypes
	  --model model         Genotype model used, default median_depth. Options
	                        kmer_count, median_depth
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
	  --min_gene_percent_covg_threshold MIN_GENE_PERCENT_COVG_THRESHOLD
	                        all genes alleles found above this percent coverage
	                        will be reported (default 100 (only best alleles
	                        reported))
	  -q, --quiet           do not output warnings to stderr
	  
#### Examples

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

	
## Make a custom probe set (for use with `mykrobe genotype`)

### Add variants to the database (for background/context)

This is optional but will make any probe sets built more robust to variation in within k-1 bases of the key variants. This will require [mongoDB](https://www.mongodb.com/) > 3.0 running in the background.

	usage: mykrobe-atlas variants add [-h] [--db_name db_name] [-f] [-q]
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

To add a VCF to the database db_name run

	mykrobe variants add --db_name :db_name sample.vcf :reference
	
Use the --method argument to specify the variant caller or pipeline used (if you'll have multiple Call Sets per sample)

	mykrobe variants add --db_name :db_name --method CORTEX sample_cortex.vcf :reference 	

### Make probes and dump-probes

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

### Examples


#### 1 Simple case - building a probe without using backgrounds

	mykrobe variants make-probes -v A1234T example-data/NC_000962.3.fasta

#### 2. 'Dumping' the Variant database

To build a ProbeSet of all non-singleton variants in the database run:

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


	 mykrobe variants dump-probes reference_set.fasta > variant_probe_set.fasta 
This will generate a probe set for each variant in the database. The resulting fasta file will look like the following:

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
Where you have a series of variants represented as a set of alleles. The reference allele followed by multiple alternate alleles. You will end up with multiple alternate alleles if there are other variants that fall within k of the target variant. 



Each variant is referenced by a `var_hash` with is the hash of ":ref:pos:alt" which is indexed in the database and can be used to query for Variant object.

See `mykrobe genotype` to use these probes to genotype a new sample.

#### 3. Building a custom probe set

`mykrobe variants  make-probes` allows you to build a probe set using Variants that are not already in the database but using the population variation to produce multiple alleles per variant. 

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

##### Build a variant probe set defined based on reference co-ordinates (1-based)



First, define your variants for which you want to build probes. Columns are 



ref/gene pos ref alt alphabet

	ref     2522798 G       T       DNA
	ref     3785555 A       G       DNA
	ref     839793  C       A       DNA
	ref     2734398 C       G       DNA
	ref     3230861 T       A       DNA
	ref     1018694 A       T       DNA 


	 mykrobe variants make-probes --db_name :db_name -f variants.txt ref.fa > variant_probe_set.fa 
##### Build a variant probe set defined based on gene co-ordinates (1-based)



You can also define your variants in terms of gene coordinates in amino acid or DNA space.

	rpoB    S431X   PROT
	rpoB    F425X   PROT
	embB    M306X   PROT
	rrs     C513X   DNA
	gyrA    D94X    PROT
	gid     P75L    PROT
	gid     V88A    PROT
	katG    S315X   PROT 
To do this you must provide a genbank file defining the position of the variants in the reference (-g (GENBANK) )

	 mykrobe variants  make-probes --db_name :db_name -f aa_variants.txt -g ref.gb  ref.fa> gene_variant_probe_set.fa 


 
# Citation 

> [Bradley, Phelim, et al. "Rapid antibiotic-resistance predictions from genome sequence data for Staphylococcus aureus and Mycobacterium tuberculosis."Nature communications 6 (2015).](http://www.nature.com/ncomms/2015/151221/ncomms10063/full/ncomms10063.html)

Please cite us if you use Mykrobe predictor in a publication



# Tests
To run tests:
    
    pip install tox
    tox

To run tests for a particular python version run, e.g. python 3.6:

    tox -e py36

