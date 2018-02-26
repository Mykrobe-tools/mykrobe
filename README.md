Master: [![Build Status](https://travis-ci.org/Phelimb/mykrobe.svg?branch=master)](https://travis-ci.com/Phelimb/mykrobe)

Dev: [![Build Status](https://travis-ci.org/Phelimb/mykrobe.svg?branch=dev)](https://travis-ci.com/Phelimb/mykrobe)

## Installation

	git clone --recursive git@github.com:Phelimb/mykrobe.git    
	python setup.py install

OR

	(sudo) pip install git+https://github.com/Phelimb/mykrobe

We recommend that you use virtualenv to install mykrobe - see instructions below:

### Install mykrobe with virtualenv (recommended but optional)

#### Install virtualenv

Follow instruction at https://virtualenv.readthedocs.org/en/latest/installation.html

#### Create virtualenv 

	virtualenv venv

#### Activate the virtualenv

	source venv/bin/activate

You can deactivate at anytime by typing `deactivate` in your terminal. 

#### Install mykrobe

	pip install git+https://github.com/Phelimb/mykrobe


## Usage

	usage: mykrobe [-h] [--version] {add,dump-probes,make-probes,genotype} ...

	optional arguments:
	  -h, --help            show this help message and exit
	  --version             mykrobe version

	[sub-commands]:
	  {add,dump-probes,make-probes,genotype}
	    add                 adds a set of variants to the mykrobe
	    dump-probes         dump a panel of variant alleles
	    make-probes         make probes from a list of variants
	    genotype            genotype a sample using a probe set

### Make probes

```
(venv)-bash-4.1$ ./mykrobe/mykrobe_main.py make-probes --help
usage: mykrobe make-probes [-h] [--db_name db_name] [-q] [-v VARIANT] [-f FILE]
                         [-g GENBANK] [-k KMER]
                         reference_filepath

positional arguments:
  reference_filepath    file path to reference

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
```

### Examples


#### 1 Simple case - building a probe without using backgrounds

	mykrobe make-probes -v A1234T example-data/NC_000962.3.fasta

#### 2. 'Dumping' the Variant database

`mykrobe dump-probes`

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


	 mykrobe dump-probes reference_set.fasta > variant_probe_set.fasta 
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

`mykrobe make-probes` allows you to build a probe set using Variants that are not already in the database. 

	 usage: mykrobe make-probes [-h] [--db_name db_name] [-q] [-v VARIANT] [-f FILE]
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


	 mykrobe make-probes --db_name :db_name -f variants.txt ref.fa > variant_probe_set.fa 
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

	 mykrobe make-probes --db_name :db_name -f aa_variants.txt -g ref.gb  ref.fa> gene_variant_probe_set.fa 


 


### Tests
To run tests:
    
    pip install tox
    tox

To run tests for a particular python version run, e.g. python 3.6:

    tox -e py36

### Citation

Please cite us if you use mykrobe in a publication:

<U>Bradley P</U>, ... , Iqbal Z.
Rapid antibiotic-resistance predictions from genome sequence data for _Staphylococcus aureus_ and _Mycobacterium tuberculosis_.
**Nature Communications** 2015 Dec 21;6:10063<BR>
PMID: [26686880](https://www.ncbi.nlm.nih.gov/pubmed/26686880)<BR>
DOI: [10.1038/ncomms10063](http://www.nature.com/ncomms/2015/151221/ncomms10063/full/ncomms10063.html)

All analysis in this paper was done with release [v0.1.3-beta](https://github.com/iqbal-lab/Mykrobe-predictor/releases/tag/v0.1.3-beta).


