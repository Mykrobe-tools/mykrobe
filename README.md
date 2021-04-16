[![Travis build Status](https://travis-ci.com/Mykrobe-tools/mykrobe.svg?branch=master)](https://travis-ci.org/Mykrobe-tools/mykrobe)

# Mykrobe

<http://www.mykrobe.com>


## Documentation

Please see the [mykrobe wiki](https://github.com/Mykrobe-tools/mykrobe/wiki) for documentation.


## Quick start

**Install**:

* bioconda - `conda install -c bioconda mykrobe`
* from source - `pip3 install . && mykrobe panels update_metadata && mykrobe panels update_species all`
* or using singularity or docker (see wiki for details)

**Run** on Mtb, making a JSON file of results:

```
mykrobe predict --sample my_sample_name \
  --species tb \
  --output out.json \
  --format json \
  --seq reads.fq.gz
```


Test reads can be obtained by running:

```
wget -O reads.fq.gz https://ndownloader.figshare.com/files/21059229
```
