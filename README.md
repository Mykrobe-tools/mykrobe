# Mykrobe
Antibiotic resistance prediction in minutes. 

Currently supports _Mycobacterium tuberculosis_, _Staphylococcus aureus_, _Shigella sonnei_, _Salmonella typhi_.

<http://www.mykrobe.com>


## Documentation

Please see the [mykrobe wiki](https://github.com/Mykrobe-tools/mykrobe/wiki) for documentation.
For _S. typhi_ also see https://github.com/katholt/genotyphi 


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
