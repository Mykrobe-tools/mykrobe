ls /
export PATH=/cygdrive/c/tools/cygwin/bin:$PATH
pwd
original_dir=$PWD
git clone --recursive -b geno_kmer_count https://github.com/phelimb/mccortex
cd mccortex/libs/htslib/
git checkout 1832d3a1b75133e55fb6abffc3f50f8a6ed5ceae
make
cd $original_dir/mccortex/
make

