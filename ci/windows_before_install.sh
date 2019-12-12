export PATH=/c/tools/cygwin/bin:$PATH
original_dir=$PWD
cd mccortex/libs/htslib/
git checkout 1832d3a1b75133e55fb6abffc3f50f8a6ed5ceae
make
cd $original_dir/mccortex/
make

