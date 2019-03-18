# Procedure for generating binaries for Windows 64 bits

Working directory: `mykrobe-atlas-cli`

1) modify `./mccortex/libs/xxHash/xxhsum.c` according to https://github.com/Cyan4973/xxHash/issues/100

2) with cygwin64: `(cd ./mccortex && make)`

3) with cygwin64: `pip install pyinstaller`

4) with cygwin64: `(cd dist && pyinstaller --workpath='./pyinstaller_build/binary_cache' --distpath='./pyinstaller_build' mykrobe_atlas_pyinstaller.spec)`

5) binaries and dependencies: `./dist/pyinstaller_build/mykrobe_atlas`

# Test dist

for f in /Users/Phelimb/Dropbox/Atlas/test_data/exemplar_seqeuence_data/*.fastq*
do
ff=`basename $f`
./mykrobe_atlas predict --mccortex31_path ./mccortex31 --force non-mycobacterium_saur.fastq tb -1 $f  --format json --guess_sequence_method --quiet > ../$ff.json 
done

for f in /Users/Phelimb/Dropbox/Atlas/test_data/exemplar_seqeuence_data/*.bam
do
ff=`basename $f`
./mykrobe_atlas predict --mccortex31_path ./mccortex31 --force non-mycobacterium_saur.fastq tb -1 $f  --format json --guess_sequence_method --quiet > ../$ff.json
done

