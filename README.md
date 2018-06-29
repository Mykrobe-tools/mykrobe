# Procedure for generating binaries for Windows 64 bits

Working directory: `Mykrobe-predictor`

1) modify `./mccortex/libs/xxHash/xxhsum.c` according to https://github.com/Cyan4973/xxHash/issues/100

2) with cygwin64: `(cd ./mccortex && make)`

3) with cygwin64: `pip install git+https://github.com/Phelimb/atlas`

4) with cygwin64: `pip install pyinstaller`

5) with cygwin64: `(cd dist && pyinstaller --workpath='./pyinstaller_build/binary_cache' --distpath='./pyinstaller_build' mykrobe_predictor_pyinstaller.spec)`

6) binaries and dependencies: `./dist/pyinstaller_build/mykrobe_predictor`