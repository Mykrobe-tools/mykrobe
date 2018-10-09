# Procedure for generating binaries for Windows 64 bits

Working directory: `mykrobe-atlas-cli`

1) modify `./mccortex/libs/xxHash/xxhsum.c` according to https://github.com/Cyan4973/xxHash/issues/100

2) with cygwin64: `(cd ./mccortex && make)`

3) with cygwin64: `pip install pyinstaller`

4) with cygwin64: `(cd dist && pyinstaller --workpath='./pyinstaller_build/binary_cache' --distpath='./pyinstaller_build' mykrobe_atlas_pyinstaller.spec)`

5) binaries and dependencies: `./dist/pyinstaller_build/mykrobe_atlas`

