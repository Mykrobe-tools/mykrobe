# -*- mode: python -*-

block_cipher = None


a = Analysis(['../src/mykrobe/cli.py'],
             pathex=['../../mykrobe-atlas-cli'],
             binaries=[],
             datas=[('C:\cygwin64\bin\*.dll', ".")],
             hiddenimports=['mykrobe'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

a.binaries += Tree('../mccortex/bin/', prefix='.')
a.datas += Tree('../src/mykrobe/data', prefix='mykrobe/data')

pyz = PYZ(a.pure, a.zipped_data,
          cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='mykrobe_atlas',
          debug=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='mykrobe_atlas')
