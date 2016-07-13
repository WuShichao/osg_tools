# -*- mode: python -*-

block_cipher = None


a = Analysis(['megaplot.py'],
             pathex=['/home/jclark308/src/lscsoft/bayeswave/trunk/postprocess'],
             binaries=None,
             datas=[('./navigate.js','.'),
             ('./BWBweb.css','.'),
             ('./secure_ajax.js','.'),
             ('./svn_info.txt','.')],
             hiddenimports=['scipy.linalg'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='megaplot',
          debug=False,
          strip=True,
          upx=True,
          console=True )
