# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['aqmeasy.py'],
    pathex=[],
    binaries=[],
    datas=[('ui/resources/*', 'ui/resources'),
           ('ui/*.py', 'ui'),
           ('ui/resources/*.py', 'ui/resources'),
           ('models/*.py', 'models'),
           ('controllers/*.py', 'controllers'),
           ('utils.py', '.')
           ],
    hiddenimports=['rdkit', 'rdkit.Chem.Descriptrs'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='aqmeasy',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['ui/resources/icons.iconset/aqme_icon.icns'],
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='aqmeasy',
)
app = BUNDLE(
    coll,
    name='aqmeasy.app',
    icon='ui/resources/icons.iconset/aqme_icon.icns',
    bundle_identifier=None,
)
