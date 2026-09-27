#!/bin/bash
set -e

ENV_NAME=molektra_env
PYTHON_VERSION=3.9

eval "$(conda shell.bash hook)"

echo "--- Erasing former environment, if any ---"
conda remove -n $ENV_NAME --all -y

echo "--- Creating environment with Python, PyQt5 and gfortran via conda ---"
conda create -n $ENV_NAME python=$PYTHON_VERSION pyqt gfortran meson ninja pkg-config -y

conda activate $ENV_NAME

echo "--- Fixing setuptools and numpy ---"
pip install "setuptools<60" wheel
pip install "numpy<2"

echo "--- Installing base dependencies ---"
pip install pyqtgraph PyOpenGL matplotlib scipy PyYAML posym

echo "--- Installing symgroupy from GitHub (PyPI tarball missing Fortran sources) ---"
TMPDIR=$(mktemp -d)
curl -L https://github.com/abelcarreras/symgroup/archive/refs/heads/master.zip \
     -o "$TMPDIR/symgroupy.zip"
unzip -q "$TMPDIR/symgroupy.zip" -d "$TMPDIR"
touch "$TMPDIR/symgroup-master/python/readme.md"
pip install "$TMPDIR/symgroup-master/python"

echo "--- Installing wfnsympy from GitHub (PyPI tarball missing include dir) ---"

TMPDIR=$(mktemp -d)            
curl -L https://github.com/abelcarreras/WFNSYM/archive/refs/heads/master.zip \
     -o "$TMPDIR/wfnsympy.zip"
unzip -q "$TMPDIR/wfnsympy.zip" -d "$TMPDIR"
pip install "$TMPDIR/WFNSYM-master/python"
rm -rf "$TMPDIR"


echo "--- Installing remaining cosymlib dependencies ---"
pip install pointgroup huckelpy

echo "--- Installing cosymlib (precompiled wheel for arm64, skip broken deps) ---"
pip install cosymlib --no-deps

echo ""
echo " Environment '$ENV_NAME' is ready."
echo " To verify:"
echo "   conda activate $ENV_NAME"
echo "   python -c 'import cosymlib; print(\"Cosymlib loaded correctly\")'"
