#!/bin/bash

# $PYTHON -m pip install -r requirements.txt
$PYTHON -m pip install .
# $PYTHON setup.py install

CMAT="${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}"
mkdir -p ${PREFIX}/bin ${CMAT}

chmod 775 bin/cmat/*
cp bin/cmat/* ${PREFIX}/bin

mv bin/ mappings/ pipelines/ ${CMAT}