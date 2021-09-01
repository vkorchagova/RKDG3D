#!/bin/bash

cd ${MFEM_ROOT}
make -j 4 DESTDIR=${MFEM_ROOT} install 
cd ${MFEM_APP_ROOT}/build
make -j 4