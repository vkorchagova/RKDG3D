#!/bin/bash

cd ${MFEM_ROOT}
make -j 4 install 
cd ${MFEM_APP_ROOT}/build
make -j 4