#!/bin/sh

export OMP_NUM_THREADS=8

../../bin/ccx_2.20_MT_SLS_umat pad_nofillet_real_impactor_ccx >> stdout.log
