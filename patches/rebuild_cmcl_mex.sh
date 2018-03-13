#!/bin/bash
# rebuild the mex from the mwrap-auto-gen .c files. Barnett 3/13/18
# ensure to use gcc 6.4.0 for consistency.

# note lib needs the direct ones too
gfortran -fPIC -Ofast -funroll-loops -march=native -c nufft*df90.f dirft*d.f dfftpack.f next235.f

mex -largeArrayDims -lrt -lm -lgfortran nufft3d.c nufft*df90.o dirft*d.o dfftpack.o next235.o
mex -largeArrayDims -lrt -lm -lgfortran nufft2d.c nufft*df90.o dirft*d.o dfftpack.o next235.o
mex -largeArrayDims -lrt -lm -lgfortran nufft1d.c nufft*df90.o dirft*d.o dfftpack.o next235.o
