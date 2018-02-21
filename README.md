# nufft-bench

Codes to generate figures and tables
comparing performance of various NUFFT codes.
Written in MATLAB/octave.

Alex Barnett and Jeremy Magland; contributions by Joakim Anden.

January 2018.

## Dependencies

Tested on: MATLAB R2016b

Some codes require: https://github.com/ahbarnett/memorygraph



## Installing various NUFFT packages

Here's how to set up the NUFFT codes so that `nufft-bench` works.
Each package directory should be placed within the parent directory
of `nufft-bench`, ie alongside `nufft-bench`.
If you want to change this location, edit
`setuppaths.m`.

### Installing NFFT

Download current source from here:

https://www-user.tu-chemnitz.de/~potts/nfft/download.php

Extract to directory `nfft-3.3.2`, which should be placed alongside `nufft-bench`
as explained above.

`cd` to `nfft-3.3.2` and run
```
> ./configure --with-matlab-arch=glnxa64 --with-matlab=/home/magland/MATLAB/R2017a --enable-all --enable-openmp
```
where `/home/magland/MATLAB/R2017a` is the path to your matlab installation,
 and set `glnxa64` to be the appropriate architecture.
```
> make
```

If you want to do single thread as well, repeat the above exercise and
omit the following flags in the configuration:
      `--enable-all` and `--enable-openmp`
Then change that directory name to `nfft-3.3.2-single-thread`

If you use different locations, edit `addnfftpath.m` to match.


### Installing FINUFFT

See https://github.com/ahbarnett/finufft


### Installing CMCL NUFFT

This the Greengard-Lee fortran code from 2004, and is single thread only.
MEX files are shipped, so no compilation is needed.

See https://cims.nyu.edu/cmcl/nufft/nufft.html


### Installing BART NUFFT

This is a new (2016) MRI reconstruction package with a MATLAB interface.
See: https://mrirecon.github.io/bart/

The NUFFT is 3D only, and is by Frank Ong and Martin Uecker.
The MATLAB interface writes data to a tempfile, then makes a system
call to an executable `bart nufft` which reads in the data.
We have to be careful not to penalize this.

We installed v 0.4.02 from
https://github.com/mrirecon/bart/releases/tag/v0.4.02

Instructions to make are in the `README`.
On our EL7 (centos) system 
we found we had to change the following lines in their `Makefile`
in order for openblas to be found:
```
OPT = -O3 -ffast-math -I/usr/include/openblas/
BLAS_L := -L$(BLAS_BASE)/lib -llapack -lblas -lopenblas
BLAS_H := -I$(BLAS_BASE)/include -I$(BLAS_BASE)/include/openblas
```
We then did
```
make NOLAPACKE=1
```
which compiles a local copy of lapacke rather than have to install `lapacke-dev`.

We found the nonuniform point periodic box by trial and error;
it is [-N/2,N/2]^d where N is the number of modes in each dimension
(we assume equal).
We found the normalization is (sqrt(N^3)/prefac) where prefac was
found by trial and error to be the constant 1.00211.
This is strange, and probably something to do with BART's kernel.

In addition, points near box edges incorrectly wrap, giving O(1) errors.
Uecker kindly sent us a patched `grid.c` which you may find in `patches`.
Replace `bart-0.4.02/src/noncart/grid.c` with this file and recompile.

See `bart_*.m`


