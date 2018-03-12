# nufft-bench

Codes to generate figures and tables
comparing performance of various NUFFT codes.
Written in MATLAB/octave.

Alex Barnett and Jeremy Magland; contributions by Joakim Anden.

January-March 2018.

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

Note: we find success with gcc 4.8.5 and 6.4.0. We find compilation
fails for gcc 7.2.0 and 7.3.0 with:

```
cc1: internal compiler error: in gimplify_modify_expr, at gimplify.c:5638
```

This gcc bug seems to be slated to be fixed in 7.4


### Installing FINUFFT

See https://github.com/ahbarnett/finufft

This is multithreaded.

### Installing CMCL NUFFT

This the Greengard-Lee fortran code from 2009-2014, and is single-thread only.
MEX files are shipped, so no compilation is needed.
We use nufft-all version 1.3.3.

See https://cims.nyu.edu/cmcl/nufft/nufft.html


### Installing BART NUFFT

This is a new (2016) MRI reconstruction package with a MATLAB interface.
See: https://mrirecon.github.io/bart/

The NUFFT is 3D only, and is by Frank Ong and Martin Uecker.
It uses all available threads.
The MATLAB interface writes data to a tempfile, then makes a system
call to an executable `bart nufft` which reads in the data.
We have to be careful to use problem sizes that don't
penalize this.

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
This choice improves accuracy from 3 to 5 digits.
This is strange, and probably derives from BART's spreading kernel.

In addition, points near box edges incorrectly wrap, giving O(1) errors.
Uecker kindly sent us a patched `grid.c` which you may find in `patches`.
Replace `bart-0.4.02/src/noncart/grid.c` with this file and recompile.

See `bart_*.m`


### Installing MIRT (Fessler) NUFFT

The complete package is at:
https://web.eecs.umich.edu/~fessler/irt/fessler.tgz

All MEX files are included; no compilation is needed.
The precomputation is very expensive (up to 100x slower than
the apply when limited to one thread).
The apply appears to be single-threaded and very fast.

