# nufft-bench

Codes to generate figures and tables
comparing performance of various NUFFT codes,
for the FINUFFT paper.
Written in MATLAB/octave.

Alex Barnett and Jeremy Magland; contributions by Joakim Anden.

January 2018.




## Installing NUFFT packages

### Installing NFFT

Download current source from here:

https://www-user.tu-chemnitz.de/~potts/nfft/download.php

Extract to directory `nfft-3.3.2` which should be placed in the
same directory as `nufft-bench`.

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

This is single thread only.

See https://cims.nyu.edu/cmcl/nufft/nufft.html

