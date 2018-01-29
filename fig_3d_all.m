% Top-level driver script for 3D expts for NUFFT.
% Barnett 1/28/18

clear
calc=1;
expt='i7'; multi=0;
nfftpres = 1;  % note: 0 no extra RAM, 1 needs 0.6kB/NUpt, 2 needs 35kB/NUpt!
               % (these are for m=6 kernel)
N=32; M=1e6;
fig_3d(expt,1,multi,N,M,nfftpres,calc)
fig_3d(expt,2,multi,N,M,nfftpres,calc)
