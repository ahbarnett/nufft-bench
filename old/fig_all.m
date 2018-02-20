% Top-level driver script for NUFFT benchmarking/
% Barnett 2/20/18

if 0
clear
calc=1;
expt='xeon'; multi=1;
nfftpres = 1;  % note: 0 no extra RAM, 1 needs 0.6kB/NUpt, 2 needs 35kB/NUpt!
               % (these are for m=6 kernel)
%N=64; M=1e7;
N=128; M=1e8;
fig_3d(expt,1,multi,N,M,nfftpres,calc)
fig_3d(expt,2,multi,N,M,nfftpres,calc)
end
