% compare all nufft algs for 3d type 1, unif rand data.
% Attempt to clean up. Barnett 2/20/18

clear
dim=3;
ty=1;
nudist=0;
multi=1;
N=64; M=1e6;
o.nfftpres = 2;
nth = 1;
if multi, nth = java.lang.Runtime.getRuntime().availableProcessors; end
outname=sprintf('results/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',dim,ty,nudist,N,M,nth,o.nfftpres);
benchallcodes(ty,dim,N,M,nudist,multi,outname,o)
load(fnam)  % overwrites the above
