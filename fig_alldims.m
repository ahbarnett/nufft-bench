% run all dims comparisons and plot on same fig
% barnett 3/2/18

multithreaded=1;
nthreads = 1; if multithreaded, nthreads = java.lang.Runtime.getRuntime().availableProcessors; end

ty=1;
nudist=1; expt='sph, i7';  % since nudist=4 unavail for dim<3
M = 1e6;
NN = 1e6;  % for 1d N=2e6 we're 1.5x slower than NFFT due to fftw - why?
o = []; o.nfftpres = 1;
calc = 1;
figure;
for dim=1:3
  N=round(NN^(1/dim))   % modes per dim
  outname=sprintf('results/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',dim,ty,nudist,N,M,nthreads,o.nfftpres);
  if calc, benchallcodes(ty,dim,N,M,nudist,multithreaded,outname,o); end
  subplot(1,3,dim); plotabench(outname,expt); drawnow
end
