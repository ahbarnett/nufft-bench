% run all dims comparisons and plot on same fig
% barnett 3/2/18


%multithreaded=0;ty=1;M=1e7; NN=1e6; nudist=1; expt='sph, xeon';
%multithreaded=1;ty=1;M=1e8; NN=1e7; nudist=1; expt='sph, xeon';

multithreaded=1;ty=1; M=1e8; NN=1e7; nudist=4; expt='sphquad, xeon';

% memorygraph( ....

nthreads = 1; if multithreaded, nthreads = java.lang.Runtime.getRuntime().availableProcessors; end

%NN = 1e6;  % for 1d N=2e6 we're 1.5x slower than NFFT due to fftw - why?
o = []; o.nfftpres = 1;
calc = 1;
figure;
for dim=1:3
  N=round(NN^(1/dim)); if mod(N,2)==1, N=N+1; end   % modes per dim, even
  outname=sprintf('results/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',dim,ty,nudist,N,M,nthreads,o.nfftpres);
  if calc, benchallcodes(ty,dim,N,M,nudist,multithreaded,outname,o); end
  subplot(1,3,dim); plotabench(outname,expt); drawnow
end


% to fix: use memroygraph. make in benchallcodes, running... show time too.

ty = 2;

figure;
for dim=1:3
  N=round(NN^(1/dim)); if mod(N,2)==1, N=N+1; end   % modes per dim, even
  outname=sprintf('results/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',dim,ty,nudist,N,M,nthreads,o.nfftpres);
  if calc, benchallcodes(ty,dim,N,M,nudist,multithreaded,outname,o); end
  subplot(1,3,dim); plotabench(outname,expt); drawnow
end



% single ============================
multithreaded=0;ty=1; M=1e7; NN=1e6; nudist=4; expt='sphquad, xeon';


nthreads = 1; if multithreaded, nthreads = java.lang.Runtime.getRuntime().availableProcessors; end

figure;
for dim=1:3
  N=round(NN^(1/dim)); if mod(N,2)==1, N=N+1; end   % modes per dim, even
  outname=sprintf('results/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',dim,ty,nudist,N,M,nthreads,o.nfftpres);
  if calc, benchallcodes(ty,dim,N,M,nudist,multithreaded,outname,o); end
  subplot(1,3,dim); plotabench(outname,expt); drawnow
end

ty = 2;
figure;
for dim=1:3
  N=round(NN^(1/dim)); if mod(N,2)==1, N=N+1; end   % modes per dim, even
  outname=sprintf('results/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',dim,ty,nudist,N,M,nthreads,o.nfftpres);
  if calc, benchallcodes(ty,dim,N,M,nudist,multithreaded,outname,o); end
  subplot(1,3,dim); plotabench(outname,expt); drawnow
end

% then add the insets w/ nudist

