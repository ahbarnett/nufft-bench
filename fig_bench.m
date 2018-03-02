% compare all nufft algs for general dim, types 1,2, general distn.
% Cleaned up. Barnett 3/2/18

clear
if 0  % compute
dim=3;
ty=1;
nudist=0; expt='unif, i7';
%nudist=4; expt='sph, i7';
multithreaded=0;
N=64; M=1e6;
o.nfftpres = 1;
nthreads = 1; if multithreaded, nthreads = java.lang.Runtime.getRuntime().availableProcessors; end
outname=sprintf('results/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',dim,ty,nudist,N,M,nthreads,o.nfftpres);
benchallcodes(ty,dim,N,M,nudist,multithreaded,outname,o)
load(outname)  % overwrites the above
else   % load a run
  %  o.nfftpres=2; benchallcodes(1,3,32,1e5,0,1,'temp',o);
  % benchallcodes(1,1,1e5,1e6,0,0,'temp',struct('nfftpres',1))
  load temp; expt='unif, i7';
  %load results/3d1_nudist4_N64_M1000000_8thr_nfftpres1; expt = 'sph, i7';
end

figure; ms=10; sy = 's';  % markersize, line symbols
legcel = {'FINUFFT','NFFT no pre'};
semilogx(rel2err(jf),run_times(jf),'ko-','markersize',ms); hold on;
semilogx(rel2err(jnp),run_times(jnp),'go-','markersize',ms);
if o.nfftpres>0
  legcel={legcel{:}, 'NFFT pre','NFFT pre (init)'};
  semilogx(rel2err(jn),run_times(jn),'gs-','markersize',ms);
  semilogx(rel2err(jn),init_times(jn),'gs:','markersize',ms);
end
if o.nfftpres>1
  legcel={legcel{:}, 'NFFT full','NFFT full (init)'};
  semilogx(rel2err(jnf),run_times(jnf),'g*-','markersize',ms);
  semilogx(rel2err(jnf),init_times(jnf),'g*:','markersize',ms);
end
if ~isempty(jc)
  legcel={legcel{:},'CMCL'};
  semilogx(rel2err(jc),run_times(jc),'bo-','markersize',ms);
end
if ~isempty(jb)
  legcel={legcel{:},'BART'};
  semilogx(rel2err(jb),run_times(jb),'co-','markersize',ms);
end
if ~isempty(jm)
  legcel={legcel{:},'Fessler pre'};
  semilogx(rel2err(jm),run_times(jm),'mo-','markersize',ms);
  semilogx(rel2err(jm),init_times(jm),'mo:','markersize',ms);
end
axis tight; v=axis; axis([1e-12 1e-1 0 v(4)])
xlabel('$\epsilon$ (relative $l_2$ error)','interpreter','latex');
ylabel('wall-clock time (s)');
legend(legcel, 'location','northeast');
plural = 's'; if nthreads==1, plural=''; end   % for text: thread(s)
title(sprintf('%dD type-%d, %s (%d thread%s): $N=%d^%d$, $M=$ %g',dim,ty,expt,nthreads,plural,N,dim,M),'interpreter','latex');


% **** questions:  why is finufft slower than CMCL, Fessler at low acc in 2D,1D?
% Why are we so slow in 2D, 1D? Indep of accuracy?
