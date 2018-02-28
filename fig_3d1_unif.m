% compare all nufft algs for 3d type 1, unif rand data.
% Attempt to clean up. Barnett 2/20/18

clear
dim=3;
ty=1;
%nudist=0; expt='unif, i7';
nudist=4; expt='sph, i7';
multi=1;
N=64; M=1e6;
o.nfftpres = 1;
nth = 1; if multi, nth = java.lang.Runtime.getRuntime().availableProcessors; end
outname=sprintf('results/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',dim,ty,nudist,N,M,nth,o.nfftpres);
benchallcodes(ty,dim,N,M,nudist,multi,outname,o)
load(outname)  % overwrites the above
 
%  o.nfftpres=2; benchallcodes(1,3,32,1e5,0,1,'temp',o);
figure; ms=10; sy = 's';  % markersize, line symbols
legcel = {'FINUFFT','NFFT no pre'};
semilogx(rel2err(jf),run_times(jf),'ko-','markersize',ms); hold on;
semilogx(rel2err(jnp),run_times(jnp),'go-','markersize',ms);
if o.nfftpres>0
  legcel={legcel{:}, 'NFFT pre','precomp NFFT'};
  semilogx(rel2err(jn),run_times(jn),'gs-','markersize',ms);
  semilogx(rel2err(jn),init_times(jn),'gs:','markersize',ms);
end
if o.nfftpres>1
  legcel={legcel{:}, 'NFFT full pre','precomp NFFT full'};
  semilogx(rel2err(jnf),run_times(jnf),'g*-','markersize',ms);
  semilogx(rel2err(jnf),init_times(jnf),'g*:','markersize',ms);
end
legcel={legcel{:},'BART','Fessler'};
semilogx(rel2err(jb),run_times(jb),'bo-','markersize',ms);
semilogx(rel2err(jm),run_times(jm),'mo-','markersize',ms);
% insert CMCL ***
axis tight; v=axis; axis([1e-12 1e-1 0 v(4)])
xlabel('$\epsilon$ (relative $l_2$ error)','interpreter','latex');
ylabel('wall-clock time (s)');
legend(legcel, 'location','northeast');
%set(gca,'yscale','log')
title(sprintf('%dD type-%d, %s (%d threads): $N=%d^%d$, $M=$ %g',dim,ty,expt,nth,N,dim,M),'interpreter','latex');
