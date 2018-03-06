function plotabench(file,expt)
% func to plot given benchmark vars in current workspace
%
% plotabench(file,expt) file is data filename, expt is a text string.

jf=0;   % int needed to avoid the IRT func jf(); see https://www.mathworks.com/help/matlab/import_export/troubleshooting-loading-variables-within-a-function.html
load(file)

%figure;
ms=10; sy = 's';  % markersize, line symbols
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
  %legcel={legcel{:},'Fessler pre (init)'};
  %semilogx(rel2err(jm),init_times(jm),'mo:','markersize',ms);
end
axis tight; v=axis; axis([1e-12 1e-1 0 v(4)])
xlabel('$\epsilon$ (relative $l_2$ error)','interpreter','latex');
ylabel('wall-clock time (s)');
legend(legcel, 'location','northeast');
plural = 's'; if nthreads==1, plural=''; end   % for text: thread(s)
title(sprintf('%dD type-%d, %s (%d thread%s): $N=%d^%d$, $M=$ %g',dim,ty,expt,nthreads,plural,N,dim,M),'interpreter','latex');
