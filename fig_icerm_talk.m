% icerm talk new figs, from existing xeon data. Barnett 3/9/18
% based on plotabench.m

% multithr
expt='xeon';
dim=3; ty=1; o.nfftpres=1;

nudists = [0 4];    % our various cases to plot
nthrs = [1 24];
Ns = [100 216];   % should really get these automatically from NN
Ms = [1e7 1e8];

for n = [1 2]                           % ....... loop over cases
  nudist = nudists(n); nthreads = nthrs(n); N = Ns(n); M = Ms(n);
  data = sprintf('results/%s/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',expt,dim,ty,nudist,N,M,nthreads,o.nfftpres);
  load(data)

  figure; ms=10;
  legcel = {'FINUFFT','NFFT no pre'};
  semilogx(rel2err(jf),run_times(jf),'k.-','markersize',ms); hold on;
  semilogx(rel2err(jnp),run_times(jnp),'ms-','markersize',ms);
  if o.nfftpres>0
    legcel={legcel{:}, 'NFFT pre'}; %'init NFFT pre'};
    semilogx(rel2err(jn),run_times(jn),'ro-','markersize',5);
    %semilogx(rel2err(jn),init_times(jn),'r:');
  end
  if ~isempty(jc)
    legcel={legcel{:},'CMCL'};
    semilogx(rel2err(jc),run_times(jc),'g.-','markersize',ms);
  end
  if ~isempty(jb)
    legcel={legcel{:},'BART'};
    semilogx(rel2err(jb),run_times(jb),'c*-','markersize',ms);
  end
  if ~isempty(jm)
    legcel={legcel{:},'Fessler pre'};
    semilogx(rel2err(jm),run_times(jm),'bd-','markersize',ms);
    %legcel={legcel{:},'Fessler pre (init)'};
    %semilogx(rel2err(jm),init_times(jm),'mo:','markersize',ms);
  end
  set(gca,'ysc','log');
  axis tight; v=axis; axis([6e-13 1e-1 0 v(4)]);
  xlabel('$\epsilon$ (relative $l_2$ error)','interpreter','latex');
  ylabel('wall-clock time (s)');
  legend(legcel, 'location','northeast');
  nudnam
  
  plural = 's'; if nthreads==1, plural=''; end   % for text: thread(s)
  title(sprintf('%dD type-%d, %s (%d thread%s): $N=%d^%d$, $M=$ %.2g, %s',dim,ty,expt,nthreads,plural,N,dim,M,nudnam),'interpreter','latex');

  axes('position',[.7 .5+(n-1)*.07 .2 .2]);         % inset (comes after title!)
  [x y z nam] = nudata(dim,nudist,5e3);
  plot3(x,y,z,'.','markersize',1);
  axis vis3d equal tight off; view(20,10);

  set(gcf,'paperposition', [0 0 5 5]);
  print('-depsc2',[data '_icerm.eps'])
end                                             % ......
