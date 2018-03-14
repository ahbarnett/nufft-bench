% playing with figures for comparison. Barnett 3/13/18
% based on fig_icerm_talk.m

% multithr
cpu='i7';
dim=2; ty=1; o.nfftpres=2;

nudists = 0;
nthrs = 8;
Ns = 1e3;
Ms = 1e7; 

for n = 1%[1 2]                           % ....... loop over cases
  nudist = nudists(n); nthreads = nthrs(n); N = Ns(n); M = Ms(n);
  data = sprintf('results/%s/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',cpu,dim,ty,nudist,N,M,nthreads,o.nfftpres);
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
  if o.nfftpres>1
    legcel={legcel{:}, 'NFFT full pre'}; %'init NFFT pre'};
    semilogx(rel2err(jnf),run_times(jnf),'r*-','markersize',5);
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
  tstr = sprintf('%dD type-%d, %s (%d thread%s): $N=%d^%d$, $M=$ %.2g, %s',dim,ty,cpu,nthreads,plural,N,dim,M,nudnam);
  title(tstr,'interpreter','latex');

  axes('position',[.7 .5+(n-1)*.07 .2 .2]);         % inset (comes after title!)
  [x y z nam] = nudata(dim,nudist,5e3);
  if dim==2
    plot(x,y,'.','markersize',1);
  elseif dim==3
    plot3(x,y,z,'.','markersize',1); axis vis3d; view(20,10);
  end
  axis equal tight off
  
  %set(gcf,'paperposition', [0 0 5 5]);
  %print('-depsc2',[data '.eps'])

  memorygraph('plot',bytes,est_times,cpu_times,cpu_usages,labelstrings,labeltimes);
  a=axes; t=title(tstr,'interpreter','latex'); a.Visible='off'; t.Visible='on';
end                                             % ......

