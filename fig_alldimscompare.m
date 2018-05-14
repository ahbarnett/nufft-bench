% playing with figures for comparison. Barnett 3/13/18
% based on fig_icerm_talk.m

% laptop warm-up...
%cpu='i7'; dim=2; ty=1; o.nfftpres=2; nudists = 0; nthrs = 8; Ns = 1e3; Ms = 1e7; 

% Tues dinnertime xeon tests...
%cpu='xeon'; dim=1; ty=1; o.nfftpres=2; nudists=4; nthrs=24; Ns=1e7; Ms=1e8;
%cpu='xeon'; dim=2; ty=1; o.nfftpres=2; nudists=4; nthrs=24; Ns=3162; Ms=1e8;

if 0 % WEd am...     all 3D figs:
cpu='xeon'; dim=3;
s = [1 1 1 1]; % use for const lists
nudists=[0 0 4 4 0 0 4 4]; tys = [1 2 1 2 1 2 1 2];
nthrs=[s 24*s]; NNs=[100*s 216*s]; Ms=[1e7*s 1e8*s];
end

if 1 % all 2d
cpu='xeon'; dim=2;
s = [1 1 1 1]; % use for const lists
nudists=[0 0 4 4 0 0 4 4]; tys = [1 2 1 2 1 2 1 2];
nthrs=[s 24*s]; NNs=[1000*s 3162*s]; Ms=[1e7*s 1e8*s];
end

if 0 % 1d
cpu='xeon'; dim=1;
s = [1 1]; % use for const lists
nudists=[0 0 0 0]; tys = [1 2 1 2];
nthrs=[s 24*s]; NNs=[1e6*s 1e7*s]; Ms=[1e7*s 1e8*s];
end

if 0 % 1d 1thr
cpu='xeon'; dim=1;
s = [1 1]; % use for const lists
nudists=[0 0]; tys = [1 2];
nthrs=s; NNs=[1e6*s]; Ms=[1e7*s];
end
  

for n = 1:numel(Ms)                        % ....... loop over cases
  nudist = nudists(n); nthreads = nthrs(n); N = NNs(n); M = Ms(n); ty=tys(n);
  o = []; o.nfftpres = 2; if dim==3 & M>1e7, o.nfftpres=1; end
  data = sprintf('results/%s/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',cpu,dim,ty,nudist,N,M,nthreads,o.nfftpres)
  load(data)
  if dim==2 & ty==1 & nthreads==24 & nudist==0, jf=jf(2:end); end  % kill glitch
  
  figure; %subplot(2,2,n);
  ms=10;
  legcel = {'FINUFFT','NFFT'};
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
  if mod(n,4)==1
    legend(legcel, 'location','best');
  end
  
  if nudist==0     % decide NU pt dist names appearing on plots
    if dim==1, nudnam='rand', elseif dim==2, nudnam='rand square'; else, nudnam='rand cube'; end
  elseif nudist==4
    if dim==2, nudnam='disc quad'; else nudnam='sph quad'; end
  end
  
  plural = 's'; if nthreads==1, plural=''; end   % for text: thread(s)
  tstr = sprintf('type-%d, %d thread%s: $N=%d^%d$, $M=$ %.2g, %s',ty,nthreads,plural,N,dim,M,nudnam);
  title(tstr,'interpreter','latex');

  % inset...
  if dim==3
  pos = [.67 .67]; if n==1, pos = [.13 .59]; elseif n==5, pos(2)= .45; end
  axes('position',[pos .25 .25]);         % inset (comes after title!)
  elseif dim==2
    pos = [.17 .55]; if nthreads==1, pos = [.17 .17]; end
    axes('position',[pos .17 .17]);         % inset (comes after title!)
  end
  if dim>1, [x y z nam] = nudata(dim,nudist,2500); end
  if dim==2
    plot(x,y,'.','markersize',.5);
  elseif dim==3
    plot3(x,y,z,'.','markersize',.5); axis vis3d; view(20,10);
  end
  if dim>1, axis equal tight off; end
  
  set(gcf,'paperposition', [0 0 4 3.5]);
  print('-depsc2',[data '.eps'])

  if 0 & o.memcpu
    memorygraph('plot',bytes,est_times,cpu_times,cpu_usages,labelstrings,labeltimes);
    a=axes; t=title(tstr,'interpreter','latex'); a.Visible='off'; t.Visible='on';
  end
end                                             % ......

