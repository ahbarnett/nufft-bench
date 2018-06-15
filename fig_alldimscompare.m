% Main comparison figures. Barnett 3/13/18, cleaned up 5/18/18.
% based on fig_icerm_talk.m

% laptop warm-up...
%cpu='i7'; dim=2; ty=1; o.nfftpres=2; nudists = 0; nthrs = 8; Ns=1e3; Ms=1e7; 
% xeon warm-up...
%cpu='xeon'; dim=1; ty=1; o.nfftpres=2; nudists=4; nthrs=24; Ns=1e7; Ms=1e8;
%cpu='xeon'; dim=2; ty=1; o.nfftpres=2; nudists=4; nthrs=24; Ns=3162; Ms=1e8;

cpu='xeon'; mt=24; % max # threads
dim=3;   % choose 1,2, or 3.

s = [1 1]; % use for const lists
Ms=kron(s,[1e7 1e8]); % same for 1d,2d
if dim==1
  nudists=[0*s 0*s]; tys = [1*s 2*s];
  nthrs=kron(s,[1 mt]); NNs=kron(s,[1e6 1e7]); 
elseif dim==2
  nudists=[4*s 4*s]; tys = [1*s 2*s];
  nthrs=kron(s,[1 mt]); NNs=kron(s,[1000 3162]); 
elseif dim==3
  s = [1 1 1 1]; % use for const lists. Note all 4 1-thr then 4 mt-thr
  nudists=kron(s,[0 4]); tys = [1 1 2 2 1 1 2 2]; 
  nthrs=[1*s mt*s]; NNs=[100*s 216*s];
  Ms=[1e7*s 1e8*s];
end

subfigs='abcdabcd';   % covers the 2 figs for 3d case
for n = 1:numel(Ms)                        % ....... loop over cases
  nudist = nudists(n); nthreads = nthrs(n); N = NNs(n); M = Ms(n); ty=tys(n);
  o = []; o.nfftpres = 2; if dim==3 & M>1e7, o.nfftpres=1; end
  data = sprintf('results/%s/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',cpu,dim,ty,nudist,N,M,nthreads,o.nfftpres)
  load(data)
  %if dim==2 & ty==1 & nthreads==24 & nudist==0, jf=jf(2:end); end  % kill glitch
  
  figure;
  ms=10;  % markersize
  legcel = {'FINUFFT','NFFT'};   % start the legend string list
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
    legcel={legcel{:},'MIRT pre'};
    semilogx(rel2err(jm),run_times(jm),'bd-','markersize',ms);
    %legcel={legcel{:},'MIRT pre (init)'};  % very long! (100x longer)
    %semilogx(rel2err(jm),init_times(jm),'mo:','markersize',ms);
  end
  set(gca,'ysc','log'); axis tight; v=axis; axis([6e-13 1e-1 v(3:4)]);
  text(2e-12,v(4)^0.9*v(3)^0.1,['(' subfigs(n) ')'],'fontsize',12);
  xlabel('$\epsilon$ (relative $\ell_2$ error)','interpreter','latex');
  ylabel('run time (s)');
  if mod(n,4)==1
    if dim==1,  v=axis; v(3)=0.8*v(3); axis(v);
      h=legend(legcel, 'location','southwest');
      v=get(h,'position'), v(1:2)=v(1:2)+[.03,.02]; set(h,'position',v);
    elseif dim==2, h=legend(legcel, 'location','north');
      v=get(h,'position'); v(1)=v(1)+0.05; set(h,'position',v);
    else, h=legend(legcel, 'location','best');
    end
  end
    
  if nudist==0     % decide NU pt dist names appearing on plots
    if dim==1, nudnam='rand', elseif dim==2, nudnam='rand square'; else, nudnam='rand cube'; end
  elseif nudist==4
    if dim==2, nudnam='disc quad'; else nudnam='sph quad'; end
  end
  
  plural = 's'; if nthreads==1, plural=''; end   % for text: thread(s)
  tstr = sprintf('type-%d, %d thread%s: $N=%d^%d$, $M=$ %.2g, %s',ty,nthreads,plural,N,dim,M,nudnam);
  title(tstr,'interpreter','latex');

  % inset...  (comes after title!)
  if dim==3
  pos = [.67 .67]; if n==1, pos = [.13 .59]; elseif n==5, pos(2)= .45; end
  axes('position',[pos .25 .25]);
  elseif dim==2
    pos = [.71 .45]; axes('position',[pos .17 .17]);
  end
  if dim>1, [x y z nam] = nudata(dim,nudist,2500); end
  if dim==2
    plot(x,y,'.','markersize',.5);
  elseif dim==3
    plot3(x,y,z,'.','markersize',.5); axis vis3d; view(20,10);
  end
  if dim>1, axis equal tight off; end
  
  set(gcf,'paperposition', [0 0 4.1 3.6]);  % control overall font size!
  print('-depsc2',[data '.eps'])

  if 0 & o.memcpu
    memorygraph('plot',bytes,est_times,cpu_times,cpu_usages,labelstrings,labeltimes);
    a=axes; t=title(tstr,'interpreter','latex'); a.Visible='off'; t.Visible='on';
  end
end                                             % ......

