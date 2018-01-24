% 3d1 multicore comparisons figs
% Barnett 11/7/17

clear
cpuname = 'i7'; multi=1; N=64; M=1e6;
if multi, thrstr='multi'; else thrstr='single'; end
fighead = sprintf('%s_%s_3d1_N%d_M%d',cpuname,thrstr,N,M);  % figure filehead
nthreads = java.lang.Runtime.getRuntime().availableProcessors;  % for title
if ~multi, nthreads=0; end
ds=0:1;  % nudist vals
nams = {'unif','sph'};  % nudist names for fig
ps=0:1;  % nfft pre vals
calc = 0;  % run expts (saves to files)
if calc
  for i = 1:numel(ds), nudist=ds(i);
    for j = 1:numel(ps), nfftpre = ps(j);
      head{i}{j} = run_perf_3d1(N,M,nudist,multi,nfftpre,cpuname);
    end
  end
  save 3d1_heads head
else, load 3d1_heads
end

% make plots
for i = 1:numel(ds), nudist=ds(i);
  np = load(head{i}{1}); % nfft no pre
  load(head{i}{2});      % nfft pre
  load 3d1_heads  % *** kill
  figure;
  semilogx(errors(ii1),run_times(ii1),'b.-'); hold on; % finufft
  semilogx(errors(ii2),run_times(ii2),'r.-');      % nfft pre_psi
  semilogx(errors(ii2),init_times(ii2),'r.:');    % nfft's precomp time
  semilogx(np.errors(np.ii2),np.run_times(np.ii2),'m.-'); % nfft no pre
  axis tight; v=axis; axis([1e-12 1e-1 0 v(4)])
  xlabel('relative l_2 error \epsilon');
  ylabel('wall-clock time (s)');
  legend({'finufft','NFFT','NFFT init','NFFT no pre'});
  title(sprintf('3D type-1, N=%d^3, M=%g (%s), %d threads',NN,M,nams{i},nthreads)); % ** edit
  set(gca,'fontsize',22); drawnow
  set(gcf,'paperposition',[0 0 4 4]);
  print([fighead '_' nams{i} '.eps'], '-depsc2');
end
