% driver to make 3d1 multicore comparisons figs
% Barnett 11/7/17

clear
expt = 'i7'; multi=1; N=64; M=1e6;
if multi, thrstr='multi'; else thrstr='single'; end
fighead = sprintf('results/%s_%s_3d1_N%d_M%d',expt,thrstr,N,M);  % fig filehead
nthreads = java.lang.Runtime.getRuntime().availableProcessors;  % for title
if ~multi, nthreads=0; end
ds=0:1;  % nudist vals
nams = {'unif','sph'};  % nudist names for fig
ps=0:1;  % nfft pre vals
calc = 0;  % 0: use files, 1: run expts (saves to files)
if calc
  for i = 1:numel(ds), nudist=ds(i);
    for j = 1:numel(ps), nfftpre = ps(j);
      head{i}{j} = run_perf_3d1(N,M,nudist,multi,nfftpre,expt);
    end
  end
  save results/3d1_heads head
else, load results/3d1_heads; end  % just the filenames

if 1 % make combined plots with both distn's, good plotting order
  figure; symbs = 's*';  % line symbols for unif, sph.
  ms = 10;   % markersize
  np = load(head{1}{1}); load(head{1}{2});  % nfft no pre, for unif (sph same)
  for i = 1:numel(ds), nudist=ds(i); sy=symbs(i);
    np = load(head{i}{1}); load(head{i}{2});  % nfft no pre & pre
    semilogx(errors(ii1),run_times(ii1),['b' sy '-'],'markersize',ms); hold on; % finufft
  end
  semilogx(errors(ii2),init_times(ii2),['r' symbs(1) ':'],'markersize',ms);   % nfft precomp
  for i = 1:numel(ds), nudist=ds(i); sy=symbs(i);
    np = load(head{i}{1}); load(head{i}{2});  % nfft no pre & pre
    semilogx(errors(ii2),run_times(ii2),['r' sy '-'],'markersize',ms);      % nfft pre_psi    
  end
  for i = 1:numel(ds), nudist=ds(i); sy=symbs(i);
    np = load(head{i}{1}); load(head{i}{2});  % nfft no pre & pre
    semilogx(np.errors(np.ii2),np.run_times(np.ii2),['g' sy '-'],'markersize',ms); % nfft no pre
  end
  axis tight; v=axis; axis([1e-12 1e-1 0 v(4)])
  xlabel('$\epsilon$ (relative $l_2$ error)','interpreter','latex');
  ylabel('wall-clock time (s)');
  legend({'FINUFFT, unif','FINUFFT, sph','NFFT init','NFFT, unif', 'NFFT, sph','NFFT no pre, unif','NFFT no pre, sph'},'location','northeast');
  title(sprintf('3D type-1, %s (%d threads): $N=%d^3$, $M=$%g',expt,nthreads,NN,M),'interpreter','latex');
else   % make plots: separate plots for each dist version
  for i = 1:numel(ds), nudist=ds(i);
    np = load(head{i}{1}); % nfft no pre
    load(head{i}{2});      % nfft pre
    figure;
    semilogx(errors(ii1),run_times(ii1),'b.-'); hold on; % finufft
    semilogx(errors(ii2),run_times(ii2),'r.-');      % nfft pre_psi
    semilogx(errors(ii2),init_times(ii2),'r.:');    % nfft's precomp time
    semilogx(np.errors(np.ii2),np.run_times(np.ii2),'m.-'); % nfft no pre
    axis tight; v=axis; axis([1e-12 1e-1 0 v(4)])
    xlabel('relative l_2 error \epsilon');
    ylabel('wall-clock time (s)');
    legend({'finufft','NFFT','NFFT init','NFFT no pre'});
    title(sprintf('3D type-1, N=%d^3, M=%g (%s), %d threads',NN,M,nams{i},nthreads));
  end
end
set(gcf,'paperposition',[0 0 5 5]);
print([fighead '.eps'], '-depsc2');
