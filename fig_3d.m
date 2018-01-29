function fig_3d(expt,ty,multi,N,M,nfftpres,calc)
% FIG_3D   Driver to make 3d comparisons fig for one type, dataset size & cpu
%
% fig_3d(expt,ty,multi,N,M) runs a set of expts for a given dataset size & cpu
%  and makes a plot. Two types of NU pt dists are tested.
%
%  expt = string giving cpu/expt name
%  ty = 1 or 2 (type of xform)
%  multi = 0 (single thread), 1 (multi-thread)
%  N = # modes in each direction (cube)
%  M = # NU pts
%  nfftpres = which NFFT PSI precompute stage to test up to (PRE_PHI always on).
%             0: no pre psi,  1: also test pre psi,  2: also test pre full psi.
%                               (needs ~ 600B/NUpt)    (needs ~ 35kB/NUpt)
%  calc = 0 (load precomputed .mat files, just replot fig), 1 (re-calc & plot)
%
% Writes multiple .mat and a single .eps to results/ directory
%
% Barnett 1/28/18

if nargin<1, expt = 'generic'; end
if nargin<2, ty=1; end
if nargin<3, multi=1; end
if nargin<4, N=64; end
if nargin<5, M=1e6; end
if nargin<6, nfftpres = 1; end
if nargin<7, calc=1; end

dim = 3;  % *** for now

% build the fig filehead...
nthreads = 1;     % used for title, not calc
if multi, nthreads = java.lang.Runtime.getRuntime().availableProcessors; end
fighead = sprintf('results/%s_%dthread_%dd%d_N%d_M%d_np%d',expt,nthreads,dim,ty,N,M,nfftpres);

% possibly do the calc...
nudists=0:1;  % nudist (types of NU distn) vals to loop over
if calc
  for i = 1:numel(nudists), nudist=nudists(i)
    head{i} = benchallcodes(ty,dim,N,M,nudist,multi,nfftpres,expt);
  end
  save([fighead '_heads'],'head');
else
  load([fighead '_heads']);        % no calc: get filenames into head{:}
end

% make combined plots with both distn's, good legend (hence plotting) order...
figure; ms=10; symbs = 's*';  % markersize, line symbols for nudists
legcel = {'FINUFFT, unif','FINUFFT, sph','NFFT no pre, unif','NFFT no pre, sph'};   % could use nudnam here if clever
for i = 1:numel(nudists), sy=symbs(i);      % FINUFFT
  load(head{i});        % crucial this doesn't wipe out variables i, j, etc!
  semilogx(rel2err(jf),run_times(jf),['b' sy '-'],'markersize',ms); hold on;
  plotcmcl = ~isempty(jc);
end
for i = 1:numel(nudists), sy=symbs(i);           % NFFT no pre
  load(head{i});
  semilogx(rel2err(jn),run_times(jn),['g' sy '-'],'markersize',ms);
end
if nfftpres>0                                  % append NFFT pre psi
  legcel = {legcel{:}, 'NFFT, precomp', 'NFFT pre, unif','NFFT pre, sph'};
  semilogx(rel2err(jnp),init_times(jnp),'r.:','markersize',ms);  % init
  for i = 1:numel(nudists), sy=symbs(i);
    load(head{i});
    semilogx(rel2err(jnp),run_times(jnp),['r' sy '-'],'markersize',ms);
  end
end
if plotcmcl
  legcel = {legcel{:}, 'CMCL, unif','CMCL, sph'};
  for i = 1:numel(nudists), sy=symbs(i);
    load(head{i});
    semilogx(rel2err(jc),run_times(jc),['c' sy '-'],'markersize',ms);
  end
end
axis tight; v=axis; axis([1e-12 1e-1 0 v(4)])
xlabel('$\epsilon$ (relative $l_2$ error)','interpreter','latex');
ylabel('wall-clock time (s)');
legend(legcel, 'location','northeast');
title(sprintf('3D type-%d, %s (%d threads): $N=%d^3$, $M=$ %g',ty,expt,nthreads,N,M),'interpreter','latex');
drawnow
set(gcf,'paperposition',[0 0 5 5]);    % output fig...
print([fighead '.eps'], '-depsc2');
