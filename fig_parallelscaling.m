% parallel scaling of finufft w/ # threads.
% Barnett 3/19/18

addpath ../finufft/matlab/
clear

dim=3; nudist=4;

N1 = 100; N2=N1;N3=N1; % # modes in each dim. fixed.

o = [];  % finufft opts
o.debug = 0;
o.fftw = 0;      % default
eps = 1e-6;


if 1 % STRONG SCALING ---------------------------------------------

M=1e7;
[x y z nudnam] = nudata(dim,nudist,M); M=numel(x);   % NU locs, any dim
c = randn(M,1)+1i*randn(M,1);

N1 = 100; N2=N1;N3=N1; % # modes in each dim

o = [];  % finufft opts
o.debug = 0;
o.fftw = 0;      % default
eps = 1e-6;

nthrs = [1 2 3 4 8 12 16 24];
nthrold = maxNumCompThreads;

ts = nan(numel(nthrs),2);
for i=1:numel(nthrs)   % run the expt
  nthr = nthrs(i);
  maxNumCompThreads(nthr);
  t = tic; [f ier]=finufft3d1(x,y,z,c,+1,eps,N1,N2,N3,o);
  ts(i,1) = toc(t);
  t = tic; [c ier]=finufft3d2(x,y,z,+1,eps,f,o);
  ts(i,2) = toc(t);
  fprintf('%d threads:\t 3d1 %.3g s \t 3d2 %.3g s\n',nthr,ts(i,1),ts(i,2))
end
maxNumCompThreads(nthrold);
disp('        nthr       time(s)')
disp([nthrs', ts])
clear x y z c f; save results/xeon/strongscaling_3d_M1e7_N100

figure;
bar(nthrs,M./nthrs'./ts);
title('\hspace{.7in}strong scaling, $\epsilon=10^{-6}$, $M=10^7$, sph quad','interpreter','latex');
xlabel('p (number of threads)'); ylabel('throughput (NU pts / thread / sec)');
legend('3D type-1', '3D type-2','location','northeast');
set(gcf,'paperposition',[0 0 4 3]);
print -depsc2 results/xeon/strongscaling_3d_M1e7_N100.eps
end


if 1 % WEAK SCALING ---------------------------------------------

nthrs = [1 2 3 4 8 12 16 24];
nthrold = maxNumCompThreads;
nthr_avail = java.lang.Runtime.getRuntime().availableProcessors;

ts = nan(numel(nthrs),2);
for i=1:numel(nthrs)   % run the expt
  nthr = nthrs(i);
  M=1e7 * nthr;
  maxNumCompThreads(nthr_avail);
  [x y z nudnam] = nudata(dim,nudist,M); M=numel(x);   % NU locs, any dim
  c = randn(M,1)+1i*randn(M,1);  
  maxNumCompThreads(nthr);
  t = tic; [f ier]=finufft3d1(x,y,z,c,+1,eps,N1,N2,N3,o);
  ts(i,1) = toc(t);
  t = tic; [c ier]=finufft3d2(x,y,z,+1,eps,f,o);
  ts(i,2) = toc(t);
  fprintf('%d threads:\t 3d1 %.3g s \t 3d2 %.3g s\n',nthr,ts(i,1),ts(i,2))
end
maxNumCompThreads(nthrold);
disp('        nthr       3d1 time(s)    3d2 time(s)')
disp([nthrs', ts])
clear x y z c f; save results/xeon/weakscaling_3d_M1e7p_N100

figure;
h = bar(nthrs,ts);
title('weak scaling, $\epsilon=10^{-6}$, \# NU pts $M=10^7 p$, sph quad','interpreter','latex');
xlabel('p (number of threads)'); ylabel('time (s)');
legend('3D type-1', '3D type-2','location','northwest');
set(gcf,'paperposition',[0 0 4 3]);
print -depsc2 results/xeon/weakscaling_3d_M1e7p_N100.eps

end
