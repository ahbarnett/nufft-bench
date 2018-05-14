% 2d type1 compare algs: finufft vs NFFT. Simple, with M sweep.
% User needs to set paths below.
% Barnett 5/4/18
clear

maxthr = java.lang.Runtime.getRuntime().availableProcessors;
nthr = maxthr; %1;
maxNumCompThreads(nthr);

tol = 1e-6;   % note 6 is even, so NFFT m hits it better than if odd.

n=256;    % modes per dim
N=n^2;
densities = 10.^(-3:2);   % set of ratios M/N

speedups = nan*densities;
for i=1:numel(densities), dens = densities(i);   % sweep over M .........
  M = round(dens*N); % choose # NU
  
  reps = ceil(1e7/(N+M));  % # trials for ~1 sec runtime
  fprintf('M/N=%.3g: M=%d, reps=%d...\n',dens,M,reps)
  
  x=pi*(2*rand(M,1)-1); y=pi*(2*rand(M,1)-1);   % NU pts in [-pi,pi]^2
  c = randn(M,1)+1i*randn(M,1);   % strengths (reused for now)

  addpath ~/numerics/finufft/matlab
  o.fftw=1;   % FFTW_MEASURE
  isign=+1;  % matches NFFT (-1 for type 2)
  f = finufft2d1(x,y,rand(M,1),+1,tol,n,n,o);   % warms up FFTW
  tic;
  for r=1:reps
    f = finufft2d1(x,y,c,isign,tol,n,n,o);
  end
  t=toc; fprintf('2d1 finufft\t\t%.3g s   \t(%.3g NUpt/s)\n',t,reps*M/t)
  
  addpath ~/numerics/nufft/nfft-3.3.2/matlab/nfft
  m=ceil(log10(1/tol)/2);       % spread width is 2m+1
  nf=2^(ceil(log2(max(n)))+1);  % fine grid: lowest 2^k that's at least 2n
  % (maybe this should be set to finufft's choice which is better for eg n=130!)
  flags=NFFT_OMP_BLOCKWISE_ADJOINT;
  flags=bitor(FFT_OUT_OF_PLACE,flags);
  flags=bitor(bitshift(uint32(1),11),flags); % NFFT_SORT_NODES
  flags=bitor(PRE_PSI,flags);
  fftw_flag=bitor(FFTW_MEASURE,1);  % FFTW_DESTROY_INPUT
  plan=nfft(2,[n;n],M,nf,nf,m,flags,fftw_flag);
  plan.x = [x y]*(1/2/pi);          % rescale NU
  nfft_precompute_psi(plan);
  plan.f = rand(M,1); nfft_adjoint(plan);   % warms up FFTW
  tic;
  for r=1:reps
    plan.f = c;      % send strengths in
    nfft_adjoint(plan);
    f2 = plan.fhat;  % get coeffs out
  end
  t2=toc;
  fprintf('2d1 NFFT pre psi\t%.3g s   \t(%.3g NUpt/s)\tspeedup:%.3g\n',t2,reps*M/t2, t2/t)
  speedups(i)=t2/t;
  delete(plan);
  f = f(:);
  fprintf('rel l2 diff btw finufft & NFFT = %.3g \t (req=%.3g)\n',norm(f-f2)/norm(f),tol)
end                                                   % ..........

figure; semilogx(densities,speedups,'+-'); xlabel('M/N');
ylabel('t_{NFFT}/t_{finufft}'); hold on; plot(densities,1+0*densities,'r:');
title(sprintf('2d1 finufft speedup vs NFFT: n=%d, nthr=%d, tol=%.3g',n,nthr,tol))
v=axis; v(3)=0; axis(v);   % y axis start at zero

%print -dpng 2d1_speedups_1e-6_n256_i7_8thr.png