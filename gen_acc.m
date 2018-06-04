function [er2 eth] = gen_acc(N,M,tols,verb)
% GEN_ACC  Sweep requested tolerance and collect all dims all types errors
%
% [er2 eth] = gen_acc(N,M,tols,verb)
% verb = 0 (no figs), 1 (main figs for rel l2), 2 (also figs for theory err)
%
% Barnett 6/2/18

if nargin<4, verb = 0; end

reftol = 1e-15;
isign   = +1;     % sign of imaginary unit in exponential
o.debug = 0;      % choose 1 for timing breakdown text output
o.nthreads = 0;   % omit, or use 0, to use default num threads.
o.fftw = 0;       % style of FFTW: 0 (ESTIMATE) vs 1 (MEASURE, slow but reuses)

nt = numel(tols);
eth = nan(9,nt);  % all "theory errors" (l_infty rel to l_1 of input)
er2 = eth;        % all l_2 errors (rel to output size)

% --------- 1D --------------------------------------------------------
fprintf('1D: using %d modes...\n',N)
x = pi*(2*rand(1,M)-1);
c = randn(1,M)+1i*randn(1,M);
fe = finufft1d1(x,c,isign,reftol,N,o); % "exact"
ty = 1;
in1 = norm(c(:),1); ou2 = norm(fe(:),2);
for t=1:nt
  fprintf('1d1 tol=%.3g...\n',tols(t));
  f = finufft1d1(x,c,isign,tols(t),N,o);
  d = f(:)-fe(:);                           % output diff
  eth(ty,t) = norm(d,Inf)/in1; er2(ty,t) = norm(d,2)/ou2;
end
if verb>1     % show that theory errors are much smaller than rel l2...
  figure; loglog(tols,[eth(1,:);er2(1,:)],'+-');hold on;loglog(tols,tols,'r--')
  title(sprintf('1d1: M=%.3g N=%.3g',M,N)); legend('l1-linf','rel l2','req')
end

x = pi*(2*rand(1,M)-1);  % new NU
f = randn(1,N)+1i*randn(1,N);
ce = finufft1d2(x,isign,reftol,f,o);
ty = 2;
in1 = norm(f(:),1); ou2 = norm(ce(:),2);
for t=1:nt
  fprintf('1d2 tol=%.3g...\n',tols(t));
  c = finufft1d2(x,isign,tols(t),f,o);
  d = c(:)-ce(:);                           % output diff
  eth(ty,t) = norm(d,Inf)/in1; er2(ty,t) = norm(d,2)/ou2;
end

x = pi*(2*rand(1,M)-1);  % new NU
c = randn(1,M)+1i*randn(1,M);
s = (N/2)*(2*rand(1,M)-1);                      % target freqs of size O(N)
ty = 3;
[fe ier] = finufft1d3(x,c,isign,reftol,s,o);
in1 = norm(c(:),1); ou2 = norm(fe(:),2);
for t=1:nt
  fprintf('1d3 tol=%.3g...\n',tols(t));
  f = finufft1d3(x,c,isign,tols(t),s,o);
  d = f(:)-fe(:);                           % output diff
  eth(ty,t) = norm(d,Inf)/in1; er2(ty,t) = norm(d,2)/ou2;
end
if verb
figure; loglog(tols,er2(1:3,:),'+-'); hold on;loglog(tols,tols,'k-')
title(sprintf('1D: M=%.3g N=%.3g',M,N)); axis tight
legend('type 1','type 2','type 3','req tol','location','northwest')
end

% --------- 2D ---------------------------------------------------------------
N1=ceil(sqrt(N)); N2=round(N/N1);           % pick Fourier mode ranges
fprintf('2D: using %d*%d modes (total %d)...\n',N1,N2,N1*N2)
x = pi*(2*rand(1,M)-1); y = pi*(2*rand(1,M)-1);
c = randn(1,M)+1i*randn(1,M);
ty = 4;
fe = finufft2d1(x,y,c,isign,reftol,N1,N2,o); % "exact"
in1 = norm(c(:),1); ou2 = norm(fe(:),2);
for t=1:nt
  fprintf('2d1 tol=%.3g...\n',tols(t));
  f = finufft2d1(x,y,c,isign,tols(t),N1,N2,o);
  d = f(:)-fe(:);                           % output diff
  eth(ty,t) = norm(d,Inf)/in1; er2(ty,t) = norm(d,2)/ou2;
end
if verb>1
  figure; loglog(tols,[eth(4,:);er2(4,:)],'+-');hold on;loglog(tols,tols,'r--')
  title(sprintf('2d1: M=%.3g N=%.3g',M,N)); legend('l1-linf','rel l2','req')
end

x = pi*(2*rand(1,M)-1); y = pi*(2*rand(1,M)-1); % new NU
f = randn(N1,N2)+1i*randn(N1,N2);
ce = finufft2d2(x,y,isign,reftol,f,o);
ty = 5;
in1 = norm(f(:),1); ou2 = norm(ce(:),2);
for t=1:nt
  fprintf('2d2 tol=%.3g...\n',tols(t));
  c = finufft2d2(x,y,isign,tols(t),f,o);
  d = c(:)-ce(:);                           % output diff
  eth(ty,t) = norm(d,Inf)/in1; er2(ty,t) = norm(d,2)/ou2;
end

x = pi*(2*rand(1,M)-1); y = pi*(2*rand(1,M)-1); % new NU
c = randn(1,M)+1i*randn(1,M);
s = (N1/2)*(2*rand(1,M)-1);                      % target freqs of size O(N1)
u = (N2/2)*(2*rand(1,M)-1);                      % target freqs of size O(N2)
ty = 6;
[fe ier] = finufft2d3(x,y,c,isign,reftol,s,u,o);
in1 = norm(c(:),1); ou2 = norm(fe(:),2);
for t=1:nt
  fprintf('2d3 tol=%.3g...\n',tols(t));
  f = finufft2d3(x,y,c,isign,tols(t),s,u,o);
  d = f(:)-fe(:);                           % output diff
  eth(ty,t) = norm(d,Inf)/in1; er2(ty,t) = norm(d,2)/ou2;
end
if verb
figure; loglog(tols,er2(4:6,:),'+-'); hold on;loglog(tols,tols,'k-')
title(sprintf('2D: M=%.3g N=%.3g',M,N)); axis tight
legend('type 1','type 2','type 3','req tol','location','northwest')
end

% --------- 3D ---------------------------------------------------------------
N1=ceil(N^(1/3)); N2=N1; N3=round(N/N1/N2);  % pick Fourier mode ranges
fprintf('3D: using %d*%d*%d modes (total %d)...\n',N1,N2,N3,N1*N2*N3)
x = pi*(2*rand(1,M)-1); y = pi*(2*rand(1,M)-1); z = pi*(2*rand(1,M)-1);
c = randn(1,M)+1i*randn(1,M);
ty = 7;
fe = finufft3d1(x,y,z,c,isign,reftol,N1,N2,N3,o); % "exact"
in1 = norm(c(:),1); ou2 = norm(fe(:),2);
for t=1:nt
  fprintf('3d1 tol=%.3g...\n',tols(t));
  f = finufft3d1(x,y,z,c,isign,tols(t),N1,N2,N3,o);
  d = f(:)-fe(:);                           % output diff
  eth(ty,t) = norm(d,Inf)/in1; er2(ty,t) = norm(d,2)/ou2;
end
if verb>1
  figure; loglog(tols,[eth(7,:);er2(7,:)],'+-');hold on;loglog(tols,tols,'r--')
  title(sprintf('3d1: M=%.3g N=%.3g',M,N)); legend('l1-linf','rel l2','req')
end

x = pi*(2*rand(1,M)-1); y = pi*(2*rand(1,M)-1); z = pi*(rand(1,M)-1); % new NU
f = randn(N1,N2,N3)+1i*randn(N1,N2,N3);
ce = finufft3d2(x,y,z,isign,reftol,f,o);
ty = 8;
in1 = norm(f(:),1); ou2 = norm(ce(:),2);
for t=1:nt
  fprintf('3d2 tol=%.3g...\n',tols(t));
  c = finufft3d2(x,y,z,isign,tols(t),f,o);
  d = c(:)-ce(:);                           % output diff
  eth(ty,t) = norm(d,Inf)/in1; er2(ty,t) = norm(d,2)/ou2;
end

x = pi*(2*rand(1,M)-1); y = pi*(2*rand(1,M)-1); z = pi*(rand(1,M)-1); % new NU
c = randn(1,M)+1i*randn(1,M);
s = (N1/2)*(2*rand(1,M)-1);                      % target freqs of size O(N1)
u = (N2/2)*(2*rand(1,M)-1);                      % target freqs of size O(N2)
v = (N3/2)*(2*rand(1,M)-1);                      % target freqs of size O(N3)
ty = 9;
[fe ier] = finufft3d3(x,y,z,c,isign,reftol,s,u,v,o);
in1 = norm(c(:),1); ou2 = norm(fe(:),2);
for t=1:nt
  fprintf('3d3 tol=%.3g...\n',tols(t));
  f = finufft3d3(x,y,z,c,isign,tols(t),s,u,v,o);
  d = f(:)-fe(:);                           % output diff
  eth(ty,t) = norm(d,Inf)/in1; er2(ty,t) = norm(d,2)/ou2;
end
if verb
figure; loglog(tols,er2(7:9,:),'+-'); hold on;loglog(tols,tols,'k-')
title(sprintf('3D: M=%.3g N=%.3g',M,N)); axis tight
legend('type 1','type 2','type 3','req tol','location','northwest')
end