% basic timing of CMCL

nj = 1e6;
ms = 20;
mt = 20;
mu = 20;
eps = 1e-6;

xj = (rand(nj,1)*2-1)*pi;
yj = (rand(nj,1)*2-1)*pi;
zj = (rand(nj,1)*2-1)*pi;
cj = randn(nj,1)+1i*randn(nj,1);
iflag = +1;

tic
fk1 = nufft3d1(nj,xj,yj,zj,cj,iflag,eps,ms,mt,mu);
t=toc;
fprintf('3d1 done in %.3g s:  %.3g NUpts/s\n',t,nj/t)

% rebuilding with -O3 -ffast-math -funroll-loops -march=native 
% got 10-20% improvement.
