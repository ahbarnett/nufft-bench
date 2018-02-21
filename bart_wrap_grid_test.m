% get convention correct for BART NUFFT and demo wrapping failure
% Barnett 2/19/18

run ~/numerics/nufft/bart-0.4.02/startup
%bart('nufft -h')

clear
n=40;    % number of uniform output pixels (aka modes) in each dimension
k=(-n/2:n/2-1)';                            % col vec of output indices
N=n^3;   % # 3D modes (apparently needed for normalization, below)

for frac = [0.5 0.8 0.9 0.95 0.98 0.99 1]   % k-space loc as fraction of n/2
  xyz = [frac*(n/2);0;0];                   % one k-space on x axis (y=z=0)
  c = randn+1i*randn;                       % strength
  cmd=sprintf('nufft -a -d %d:%d:%d',n,n,n);
  F = bart(cmd,xyz,c);                      % do the NUFFT
  diff = sqrt(N)*F(:,1,1) - c*exp(pi*2i*k*xyz(1)/n);
  fprintf('k-space pt loc %g*n/2: \t max err = %.3g\n',frac,max(abs(diff)))
end

% 2/20/18: Got response from Martin Uecker which included new grid.c
% I inserted into src/noncart/grid.c
% and it fixed the wrap-around issue - good!

