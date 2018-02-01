function [x y z name] = nudata(dim,nudist,M)
% NUDATA.  make non-uniform data for NUFFT tests, lying in [-pi,pi]^d.
%
% [x y z name] = nudata(dim,nudist,M)
%
% Inputs: dim=1,2, or 3.
%  nudist = 0 (unif in cube), 1 (spherical)
%  M = number of NU pts
%
% Outputs: the NU locations are
%  in x (for 1D), x y (for 2D), x y z (for 3D). unused outputs are empty.
%  name = string abbrev of the distribution.
%
% Without args, does 3d test and makes EPS figs

% Barnett 11/7/17, changed interface to include name, 1/28/18
if nargin==0, test_nudata; return; end

y = []; z = [];
x = pi*(2*rand(M,1)-1);

if nudist==0  % quasi-uniform
  name = 'unif';
  if dim>1, y = pi*(2*rand(M,1)-1); end
  if dim>2, z = pi*(2*rand(M,1)-1); end

elseif nudist==1  % radial (unif in 1D, 1/r density in 2D, 1/r^2 in 3D)
  name = 'sph';
  if dim>1, r = pi*rand(M,1); phi = 2*pi*rand(M,1); end
  if dim==2
    x = r.*cos(phi);
    y = r.*sin(phi);
  elseif dim==3
    costh = 2*rand(M,1)-1; sinth = sqrt(1-costh.^2);
    x = r.*sinth.*cos(phi);
    y = r.*sinth.*sin(phi);
    z = r.*costh;
  end

elseif nudist==2    % reprod code from expts of 11/7/17 used for FSU talk:
  name = 'sphwrong';
  if dim==3
    theta=rand(M,1)*2*pi;
    phi=acos(rand(M,1));   % ahb attempt for unif on S2, but was wrong!
    rad=rand(M,1)*pi;
    x=rad.*cos(theta).*cos(phi);
    y=rad.*sin(theta).*cos(phi);
    z=rad.*sin(phi);
    %temp = x; x=z; z= temp;  % swap, made nfft faster, strangely
  end
  
elseif nudist==3  % unif in small cuboid
  name = 'x-flat pancake';      % only uses 1-2 cores for spreading!
  x = x/10;
  %x(M/10:end) = x(M/10:end)/10;  % some unif, some pancake...
  if dim>1, y = pi*(2*rand(M,1)-1); end
  if dim>2, z = pi*(2*rand(M,1)-1); end

elseif nudist==4  % deterministic quadr scheme pts, no weights
  name = 'sphquad';
  if dim==3
    m = round(M^(1/3)); if mod(m,2)==1, m=m+1; end  % m even, # grid pts/dim
    rad = pi; %/2;  % radius, arb, <=pi otherwise exceeds cube walls
    mug=gauss(m);  % mu from -1 to 1, z for unit sph.  col vec
    mp = 2*m;      % # phi pts (equator)
    pg = 2*pi*(1:mp)'/mp;
    mr = m/2;      % # radial pts
    rg = rad*(1+gauss(mr))/2;  % from [0,rad], would have r^2 metric
    sthg = sqrt(1-mug.^2);
    sth = kron(kron(ones(mp,1),sthg),ones(mr,1));  % all sin th's
    cth = kron(kron(ones(mp,1),mug),ones(mr,1));   % all cos th's
    cpg = cos(pg); spg = sin(pg);
    r = kron(ones(mp*m,1),rg);
    x = r.*sth.*kron(cpg,ones(m*mr,1));
    y = r.*sth.*kron(spg,ones(m*mr,1));
    z = r.*cth;
  end

else, error('nudist > 4 not yet implemented');
end

%%%%%%%%%%%
function test_nudata
for nudist=0:4
  [x y z nam] = nudata(3,nudist,5e3);
  figure; plot3(x,y,z,'.'); title(nam);
  axis vis3d equal tight off; view(20,10);
  set(gcf,'paperposition',[0 0 4 4]);
  print(sprintf('results/nudist_%d.eps',nudist),'-depsc2')
end


