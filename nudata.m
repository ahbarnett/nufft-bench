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
    %th = rand(M,1)*pi; costh=cos(th); sinth=sin(th); % jfm old, not unif on S2
    x = r.*sinth.*cos(phi);
    y = r.*sinth.*sin(phi);
    z = r.*costh;
  end
  
else, error('nudist > 1 not yet implemented');
         % *** add deterministic quadr scheme? or adaptive
end

%%%%%%%%%%%
function test_nudata
for nudist=0:1
  [x y z nam] = nudata(3,nudist,5e3);
  figure; plot3(x,y,z,'.'); title(nam);
  axis vis3d equal tight off; view(20,10);
  set(gcf,'paperposition',[0 0 4 4]);
  print(sprintf('results/nudist_%d.eps',nudist),'-depsc2')
end


