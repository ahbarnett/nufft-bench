% strong scaling w/ # threads
% Barnett 3/19/18


addpath ../finufft/matlab/
clear

dim=3; nudist=1;
M=1e7;
[x y z nudnam] = nudata(dim,nudist,M);    % NU locs, any dim
