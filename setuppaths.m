function setuppaths(multithreaded)
% set up all paths for NUFFT comparisons
% Barnett 1/22/18

% NFFT
addnfftpath(multithreaded)

% FINUFFT
addpath ~/numerics/finufft/matlab/

% CMCL
addpath ~/numerics/nufftall-1.33/
