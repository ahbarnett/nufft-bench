function setuppaths(multithreaded)
% set up all paths for NUFFT comparisons. Also see README.md
% Barnett 2/20/18

if nargin<1, multithreaded=1; end    % default

% NFFT
addnfftpath(multithreaded)

% FINUFFT
addpath ../finufft/matlab/

% CMCL
addpath ../nufftall-1.33/

% BART
run ../bart-0.4.02/startup

% MIRT (Fessler)
irtdir=[fileparts(mfilename('fullpath')) '/../irt'];
run([irtdir '/setup'])
clear irtdir
