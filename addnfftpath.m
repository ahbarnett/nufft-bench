function addnfftpath(multithreaded)
% ADDNFFTPATH.  adds NFFT to matlab path. Must match your NFFT installation.

% Barnett 11/7/17
% This should match your nfft installation; see README.md

wd = fileparts(mfilename('fullpath'));   % directory of this script
base = [wd '/..'];                       % dir above the nfft installations
mul = [base '/nfft-3.3.2/matlab/nfft'];
sng = [base '/nfft-3.3.2-single-thread/matlab/nfft'];
warning('off','MATLAB:rmpath:DirNotFound');
if multithreaded
  addpath(mul)
  rmpath(sng)
else
  addpath(sng)
  rmpath(mul)
end
