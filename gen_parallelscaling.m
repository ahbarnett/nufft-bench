% parallel scaling of finufft w/ # threads.
% Barnett 3/19/18. added reps 5/15/18, split off fig plot.

addpath ../finufft/matlab/
clear
dim=3; nudist=4;
N1 = 100; N2=N1;N3=N1; % # modes in each dim. fixed.
dir = 'results/xeon/';  % output

for cc=2 % 1:2    % ................ loop over cases ......................
if cc==1
  baseM = 3e6;   % high acc
  eps = 1e-12;
elseif cc==2
  baseM = 2e7;   % low acc
  eps = 1e-3;
end

o = [];  % finufft opts
o.debug = 1; o.chkbnds = 0; o.fftw = 0;      % default
reps=3;   % best of
nthrs = [1 2 3 4 8 12 16 24];  % thread numbers to try

nthrold = maxNumCompThreads;
nthr_avail = java.lang.Runtime.getRuntime().availableProcessors;

% STRONG SCALING ---------------------------------------------
M=baseM;
[x y z nudnam] = nudata(dim,nudist,M); M=numel(x);   % NU locs, any dim
c = randn(M,1)+1i*randn(M,1);
ts = nan(numel(nthrs),2);
for i=1:numel(nthrs)   % run the expt
  nthr = nthrs(i);
  maxNumCompThreads(nthr);
  for r=1:reps
    t = tic; [f ier]=finufft3d1(x,y,z,c,+1,eps,N1,N2,N3,o);
    tlist(r) = toc(t);
  end
  ts(i,1) = min(tlist);
  for r=1:reps
    t = tic; [c ier]=finufft3d2(x,y,z,+1,eps,f,o);
    tlist(r) = toc(t);
  end
  ts(i,2) = min(tlist);
  fprintf('%d threads:\t 3d1 %.3g s \t 3d2 %.3g s\n',nthr,ts(i,1),ts(i,2))
end
maxNumCompThreads(nthrold);
nam = [dir sprintf('strongscaling_3d_case%d',cc)];  % file head
clear x y z c f;
save(nam);
disp('        nthr        par efficiencies')
disp([nthrs', ts(1,:) ./ (ts.*nthrs')])

% WEAK SCALING ---------------------------------------------
ts = nan(numel(nthrs),2);
for i=1:numel(nthrs)   % run the expt
  nthr = nthrs(i);
  M= baseM * nthr;
  maxNumCompThreads(nthr_avail);
  [x y z nudnam] = nudata(dim,nudist,M); M=numel(x);   % NU locs, any dim
  c = randn(M,1)+1i*randn(M,1);  
  maxNumCompThreads(nthr);
  for r=1:reps
    t = tic; [f ier]=finufft3d1(x,y,z,c,+1,eps,N1,N2,N3,o);
    tlist(r) = toc(t);
  end
  ts(i,1) = min(tlist);
  for r=1:reps
    t = tic; [c ier]=finufft3d2(x,y,z,+1,eps,f,o);
    tlist(r) = toc(t);
  end
  ts(i,2) = min(tlist);
  fprintf('%d threads:\t 3d1 %.3g s \t 3d2 %.3g s\n',nthr,ts(i,1),ts(i,2))
end
maxNumCompThreads(nthrold);
nam = [dir sprintf('weakscaling_3d_case%d',cc)];  % file head
clear x y z c f;
save(nam);
disp('        nthr        par efficiencies')
disp([nthrs', ts(1,:)./ts])

end  % ............. end cases ...........
