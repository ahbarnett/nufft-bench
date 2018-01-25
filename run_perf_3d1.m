function outfile = run_perf_3d1(NN,M,nudist,multithreaded,nfftpre,cpuname,verb)
% FIG_PERF_3D1.  compare NUFFT 3D type-1 performance of various codes
%
% fig_perf_3d1(NN,M,nudist,multithreaded,nfftpre,cpuname,verb), with each input
%  argument optional, prints summary data, maybe plots a figure, and writes
%  a results file.
%
% Inputs:
%  NN = # modes in each direction. M = # NU pts.
%  nudist = 0 (quasi-unif), 1 (radial, ie 1/r^2 dist in 3D),
%  multithreaded = 0 (single core) or 1 (all cores; user should control ahead
%                  of time).
%  nfftpre = 0 (no pre psi), 1 (pre psi), or 2 (full pre psi)
%  cpuname = text string giving cpu name (used to name the output file) 
%  verb = 0 (no figure), 1 (make figure)
% Output:
%  outfile = file name written to (without .mat suffix)
%
% Magland; Barnett edit of 1/23/18

if nargin<1 || isempty(NN), NN=64; end    % defaults
if nargin<2 || isempty(M), M=1e5; end
if nargin<3 || isempty(nudist), nudist=0; end
if nargin<4 || isempty(multithreaded), multithreaded=1; end
if nargin<5 || isempty(nfftpre), nfftpre = 1; end
if nargin<6, cpuname='unknowncpu'; end
if nargin<7, verb=0; end

setuppaths(multithreaded);
nthreads = java.lang.Runtime.getRuntime().availableProcessors;
%https://stackoverflow.com/questions/8311426/how-can-i-query-the-number-of-physical-cores-from-matlab
if ~multithreaded, nthreads=0; end
N1=NN; N2=NN; N3=NN;  % the output mode cube

% FFTW opts for everybody...
fftw_meas=1;   % 0 for ESTIMATE, 1 for MEASURE, in all algs that use FFTW

rng(1); [x y z] = nudata(nudist,M);  % gen data
data_in = randn(M,1) + 1i*randn(M,1);   % nu pt strengths

title0=sprintf('Type 1, %dx%dx%d, %g sources, %d threads',NN,NN,NN,M,nthreads);
if multithreaded, thrstr='multi'; else thrstr='single'; end
if nudist, nudstr='sph'; else nudstr='unif'; end
nfftstr=''; if nfftpre==0, nfftstr='_nopre'; elseif nfftpre==2, nfftstr='_full';end

outfile=sprintf('results/%s_%s_3d1_N%d_M%d%s%s',cpuname,thrstr,NN,M,nudstr,nfftstr);

ALGS={};  % build the ALGS struct encoding all runs ========================

%"truth"
eps=1e-14;
ALG=struct;
ALG.algtype=0;
ALG.name=sprintf('truth',eps);
ALG.algopts.eps=eps;
ALG.algopts.opts.nthreads=0;  % all threads 
ALG.algopts.opts.spread_sort=1;
ALG.algopts.isign=1;
ALG.algopts.opts.fftw=fftw_meas;   % serves to precompute FFTW for finufft tests
ALG.init=@dummy_init; ALG.run=@run_finufft3d1;
ALGS{end+1}=ALG;

% finufft
epsilons=10.^(-(2:12));      % range of accuracies
finufft_algopts.opts.nthreads=nthreads;
finufft_algopts.opts.spread_sort=1;               % whether to sort
finufft_algopts.opts.fftw=fftw_meas;
finufft_algopts.isign=1;
finufft_algopts.opts.debug=0;
for ieps=1:length(epsilons)
    eps=epsilons(ieps);
    ALG=struct;
    ALG.algtype=1;
    ALG.name=sprintf('finufft(%g)',eps);
    ALG.algopts=finufft_algopts;
    ALG.algopts.eps=eps;
    ALG.init=@dummy_init; ALG.run=@run_finufft3d1;
    ALGS{end+1}=ALG;
end;

% nfft
nfft_algopts.fftw_measure=fftw_meas;  % opts for all accs
nfft_algopts.reshape=0;
nfft_algopts.precompute_phi=1;  % 0 or 1.
nfft_algopts.precompute_psi=nfftpre;
if 1 % create nfft dummy run to initialize with fftw_measure
    ALG=struct;
    ALG.algtype=0;
    ALG.name=sprintf('nfft(1) - dummy');
    ALG.algopts=nfft_algopts;
    ALG.algopts.m=1;
    ALG.init=@init_nfft; ALG.run=@run_nfft;
    ALGS{end+1}=ALG;
end;
m_s=1:6;   % NFFT range of accuracies (kernel integer half-widths)
for im=1:length(m_s)
    ALG=struct;
    ALG.algtype=2;
    ALG.name=sprintf('nfft(%d)',m_s(im));
    ALG.algopts=nfft_algopts;
    ALG.algopts.m=m_s(im);
    ALG.init=@init_nfft; ALG.run=@run_nfft;
    ALGS{end+1}=ALG;
end;

% cmcl
if multithreaded==0
  cmcl_algopts.isign=1;
  for ieps=1:length(epsilons)
    eps=epsilons(ieps);
    ALG=struct;
    ALG.algtype=3;
    ALG.name=sprintf('CMCL(%g)',eps);
    ALG.algopts=cmcl_algopts;
    ALG.algopts.eps=eps;
    ALG.init=@dummy_init; ALG.run=@run_nufft3d1;
    ALGS{end+1}=ALG;
  end
end
% ========================================================================

% do all expts
results=run_algs(ALGS,x,y,z,data_in,N1,N2,N3);
%print_accuracy_comparison_and_timings(ALGS,results);

errors=[];
init_times=[];
run_times=[];
total_times=[];
algtypes=[];
X0=results{1}.data_out;
fprintf('\n\n%15s %15s %15s %15s %15s\n','name','init_time(s)','run_time(s)','tot_time(s)','err');
for j=1:length(ALGS)            % generate error values, text output
    X1=results{j}.data_out;
    %m0=(max(abs(X0(:)))+max(abs(X0(:))))/2;
    %max_diff0=max(abs(X0(:)-X1(:)))/m0;
    %avg_diff0=mean(abs(X0(:)-X1(:)))/m0;
    numer0=norm(X0(:)-X1(:));              % l2 norm
    denom0=norm(X0(:));
    errors(j)=numer0/denom0;
    init_times(j)=results{j}.init_time;
    run_times(j)=results{j}.run_time;
    total_times(j)=results{j}.tot_time;
    algtypes(j)=ALGS{j}.algtype;
    fprintf('%15s %15g %15g %15g %15g\n',ALGS{j}.name,init_times(j),run_times(j),total_times(j),errors(j));
end
ii1=find(algtypes==1);  % finufft
ii2=find(algtypes==2);  % nfft
ii3=find(algtypes==3);  % cmcl
clear x y z data_in X0 X1 results
save(outfile);  % save just errors to disk (small)

if verb
  figure;  % just for checking
  semilogx(errors(ii1),run_times(ii1),'b.-');
  hold on;
  semilogx(errors(ii2),run_times(ii2),'r.-');
  semilogx(errors(ii2),init_times(ii2),'r.--');
  if ~isempty(ii3), semilogx(errors(ii3),run_times(ii3),'g.-'); end
  set(gca,'fontsize',22);
  v=axis; axis([1e-12 1e-1 0 v(4)])
  xlabel('Error'); ylabel('Run time (s)');
  legend({'finufft','NFFT','NFFT init','CMCL'}); % ignores unused ones
  title(title0); drawnow
end

%keyboard
%axis([1e-12 1e-2 0 60])
%axis([1e-12 1e-2 0 28])
%set(gcf,'paperposition',[0 0 4 4]); print -depsc2 talk/multitime.eps
%%%%%%%%%%%%%%%%%%%% done


function results=run_algs(ALGS,x,y,z,data_in,N1,N2,N3)
% run list of expts specified in ALGS struct
results={};
for j=1:length(ALGS)
    ALG=ALGS{j};
    RESULT=struct;
    
    fprintf('(%d/%d) Initializing %s...\n',j,length(ALGS),ALG.name);
    tA=tic;
    init_data=ALG.init(ALG.algopts,x,y,z,N1,N2,N3);
    RESULT.init_time=toc(tA);
    
    fprintf('Running %s...\n',ALG.name);
    tA=tic;
    RESULT.data_out=ALG.run(ALG.algopts,x,y,z,data_in,N1,N2,N3,init_data);
    RESULT.run_time=toc(tA);
    
    RESULT.tot_time=RESULT.init_time+RESULT.run_time;
    
    fprintf('init_time=%.3f, run_time=%.3f, tot_time=%.3f\n',RESULT.init_time,RESULT.run_time,RESULT.tot_time);
    
    results{j}=RESULT;
end;

function print_accuracy_comparison_and_timings(ALGS,results)
MM=length(ALGS);

max_diffs=zeros(MM,MM);
avg_diffs=zeros(MM,MM);
for j1=1:MM
    X1=results{j1}.data_out;
    for j2=1:MM
        X2=results{j2}.data_out;
        m0=(max(abs(X1(:)))+max(abs(X2(:))))/2;
        max_diffs(j1,j2)=max(abs(X1(:)-X2(:)))/m0;
        avg_diffs(j1,j2)=mean(abs(X1(:)-X2(:)))/m0;
    end;
end;

fprintf('Max. differences:\n');
fprintf('%15s ','');
for j2=1:MM
    fprintf('%15s ',ALGS{j2}.name);
end;
fprintf('\n');
for j1=1:MM
    fprintf('%15s ',ALGS{j1}.name);
    for j2=1:MM
        fprintf('%15g ',max_diffs(j1,j2));
    end;
    fprintf('\n');
end;
fprintf('\n');
fprintf('Avg. differences:\n');
fprintf('%15s ','');
for j2=1:MM
    fprintf('%15s ',ALGS{j2}.name);
end;
fprintf('\n');
for j1=1:MM
    fprintf('%15s ',ALGS{j1}.name);
    for j2=1:MM
        fprintf('%15g ',avg_diffs(j1,j2));
    end;
    fprintf('\n');
end;
fprintf('\n');
%fprintf('%15s %15s %15s %15s\n','Algorithm','Init time (s)','Run time (s)','RAM (GB)');
fprintf('%15s %15s %15s %15s\n','Algorithm','Init time (s)','Run time (s)','Tot time (s)');
for j1=1:MM
    %fprintf('%15s %15.3f %15.3f %15.3f\n',ALGS{j1}.name,results{j1}.init_time,results{j1}.run_time,results{j1}.memory_used);
    fprintf('%15s %15.3f %15.3f %15.3f\n',ALGS{j1}.name,results{j1}.init_time,results{j1}.run_time,results{j1}.tot_time);
end;
fprintf('\n');

function init_data=dummy_init(algopts,x,y,z,N1,N2,N3)
init_data=struct;



%%%%%%%%%%%% WRAPPERS FOR THE VARIOUS NUFFT ALGORITHMS: %%%%%%%%%%%%%%%%%%%

% finufft
function data_out=run_finufft3d1(algopts,x,y,z,data_in,N1,N2,N3,init_data)
data_out=finufft3d1(x,y,z,data_in,algopts.isign,algopts.eps,N1,N2,N3,algopts.opts);

% cmcl
function data_out=run_nufft3d1(algopts,x,y,z,data_in,N1,N2,N3,init_data)
M=length(x);
data_out=nufft3d1(M,x,y,z,data_in,algopts.isign,algopts.eps,N1,N2,N3);
data_out = M*data_out;    % cmcl normalization


% nfft
function init_data=init_nfft(algopts,x,y,z,N1,N2,N3)
M=length(x);
N=[N1;N2;N3];       % note N is a vector!
ticA=tic;   % set up the plan...
%plan=nfft(3,N,M);  % default interface gives 13 digits, uses full kernel width
n=2^(ceil(log(max(N))/log(2))+1); % lowest power of 2 that's at least 2N
flags=NFFT_OMP_BLOCKWISE_ADJOINT;         % choose more NFFT flags
flags=bitor(FFT_OUT_OF_PLACE,flags);
flags=bitor(bitshift(uint32(1),11),flags); % NFFT_SORT_NODES
if algopts.precompute_psi==1
  flags=bitor(PRE_PSI,flags);
elseif algopts.precompute_psi==2
  flags=bitor(PRE_FULL_PSI,flags);   % uses 35 GB for M=1e6, m=6 !!
end
if (algopts.precompute_phi)
  flags=bitor(PRE_PHI_HUT,flags);
end;
fftw_flag=FFTW_ESTIMATE;
if (algopts.fftw_measure)
  fftw_flag=FFTW_MEASURE;
end;
fftw_flag = bitor(fftw_flag, 1);    % FFTW_DESTROY_INPUT
plan=nfft(3,N,M,n,n,n,algopts.m,flags,fftw_flag); % use of nfft_init_guru
fprintf('  Creating nfft plan: %g s\n',toc(ticA));

ticA=tic;
xyz=[x y z];
fprintf('  Concatenating locations: %g s\n',toc(ticA));

ticA=tic;
plan.x=(xyz)/(2*pi); % set nodes in plan
fprintf('  Setting nodes in plan: %g s\n',toc(ticA));

ticA=tic;
nfft_precompute_psi(plan);                     % precomputations
fprintf('  time for nfft_precompute_psi: %g s\n',toc(ticA));

init_data.plan=plan;

function X=run_nfft(algopts,x,y,z,d,N1,N2,N3,init_data)
plan=init_data.plan;
M=length(x);
plan.f=d;
nfft_adjoint(plan);  % type 1
if algopts.reshape
    X=reshape(plan.fhat,[N1,N2,N3]);            % ahb removed fac 1/M, 10/24/17
else
    X=plan.fhat;
end;
