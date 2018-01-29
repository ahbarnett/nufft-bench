function outfile = benchallcodes(ty,dim,N,M,nudist,multithreaded,nfftpres,name)
% BENCHALLCODES   benchmark various NUFFT codes, for fixed dim, type, NU distn.
%
% benchallcodes(type,dim,N,M,nudist,multithreaded,nfftpres,name,verb), with each
%  input argument optional, prints summary data, writes a results .mat file.
%
% Inputs:
%  type = 1 or 2 (in latter case, isign=-1 to match nfft).
%  dim = 1,2, or 3. Spatial dimension. *** only dim=3 implemented for now.
%  N = # modes in each direction.
%  M = # NU pts.
%  nudist = 0 (quasi-unif in cube), 1 (radial, ie 1/r^2 density)
%  multithreaded = 0 (single core) or 1 (all cores; user could control ahead
%                  of time)
%  nfftpres = which NFFT PSI precompute stage to test up to (PRE_PHI always on).
%             0: no pre psi,  1: also test pre psi,  2: also test pre full psi.
%                               (needs ~ 600B/NUpt)    (needs ~ 35kB/NUpt)
%  name = text string giving cpu/expt name (used to name the output file)
% Output:
%  outfile = filename written to (without .mat suffix)
%
% Codes currently benchmarked:                   run indices for each code:
%  FINUFFT                                       jf
%  NFFT w/ above precompute opts                 jn, jnp, jnf
%  CMCL (only if single thread)                  jc
%
% Magland; Barnett edit of 1/23/18, unify 1/28/18

if nargin<1, ty=1; end                             % defaults
if nargin<2 || isempty(dim), dim=3; end
if nargin<3 || isempty(N), N=32; end
if nargin<4 || isempty(M), M=1e5; end
if nargin<5 || isempty(nudist), nudist=0; end
if nargin<6 || isempty(multithreaded), multithreaded=1; end
if nargin<7 || isempty(nfftpres), nfftpres = 1; end
if nargin<8, name='generic'; end

setuppaths(multithreaded);
nthreads = 1;
if multithreaded, nthreads = java.lang.Runtime.getRuntime().availableProcessors; end
%https://stackoverflow.com/questions/8311426/how-can-i-query-the-number-of-physical-cores-from-matlab

N1=N; Ns = N1; if dim>1, N2=N; Ns = [N1 N2]; end
if dim>2, N3=N; Ns = [N1 N2 N3]; end  %   output mode box, and list of N's

% FFTW opts for everybody...
fftw_meas=1;   % 0 for ESTIMATE, 1 for MEASURE, in all algs that use FFTW

isign = 1; if ty==2, isign=-1; end   % overall isign in others to match NFFT's defn

rng(1);
[x y z nudnam] = nudata(dim,nudist,M);    % NU locs
if ty==1
  data_in = randn(M,1) + 1i*randn(M,1);   % nu pt strengths
else
  data_in = randn([Ns 1]) + 1i*randn([Ns 1]);  % Fourier coeffs, handles all dim
end

% pick outfile name...
outfile=sprintf('results/%s_%dthread_%dd%d_N%d_M%d%s_np%d',name,nthreads,dim,ty,N,M,nudnam,nfftpres);

ALGS={};  % build the ALGS struct encoding all runs ========================
% *** need to fix rest of code for dim=2,1...

% "truth" (use finufft at high acc)
eps=1e-14;
ALG=struct;
ALG.algtype=0;
ALG.name=sprintf('truth');
ALG.algopts.eps=eps;
ALG.algopts.opts.nthreads=0;       % all threads 
ALG.algopts.opts.spread_sort=1;
ALG.algopts.isign=isign;
ALG.algopts.opts.fftw=fftw_meas;   % serves to precompute FFTW for finufft tests
ALG.init=@dummy_init;
if ty==1, ALG.run=@run_finufft3d1; else, ALG.run=@run_finufft3d2; end
ALGS{end+1}=ALG;

% finufft
epsilons=10.^(-(2:12));      % range of accuracies
finufft_algopts.opts.nthreads=nthreads;
finufft_algopts.opts.spread_sort=1;               % whether to sort
finufft_algopts.opts.fftw=fftw_meas;
finufft_algopts.isign=isign;
finufft_algopts.opts.debug=0;
for ieps=1:length(epsilons)
  eps=epsilons(ieps);
  ALG=struct;
  ALG.algtype=1;
  ALG.name=sprintf('finufft(%g)',eps);
  ALG.algopts=finufft_algopts;
  ALG.algopts.eps=eps;
  ALG.init=@dummy_init;
  if ty==1, ALG.run=@run_finufft3d1; else, ALG.run=@run_finufft3d2; end
  ALGS{end+1}=ALG;
end;

% nfft
nfft_algopts.fftw_measure=fftw_meas;  % opts for all accs
nfft_algopts.reshape=0;       % only relevant for type-1
nfft_algopts.precompute_phi=1;  % 0 or 1.
nfft_algopts.precompute_psi=1;
if 1 % create nfft dummy run to initialize with fftw_measure
    ALG=struct;
    ALG.algtype=0;
    ALG.name=sprintf('nfft(1) - dummy');
    ALG.algopts=nfft_algopts;
    ALG.algopts.m=1;
    ALG.init=@init_nfft;
    if ty==1, ALG.run=@run_nfft_adj; else, ALG.run=@run_nfft; end
    ALGS{end+1}=ALG;
end;
m_s=1:6;   % NFFT range of accuracies (kernel integer half-widths)
for pretype = 0:nfftpres
  nfft_algopts.precompute_psi=pretype;
  for im=1:length(m_s)
    ALG=struct;
    ALG.algtype=2 + 0.25*pretype;          % silly way to encode which setting
    ALG.name=sprintf('nfft(%d,pre=%d)',m_s(im),pretype);
    ALG.algopts=nfft_algopts;
    ALG.algopts.m=m_s(im);
    ALG.init=@init_nfft;
    if ty==1, ALG.run=@run_nfft_adj; else, ALG.run=@run_nfft; end
    ALGS{end+1}=ALG;
  end;
end

% cmcl
if ~multithreaded
  epsilons=10.^(-(1:11));      % range of target accuracies
  cmcl_algopts.isign=isign;
  for ieps=1:length(epsilons)
    eps=epsilons(ieps);
    ALG=struct;
    ALG.algtype=3;
    ALG.name=sprintf('CMCL(%g)',eps);
    ALG.algopts=cmcl_algopts;
    ALG.algopts.eps=eps;
    ALG.init=@dummy_init;
    if ty==1, ALG.run=@run_nufft3d1; else, ALG.run=@run_nufft3d2; end
    ALGS{end+1}=ALG;
  end
end
% ========================================================================

% do all expts!
results=run_algs(ALGS,x,y,z,data_in,N1,N2,N3);

% text reporting...
fprintf('\n\n%15s %15s %15s %15s %15s\n','name','init_time(s)','run_time(s)','rel2err','therr');
X0=results{1}.data_out;         % truth
outputl1 = sum(abs(X0(:)));
outputl2 = norm(X0(:));
inputl1 = sum(abs(data_in(:)));  % appears in theoretical err estimate.
for j=1:length(ALGS)            % generate error values, text output
  X1=results{j}.data_out;
  rel2err(j)=norm(X0(:)-X1(:))/outputl2;
  rel1err(j)=norm(X0(:)-X1(:),1)/outputl1;
  maxerr(j)=max(abs(X0(:)-X1(:)));
  theoryerr(j)=maxerr(j)/inputl1;      % theoretical eps
  init_times(j)=results{j}.init_time;
  run_times(j)=results{j}.run_time;
  total_times(j)=results{j}.tot_time;
  algtypes(j)=ALGS{j}.algtype;
  fprintf('%15s %15g %15g %15g %15g\n',ALGS{j}.name,init_times(j),run_times(j),rel2err(j),theoryerr(j));
end
jf =find(algtypes==1);  % finufft
jn =find(algtypes==2);  % nfft
jnp=find(algtypes==2.25);  % nfft no pre (see silly encoding above)
jnf=find(algtypes==2.5);  % nfft full pre
jc =find(algtypes==3);  % cmcl

clear i j x y z data_in X0 X1 results % kill the big stuff
save(outfile);                    % save just error stats to disk (small)
%%%%%%%%%%%%%%%%%%%% done


function results=run_algs(ALGS,x,y,z,data_in,N1,N2,N3)
% jfm. run list of expts specified in ALGS struct
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


%%%%%%%%%%%% WRAPPERS FOR THE VARIOUS 3D NUFFT CODES %%%%%%%%%%%%%%%%%%%

function init_data=dummy_init(algopts,x,y,z,N1,N2,N3)
init_data=struct;


% finufft
function data_out=run_finufft3d1(algopts,x,y,z,data_in,N1,N2,N3,init_data)
[data_out ier]=finufft3d1(x,y,z,data_in,algopts.isign,algopts.eps,N1,N2,N3,algopts.opts);

function data_out=run_finufft3d2(algopts,x,y,z,data_in,N1,N2,N3,init_data)
[data_out ier]=finufft3d2(x,y,z,algopts.isign,algopts.eps,data_in,algopts.opts);


% cmcl
function data_out=run_nufft3d1(algopts,x,y,z,data_in,N1,N2,N3,init_data)
M=length(x);
[data_out ier]=nufft3d1(M,x,y,z,data_in,algopts.isign,algopts.eps,N1,N2,N3);
data_out = M*data_out;    % cmcl normalization

function data_out=run_nufft3d2(algopts,x,y,z,data_in,N1,N2,N3,init_data)
M=length(x);
[data_out ier]=nufft3d2(M,x,y,z,algopts.isign,algopts.eps,N1,N2,N3,data_in);


% nfft
function init_data=init_nfft(algopts,x,y,z,N1,N2,N3)
M=length(x);
N=[N1;N2;N3];       % outgoing #s of modes. note N is a vector!
ticA=tic;   % set up the plan...
%plan=nfft(3,N,M);  % default interface gives 13 digits, uses full kernel width
%... or use guru options, as commented in nfft-3.3.2/matlab/nfft/test_nfft3d.m :
n=2^(ceil(log(max(N))/log(2))+1); % lowest power of 2 that's at least 2N
flags=NFFT_OMP_BLOCKWISE_ADJOINT;         % choose more NFFT flags
flags=bitor(FFT_OUT_OF_PLACE,flags);
flags=bitor(bitshift(uint32(1),11),flags); % NFFT_SORT_NODES
if algopts.precompute_psi==1
  flags=bitor(PRE_PSI,flags);        % uses 60 GB for M=1e8, m=6
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
plan.x=[x y z]/(2*pi);   % set M*3 array of nodes in plan
fprintf('  Setting nodes in plan: %g s\n',toc(ticA));

ticA=tic;
nfft_precompute_psi(plan);                     % precomputations
fprintf('  time for nfft_precompute_psi: %g s\n',toc(ticA));

init_data.plan=plan;

function X=run_nfft_adj(algopts,x,y,z,d,N1,N2,N3,init_data)
% type 1 ("adjoint")
plan=init_data.plan;   % already carries NU locs (x y z args unused)
M=length(x);
plan.f=d;
nfft_adjoint(plan);  % type 1
if algopts.reshape
    X=reshape(plan.fhat,[N1,N2,N3]);            % ahb removed fac 1/M, 10/24/17
else
    X=plan.fhat;
end;

function X=run_nfft(algopts,x,y,z,d,N1,N2,N3,init_data)
% type 2 ("non-adjoint")
plan=init_data.plan;   % already carries NU locs (x y z args unused)
M=length(x);
plan.fhat=d(:);     % unwrap to single col
nfft_trafo(plan);   % type 2
X=plan.f;
