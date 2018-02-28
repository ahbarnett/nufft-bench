function benchallcodes(ty,dim,N,M,nudist,multithreaded,outname,o)
% BENCHALLCODES   benchmark various NUFFT codes, for fixed dim, type, NU distn.
%
% benchallcodes(type,dim,N,M,nudist,multithreaded,filename,opts), with each
%  input argument optional, prints summary data, writes to filename.mat
%
% Inputs:
%  type = 1 ("adjoint") or 2 ("forward"; here isign=-1 to match nfft).
%  dim = 1, 2 or 3. Spatial dimension.
%  N = # modes in each dimension (same for now)
%  M = # NU pts.
%  nudist = 0 (quasi-unif in cube), 1,...4 etc; see NUDIST.
%  multithreaded = 0 (single core) or 1 (all logical cores; user could control
%                  ahead of time)
%  outname = string giving output file (excluding .mat suffix)
%  opts has fields controlling:
%   nfftpres = which NFFT PSI precompute stage to test up to (PRE_PHI always on)
%              0: no pre psi,  1: also test pre psi,  2: also test pre full psi.
%                                (3d needs ~ 600B/NUpt)  (needs ~ 35kB/NUpt)
%
% Codes currently benchmarked:                   outname indices for each code:
%  FINUFFT                                       jf
%  NFFT w/ the above precompute opts             jn, jnp, jnf
%  CMCL (only if single thread)                  jc
%  BART (only 3D)                                jb
%  MIRT                                          jm
%
% Called without arguments, does a self-test.
%
% Magland; Barnett edit of 1/23/18, unify 1/28/18, all dims 2/20/18
%
% Issues: 1) BART shows large init time if nfftpres=2, which is fake (actually 0)
%          => ignore BART's init time.
%         2) all data_out copies from all algs stored in results. can get big.
%         if M=1e8 (it's 1.6GB for type 2, per alg).

if nargin==0, benchallcodes(1); return; load temp; end        % crude self-test
if nargin<1, ty=1; end                             % defaults (for self-test)
if nargin<2 || isempty(dim), dim=3; end
if nargin<3 || isempty(N), N=32; end
if nargin<4 || isempty(M), M=1e5; end
if nargin<5 || isempty(nudist), nudist=0; end
if nargin<6 || isempty(multithreaded), multithreaded=1; end
if nargin<7, outname='temp'; end
if nargin<8, o = []; end
if ~isfield(o,'nfftpres'), o.nfftpres = 0; end  % 2

setuppaths(multithreaded);
nthreads = 1;
if multithreaded
  nthreads = java.lang.Runtime.getRuntime().availableProcessors;
end
%https://stackoverflow.com/questions/8311426/how-can-i-query-the-number-of-physical-cores-from-matlab

N1=N; N2=1; N3=1; if dim>1, N2=N; end, if dim>2, N3=N; end  % set up N1,N2,N3

% FFTW opts for everybody...
fftw_meas=1;   % 0 for ESTIMATE, 1 for MEASURE, in all algs that use FFTW

isign = 1; if ty==2, isign=-1; end   % overall isign in others to match NFFT

rng(1);
[x y z nudnam] = nudata(dim,nudist,M);    % NU locs, any dim
M = numel(x);    % since M may differ from requested
if ty==1
  data_in = randn(M,1) + 1i*randn(M,1);   % nu pt strengths
else
  Ns = [N1 N2 N3]; Ns = Ns(1:dim);
  data_in = randn([Ns 1]) + 1i*randn([Ns 1]);  % Fourier coeffs, handles all dim
end

ALGS={};  % build the ALGS struct encoding all runs ========================

% "truth" (use finufft at high acc, all threads)
eps=1e-14;
ALG=struct;
ALG.algtype=0;
ALG.name=sprintf('truth(%dd%d)',dim,ty);
ALG.algopts.eps=eps;
ALG.algopts.opts.nthreads=0;       % all threads 
ALG.algopts.opts.spread_sort=1;
ALG.algopts.isign=isign;
ALG.algopts.opts.fftw=fftw_meas;   % serves to precompute FFTW for finufft tests
ALG.init=@dummy_init;
ALG.run = str2func(sprintf('run_finufft%dd%d',dim,ty));
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
  ALG.run = ALGS{end}.run;   % same wrapper as above
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
for pretype = 0:o.nfftpres
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
    ALG.run = str2func(sprintf('run_nufft%dd%d',dim,ty));
    ALGS{end+1}=ALG;
  end
end

% bart
% Let's check that writing to SSD temp file that it does is insignificant:
% xyz = rand(3,1e7); tic; writecfl(tempname,xyz); toc    % 1 sec
% BART throughput is 5e5 to 1e6 NU pts/sec, so we're at 10% level.
if dim==3
  ALG=struct;
  ALG.algtype=4;
  ALG.name=sprintf('BART');
  ALG.algopts=[];
  ALG.init=@dummy_init;
  if ty==1, ALG.run=@run_bart_adj; else, ALG.run=@run_bart; end
  ALGS{end+1}=ALG;
end

% mirt/fessler
if dim==3
  Js = 2:2:8;  % set of kernels widths to try. >8 doesn't help?
  for i=1:numel(Js)
    ALG=struct;
    ALG.algtype=5;
    ALG.name=sprintf('Fessler(J=%d)',Js(i));
    ALG.algopts=[];
    ALG.algopts.J = Js(i);      % kernel width in fine grid pts
    ALG.algopts.oversamp = 2.0;   % "R", oversampling ratio
    ALG.init=@init_fessler;
    if ty==1, ALG.run=@run_fessler_adj; else, ALG.run=@run_fessler; end
    ALGS{end+1}=ALG;
  end
end
% =============================================

% do all expts, storing all the data_out in results (can be big)...
results=run_algs(ALGS,x,y,z,data_in,N1,N2,N3);

% copy out results to various arrays, while also summarizing to text table...
disp('computing error metrics...')
txttbl{1} = sprintf('\n\n%15s %15s %15s %15s %15s','name','init_time(s)','run_time(s)','rel2err','therr');
X0=results{1}.data_out;         % truth
outputl1 = sum(abs(X0(:)));
outputl2 = norm(X0(:));
inputl1 = sum(abs(data_in(:)));  % appears in theoretical err estimate.
for j=1:length(ALGS)            % generate error values, text output
  X1=results{j}.data_out;
  rel2err(j)=norm(X0(:)-X1(:))/outputl2;
  rel1err(j)=norm(X0(:)-X1(:),1)/outputl1;
  maxerr(j)=max(abs(X0(:)-X1(:)));
  theoryerr(j)=maxerr(j)/inputl1;      % eps for which there's rigorous estimate
  init_times(j)=results{j}.init_time;
  run_times(j)=results{j}.run_time;
  total_times(j)=results{j}.tot_time;
  algtypes(j)=ALGS{j}.algtype;
  txttbl{j+1}=sprintf('%15s %15g %15g %15g %15g',ALGS{j}.name,init_times(j),run_times(j),rel2err(j),theoryerr(j));
end
for j=1:numel(txttbl), disp(txttbl{j}); end   % stdout

% these sets of indices for each code are useful for later plotting...
jf =find(algtypes==1);  % finufft
jn =find(algtypes==2);  % nfft
jnp=find(algtypes==2.25);  % nfft no pre (see silly encoding above)
jnf=find(algtypes==2.5);  % nfft full pre
jc =find(algtypes==3);  % cmcl
jb =find(algtypes==4);  % bart
jm =find(algtypes==5);  % mirt

clear i j x y z data_in X0 X1 results % kill the big stuff
save(outname)                        % save just error stats to disk (small)
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



%%%%%%%%%%%% UNIFORM WRAPPERS TO THE VARIOUS 3D NUFFT CODES %%%%%%%%%%%%%%%%%%%

function init_data=dummy_init(algopts,x,y,z,N1,N2,N3)
init_data=struct;


% finufft -------------------------------------------------------------------
function data_out=run_finufft3d1(algopts,x,y,z,data_in,N1,N2,N3,init_data)
[data_out ier]=finufft3d1(x,y,z,data_in,algopts.isign,algopts.eps,N1,N2,N3,algopts.opts);

function data_out=run_finufft3d2(algopts,x,y,z,data_in,N1,N2,N3,init_data)
[data_out ier]=finufft3d2(x,y,z,algopts.isign,algopts.eps,data_in,algopts.opts);

function data_out=run_finufft2d1(algopts,x,y,z,data_in,N1,N2,N3,init_data)
[data_out ier]=finufft2d1(x,y,data_in,algopts.isign,algopts.eps,N1,N2,algopts.opts);

function data_out=run_finufft2d2(algopts,x,y,z,data_in,N1,N2,N3,init_data)
[data_out ier]=finufft2d2(x,y,algopts.isign,algopts.eps,data_in,algopts.opts);

function data_out=run_finufft1d1(algopts,x,y,z,data_in,N1,N2,N3,init_data)
[data_out ier]=finufft1d1(x,data_in,algopts.isign,algopts.eps,N1,algopts.opts);

function data_out=run_finufft1d2(algopts,x,y,z,data_in,N1,N2,N3,init_data)
[data_out ier]=finufft1d2(x,algopts.isign,algopts.eps,data_in,algopts.opts);


% cmcl ---------------------------------------------------------------------
function data_out=run_nufft3d1(algopts,x,y,z,data_in,N1,N2,N3,init_data)
M=length(x);
[data_out ier]=nufft3d1(M,x,y,z,data_in,algopts.isign,algopts.eps,N1,N2,N3);
data_out = M*data_out;    % cmcl normalization

function data_out=run_nufft3d2(algopts,x,y,z,data_in,N1,N2,N3,init_data)
M=length(x);
[data_out ier]=nufft3d2(M,x,y,z,algopts.isign,algopts.eps,N1,N2,N3,data_in);

function data_out=run_nufft2d1(algopts,x,y,z,data_in,N1,N2,N3,init_data)
M=length(x);
[data_out ier]=nufft2d1(M,x,y,data_in,algopts.isign,algopts.eps,N1,N2);
data_out = M*data_out;    % cmcl normalization

function data_out=run_nufft2d2(algopts,x,y,z,data_in,N1,N2,N3,init_data)
M=length(x);
[data_out ier]=nufft2d2(M,x,y,algopts.isign,algopts.eps,N1,N2,data_in);

function data_out=run_nufft1d1(algopts,x,y,z,data_in,N1,N2,N3,init_data)
M=length(x);
[data_out ier]=nufft1d1(M,x,data_in,algopts.isign,algopts.eps,N1);
data_out = M*data_out;    % cmcl normalization

function data_out=run_nufft1d2(algopts,x,y,z,data_in,N1,N2,N3,init_data)
M=length(x);
[data_out ier]=nufft3d2(M,x,algopts.isign,algopts.eps,N1,data_in);


% nfft --------------------------------------------------------------------
function init_data=init_nfft(algopts,x,y,z,N1,N2,N3)
% dim is conveyed by N2, N3 not being 1.
M=length(x);
dim=1; N=N1;
if N2>1, dim=2; N=[N1;N2]; end
if N3>1, dim=3; N=[N1;N2;N3]; end    % #s of modes per dim. note N is a vector!
ticA=tic;   % set up the plan...
%plan=nfft(dim,N,M);  % default interface gives 13 digits, uses full kernel width
%... or use guru options, as commented in nfft-3.3.2/matlab/nfft/test_nfft3d.m :
n=2^(ceil(log2(max(N)))+1);      % lowest power of 2 that's at least 2N
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
ticA=tic;
% set NU pts, and use of nfft_init_guru...
sc = 1/(2*pi);   % NU pt scale factor so that periodic domain is [-.5,.5]^d
if dim==1
  plan=nfft(dim,N,M,n,algopts.m,flags,fftw_flag);
  nodes = x*sc;
elseif dim==2
  plan=nfft(dim,N,M,n,n,algopts.m,flags,fftw_flag);
  nodes = [x y]*sc;
else
  plan=nfft(dim,N,M,n,n,n,algopts.m,flags,fftw_flag);
  nodes = [x y z]*sc;
end
plan.x = nodes;   % nodes is M*d. set_x is a method, causes a RAM spike
fprintf('  Creating nfft plan (dim=%d) & setting nodes: %g s\n',dim,toc(ticA));

ticA=tic;
nfft_precompute_psi(plan);                     % precomputations
fprintf('  time for nfft_precompute_psi: %g s\n',toc(ticA));

init_data.plan=plan;

function X=run_nfft_adj(algopts,x,y,z,d,N1,N2,N3,init_data)
% type 1 ("adjoint"), any dim
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
% type 2 ("non-adjoint"), any dim
plan=init_data.plan;   % already carries NU locs (x y z args unused)
M=length(x);
plan.fhat=d(:);     % unwrap to single col
nfft_trafo(plan);   % type 2
X=plan.f;


% bart ---------------------------------------------------------------------
function X=run_bart_adj(algopts,x,y,z,d,N1,N2,N3,init_data)
% type 1 ("adjoint"), 3D only. Writes and reads from tempfile. See README.md
% Barnett 2/20/18.  algopts unused.  matches modeord=0 in finufft.
% untested for N1 neq N2 or neq N3.
weirdprefac = 1.00211;                    % expt prefactor best for BART!
prefac = sqrt(N1*N2*N3) / weirdprefac;
xyz = [x*(N1/2/pi) y*(N2/2/pi) z*(N3/2/pi)]';     % NU locs, size d*M
clear x y z
cmd=sprintf('nufft -a -d %d:%d:%d',N1,N2,N3);
X = prefac * bart(cmd,xyz,d.');               % NB strengths must be row vec

function X=run_bart(algopts,x,y,z,d,N1,N2,N3,init_data)
% type 2 ("fwd"), 3D only. Writes and reads from tempfile. See README.md
% Barnett 2/20/18.  algopts unused.  matches modeord=0 in finufft.
% untested for N1 neq N2 or neq N3.
weirdprefac = 1.00211;
prefac = sqrt(N1*N2*N3) / weirdprefac;
xyz = [x*(N1/2/pi) y*(N2/2/pi) z*(N3/2/pi)]';     % NU locs, size d*M
clear x y z
cmd=sprintf('nufft -d %d:%d:%d',N1,N2,N3);
X = prefac * bart(cmd,xyz,d);


% mirt ---------------------------------------------------------------------
function init_data=init_fessler(algopts,x,y,z,N1,N2,N3)
J=algopts.J;  % kernel width
N = [N1 N2 N3];
xyz = [x y z];     % NU locs in [-pi,pi]^d,  size M*d
init_data=nufft_init(xyz,N,[J,J,J],algopts.oversamp*N,N/2);

function X=run_fessler_adj(algopts,x,y,z,d,N1,N2,N3,init_data)
X=nufft_adj(d,init_data);

function X=run_fessler(algopts,x,y,z,d,N1,N2,N3,init_data)
X=nufft(d,init_data);
