% xeon big 3D expts for NUFFT.
% Barnett 2/20/18 updated for memorygraph.

clear
multi=1;            % multithreaded test
nthreads=24;        % twice the physical # cores. for ccb desktop
setuppaths(multi);
opts.dt=0.2;memorygraph('start',opts);

dim=3; N=128; M=1e8;   % data size
nudist=4;              % NU dist

t=tic;
[x y z] = nudata(dim,nudist,M); M=numel(x);  % M may not be exactly as requested
d = randn(M,1) + 1i*randn(M,1);   % nu pt strengths (rand single-core)
fprintf('making data: %g s\n',toc(t));
pause(0.3);                        % get the RAM usage...
[b,t,~,c] = memorygraph('get'); ramindata = b(end);
memorygraph('label','start FINUFFT');

maxNumCompThreads(nthreads);  % makes some matlab things run faster; FINUFFT uses
opts.nthreads=nthreads;  % use all threads
opts.fftw=0;    % ESTIMATE
opts.debug=2;   % 2 also prints spread info
opts.chkbnds=0;  % saves 6 sec @ M=1e9
% Sorting time was 56s @ M=1e9 (single-core), now 7 sec (multi-core). @4x4x16.
t=tic;
[Ff ier] = finufft3d1(x,y,z,d,+1,1e-6,N,N,N,opts);   % do type 1
tf=toc(t);
fprintf('  time for finufft: %g s\n',tf);
[b,t,~,c] = memorygraph('get'); MRf = max(b)-ramindata;
fprintf('FINUFFT max RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',MRf/1e9,MRf/M)

m=3;    % NFFT: kernel half-wid (m=3 for tol 1e-6)
n=2^(ceil(log2(N))+1);                % lowest power of 2 that's at least 2N
flags=NFFT_OMP_BLOCKWISE_ADJOINT;         % advanced NFFT flags
flags=bitor(FFT_OUT_OF_PLACE,flags);
flags=bitor(bitshift(uint32(1),11),flags); % NFFT_SORT_NODES
flags=bitor(PRE_PHI_HUT,flags);
fftw_flag = FFTW_ESTIMATE;         % careful...
fftw_flag = bitor(fftw_flag, 1);    % FFTW_DESTROY_INPUT
plan=nfft(3,[N;N;N],M,n,n,n,m,flags,fftw_flag); %   choose no pre
disp('nfft plan done')
xyz = [x y z]/(2*pi); clear x y z  % since want save RAM
memorygraph('label','start NFFT plan setx');
plan.x=xyz; clear xyz   % method sets M*3 node array in plan - a spike in RAM!
disp('nfft set_x done (it''s a method!)')
[b,t,~,c,ls,lt] = memorygraph('get'); MRc = max(b(t>lt(1)))-ramindata;
fprintf('nfft set_x max RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',MRc/1e9,MRc/M)
t=tic; nfft_precompute_psi(plan);       % precomputations (here does nothing)
fprintf('nfft_precomp done in %.3g s\n',toc(t))
plan.f=d;  clear d     % be honest about RAM usage
memorygraph('label','start NFFT no-pre');
nfft_adjoint(plan);    % do type 1
tn=toc(t); fprintf('tot time for nfft_adj(no pre): %g s\n',tn);
[b,t,~,c] = memorygraph('get'); MRn = max(b(t>lt(2)))-ramindata;
fprintf('NFFT(no pre) max RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',MRn/1e9,MRn/M)
Fn = plan.fhat;       % output
%nfft_finalize(plan); % fails!

if 0  % also include NFFT pre psi?
flags=bitor(PRE_PSI,flags);        % PRE PSI
plan2=nfft(3,[N;N;N],M,n,n,n,m,flags,fftw_flag); % use of nfft_init_guru
plan2.x=plan.x; plan2.f=plan.f; clear plan  % be honest about RAM
t=tic; nfft_precompute_psi(plan2);                     % precomputations
tnpre=toc(t);
fprintf('  time for nfft_precompute_psi: %g s\n',tnpre);
t=tic; nfft_adjoint(plan2);  % type 1
tnp=toc(t); fprintf('  time for nfft_adj(pre psi): %g s\n',tnp);
Fnp = plan2.fhat;       % out
clear plan2
else
  clear plan
end

% plot
pause(0.3); [b,t,~,c,ls,lt] = memorygraph('get'); %memorygraph('plot');
memorygraph('done');
brel = b-ramindata;   % show relative to input data RAM need
figure; subplot(2,1,1); plot(t,brel/M,'.-'); xlabel('t (s)');
ylabel('RAM (bytes/NUpt)'); axis tight;
vline(lt,'r',ls); text(0,0,'making data 40B/NUpt');
subplot(2,1,2); plot(t,c/100,'.-'); xlabel('t (s)');
ylabel('CPU usage (threads)'); axis tight; vline(lt,'r',ls); drawnow
total_matlab_memory


fprintf('  rel2diff (f vs n) = %.3g\n', norm(Fn(:)-Ff(:))/norm(Ff(:)))
%fprintf('  rel2diff (np vs n) = %.3g\n\n', norm(Fn(:)-Fnp(:))/norm(Fn(:)))
fprintf('  time NFFT(no pre) / FINUFFT =      %.3g\n',tn/tf)
%fprintf('  time NFFT(pre psi) / FINUFFT =     %.3g\n',tnp/tf)
%fprintf('  time NFFT(pre psi tot) / FINUFFT = %.3g\n',(tnpre+tnp)/tf)

%results nudist=1: M=5e8: both us & nff nopre use 22B/NUpt
% we're 4.1x faster.
