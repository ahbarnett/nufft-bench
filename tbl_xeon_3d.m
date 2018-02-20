% xeon big 3D expts for NUFFT.
% Barnett 2/2/18

clear
multi=1; setuppaths(multi);
opts.dt=0.2;memorygraph('start',opts);

N=256; M=1e9;
nudist=4; dim=3;
t=tic;
[x y z] = nudata(dim,nudist,M);
d = randn(M,1) + 1i*randn(M,1);   % nu pt strengths (rand single-core)
fprintf('making data: %g s\n',toc(t));
pause(0.3);  % check the RAM usage
[b,t,~,c] = memorygraph('get'); ramindata = b(end); % 40B/NUpt = 5 doubles
tstamp0 = numel(b); % save how may timepts in MG up to now

maxNumCompThreads(24);  % makes some matlab things run faster; FINUFFT uses this
opts.nthreads=24;  % use all threads
opts.fftw=0;  % ESTIMATE
opts.debug=2;   % 2 also prints spread info
opts.chkbnds=0;  % saves 6 sec @ M=1e9.  Sorting time 56s @ M=1e9 (single-core)
% sorting: 28s @ M=5e8. 4.7s @ M=1e8.
t=tic; [Ff ier] = finufft3d1(x,y,z,d,+1,1e-6,N,N,N,opts);
tf=toc(t);
fprintf('  time for finufft: %g s\n',tf);
[b,t,~,c] = memorygraph('get'); MRf = max(b)-ramindata;
fprintf('FINUFFT max RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',MRf/1e9,MRf/M)
tstamp1 = numel(b); % save how may timepts in MG up to now

m=3;  % NFFT: kernel half-wid
n=2^(ceil(log2(N))+1);                % lowest power of 2 that's at least 2N
flags=NFFT_OMP_BLOCKWISE_ADJOINT;         % advanced NFFT flags
flags=bitor(FFT_OUT_OF_PLACE,flags);
flags=bitor(bitshift(uint32(1),11),flags); % NFFT_SORT_NODES
flags=bitor(PRE_PHI_HUT,flags);
fftw_flag = FFTW_ESTIMATE;         % careful...
fftw_flag = bitor(fftw_flag, 1);    % FFTW_DESTROY_INPUT
plan=nfft(3,[N;N;N],M,n,n,n,m,flags,fftw_flag); %   no pre
disp('nfft plan done')
xyz = [x y z]/(2*pi); clear x y z  % since pl
plan.x=xyz; clear xyz   % set M*3 array of nodes in plan, min RAM
disp('nfft set_x done (it''s a method!)')
[b,t,~,c] = memorygraph('get'); MRc = max(b(tstamp1:end))-ramindata;
tstamp2 = numel(b); % save how may timepts in MG up to now
fprintf('nfft set_x max RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',MRc/1e9,MRc/M)
t=tic; nfft_precompute_psi(plan);       % precomputations (here does nothing)
fprintf('nfft_precomp done in %.3g\n',toc(t))
plan.f=d;  clear d     % be honest about RAM usage
nfft_adjoint(plan);  % type 1
tn=toc(t); fprintf('tot time for nfft_adj(no pre): %g s\n',tn);
[b,t,~,c] = memorygraph('get'); MRn = max(b(tstamp2:end))-ramindata;
tstamp3 = numel(b); % save how may timepts in MG up to now
fprintf('NFFT(no pre) max RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',MRn/1e9,MRn/M)
Fn = plan.fhat;       % out
%nfft_finalize(plan); % fails!

if 0
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

pause(0.3); [b,t,~,c] = memorygraph('get');
memorygraph('done');
b = b-ramindata;   % work relative to input data RAM need
figure; subplot(2,1,1); plot(t,b/M,'.-'); xlabel('t (s)');
ylabel('RAM (bytes/NUpt)'); axis tight; vline(opts.dt*[tstamp0 tstamp1 tstamp2]);
text(0,0,'making data 40B/NUpt');
text(opts.dt*tstamp0,0,'start FINUFFT');
text(opts.dt*tstamp1,0,'start NFFT plan set_x');
text(opts.dt*tstamp2,0,'start NFFT no-pre');
subplot(2,1,2); plot(t,c/100,'.-'); xlabel('t (s)');
ylabel('CPU usage (threads)'); axis tight; drawnow
total_matlab_memory

fprintf('  rel2diff (f vs n) = %.3g\n', norm(Fn(:)-Ff(:))/norm(Ff(:)))
%fprintf('  rel2diff (np vs n) = %.3g\n\n', norm(Fn(:)-Fnp(:))/norm(Fn(:)))
fprintf('  time NFFT(no pre) / FINUFFT =      %.3g\n',tn/tf)
%fprintf('  time NFFT(pre psi) / FINUFFT =     %.3g\n',tnp/tf)
%fprintf('  time NFFT(pre psi tot) / FINUFFT = %.3g\n',(tnpre+tnp)/tf)

%results nudist=1: M=5e8: both us & nff nopre use 22B/NUpt
% we're 4.1x faster.
