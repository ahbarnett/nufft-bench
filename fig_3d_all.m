% Top-level driver script for 3D expts for NUFFT.
% Barnett 1/28/18

if 0
clear
calc=1;
expt='xeon'; multi=1;
nfftpres = 1;  % note: 0 no extra RAM, 1 needs 0.6kB/NUpt, 2 needs 35kB/NUpt!
               % (these are for m=6 kernel)
%N=64; M=1e7;
N=128; M=1e8;
fig_3d(expt,1,multi,N,M,nfftpres,calc)
fig_3d(expt,2,multi,N,M,nfftpres,calc)
end

% ------------------------------------------------------------------
% single entry in table - how meas RAM?
clear
memorygraph('start',struct('dt',0.1));

dim=3; N=128; M=1e8;
nudist=4;
[x y z] = nudata(dim,nudist,M); M=numel(x); % in case changed
d = randn(M,1) + 1i*randn(M,1);   % nu pt strengths
[b,t,~,c] = memorygraph('get');    % check the RAM usage
ramindata = b(end);

tol=1e-6; %1e-12;
opts.nthreads=8;  % FINUFFT: all threads
opts.fftw=0;  % ESTIMATE
opts.chkbnds=0;  % saves 6 sec @ M=1e9
opts.debug=2;
t=tic; [Ff ier] = finufft3d1(x,y,z,d,+1,tol,N,N,N,opts);
tf=toc(t);
fprintf('  time for finufft: %g s\n',tf);
[b,t,~,c] = memorygraph('get'); MRf = max(b)-ramindata;
fprintf('FINUFFT max RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',MRf/1e9,MRf/M)

m=3; %6;  % NFFT: kernel half-wid
n=2^(ceil(log2(N))+1);                % lowest power of 2 that's at least 2N
flags=NFFT_OMP_BLOCKWISE_ADJOINT;         % advanced NFFT flags
flags=bitor(FFT_OUT_OF_PLACE,flags);
flags=bitor(bitshift(uint32(1),11),flags); % NFFT_SORT_NODES
flags=bitor(PRE_PHI_HUT,flags);
fftw_flag = FFTW_ESTIMATE;         % careful...
fftw_flag = bitor(fftw_flag, 1);    % FFTW_DESTROY_INPUT
plan=nfft(3,[N;N;N],M,n,n,n,m,flags,fftw_flag); %   no pre
plan.x=[x y z]/(2*pi); clear x y z  % set M*3 array of nodes in plan, honest RAM
t=tic; nfft_precompute_psi(plan);                     % precomputations
plan.f=d;  clear d     % be honest about RAM usage
nfft_adjoint(plan);  % type 1
tn=toc(t); fprintf('  time for nfft_adj(no pre): %g s\n',tn);
[b,t,~,c] = memorygraph('get'); MRn = max(b)-ramindata;
fprintf('NFFT(no pre) max RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',MRn/1e9,MRn/M)
Fn = plan.fhat;       % out
%nfft_finalize(plan); % fails

flags=bitor(PRE_PSI,flags);        % PRE PSI
plan2=nfft(3,[N;N;N],M,n,n,n,m,flags,fftw_flag); % use of nfft_init_guru
plan2.x=plan.x; plan2.f=plan.f; clear plan  % be honest about RAM
t=tic; nfft_precompute_psi(plan2);                     % precomputations
tnpre=toc(t);
fprintf('  time for nfft_precompute_psi: %g s\n',tnpre);
t=tic; nfft_adjoint(plan2);  % type 1
tnp=toc(t); fprintf('  time for nfft_adj(pre psi): %g s\n',tnp);
[b,t,~,c] = memorygraph('get'); MRnp = max(b)-ramindata;
fprintf('NFFT(pre psi) max RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',MRnp/1e9,MRnp/M)
Fnp = plan2.fhat;       % out
clear plan2

pause(0.3); [b,t,~,c] = memorygraph('get');
memorygraph('done');
b = b-ramindata;   % work relative to input data RAM need
%fprintf('overall peak RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',max(b)/1e9,max(b)/M)
figure; subplot(2,1,1); plot(t,b/M,'.-'); xlabel('t (s)');
ylabel('RAM (bytes/NUpt)'); axis tight; subplot(2,1,2);
plot(t,c/100,'.-'); xlabel('t (s)'); ylabel('CPU usage (threads)');
axis tight; drawnow
%total_matlab_memory  % obsolete

fprintf('  rel2diff (f vs n) = %.3g\n', norm(Fn(:)-Ff(:))/norm(Ff(:)))
fprintf('  rel2diff (np vs n) = %.3g\n\n', norm(Fn(:)-Fnp(:))/norm(Fn(:)))
fprintf('  time NFFT(no pre) / FINUFFT =      %.3g\n',tn/tf)
fprintf('  time NFFT(pre psi) / FINUFFT =     %.3g\n',tnp/tf)
fprintf('  time NFFT(pre psi tot) / FINUFFT = %.3g\n',(tnpre+tnp)/tf)

