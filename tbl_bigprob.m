% Compare speed and RAM: finufft vs NFFT, big 3d1 prob, fixed acc, meas RAM.
% User needs to set paths below.
% based on tbl_xeon_3d and speedup_vs_density2d1.
% Barnett 5/14/18
clear

multi=1;            % multithreaded test
dim=3;
N=256; M=5e7;   % data size. N=power of 2, so Nf same in finufft as NFFT
%N=512; M=5e8;   % data size
nudist=4;              % NU dist
tol = 1e-6;   % note NFFT m hits even # digits better than if odd.

setuppaths(multi);
opts.dt=0.1; memorygraph('start',opts);
ti=tic;
[x y z] = nudata(dim,nudist,M); M=numel(x);  % M may not be exactly as requested
d = randn(M,1) + 1i*randn(M,1);   % nu pt strengths (rand single-core)
fprintf('making data: %g s\n',toc(ti));
pause(0.3);                        % get the RAM usage...
[b,t,~,c] = memorygraph('get'); ramindata = b(end);
memorygraph('label','start more accurate FINUFFT');

maxthr = java.lang.Runtime.getRuntime().availableProcessors;
if multi, nthr = maxthr; else, nthr=1; end
maxNumCompThreads(nthr);

opts.nthreads=maxthr;  % use all threads for exact
opts.fftw=0;    % ESTIMATE
opts.chkbnds=0;
disp('higher-acc finufft (truth)...'); ti=tic;
[Fe ier] = finufft3d1(x,y,z,d,+1, 1e-2*tol, N,N,N,opts);
fprintf('  time for truth: %g s\n',toc(ti));

opts.nthreads=nthr;
opts.debug=2;   % 2 also prints spread info
memorygraph('label','start FINUFFT');
ti=tic;
[Ff ier] = finufft3d1(x,y,z,d,+1,tol,N,N,N,opts);   % do type 1
tf=toc(ti);
fprintf('  time for finufft: %g s\n',tf);
[b,t,~,c,ls,lt] = memorygraph('get'); MRf = max(b(t>lt(2)))-ramindata;
fprintf('FINUFFT max RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',MRf/1e9,MRf/M)
ef = norm(Ff(:)-Fe(:))/norm(Fe(:)); clear Ff  % rel l2 error
tbl(1,:) = [ef nan tf MRf];  % save to table: err, pre time, run time, RAM

m=ceil(log10(1/tol)/2);       % set NFFT width param; spread width is 2m+1
nf=2^(ceil(log2(N))+1);       % fine grid: lowest 2^k that's at least 2N

% NFFT no pre (apart from PHI)
flags=NFFT_OMP_BLOCKWISE_ADJOINT;
flags=bitor(FFT_OUT_OF_PLACE,flags);
flags=bitor(bitshift(uint32(1),11),flags); % NFFT_SORT_NODES
flags=bitor(PRE_PHI_HUT,flags);
fftw_flag = FFTW_ESTIMATE;         % careful...
fftw_flag=bitor(FFTW_MEASURE,1);  % FFTW_DESTROY_INPUT
xyz = [x y z]/(2*pi); clear x y z  % since want save RAM; don't count its time
memorygraph('label','start NFFT no-pre plan');
ti=tic;
plan=nfft(3,[N;N;N],M,nf,nf,nf,m,flags,fftw_flag); % use of nfft_init_guru
disp('nfft (no pre) plan done')
plan.x=xyz; clear xyz   % method sets M*3 node array in plan - a spike in RAM!
fprintf('nfft set_x done (it''s a method!), in %.3g s\n',toc(ti))
[b,t,~,c,ls,lt] = memorygraph('get'); MRc = max(b(t>lt(3)))-ramindata;
if isempty(MRc), MRc=0; end  % rescue if bad
fprintf('nfft set_x max RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',MRc/1e9,MRc/M)
nfft_precompute_psi(plan);  % here does nothing apart from pre_phi, cheap
tnp = toc(ti);
fprintf('nfft_precomp done in %.3g s\n',tnp)
plan.f=d; clear d    % be honest about RAM usage

memorygraph('label','start NFFT no-pre');
nfft_adjoint(plan);    % do type 1
tn=toc(ti); fprintf('tot time for nfft_adj(no pre): %g s\n',tn);
[b,t,~,c,ls,lt] = memorygraph('get'); MRn = max(b(t>lt(3)))-ramindata;
fprintf('NFFT(no pre) max RAM use above input data = %.3g GB (%.3g bytes/NUpt)\n',MRn/1e9,MRn/M)
Fn = plan.fhat;       % output
en = norm(Fn(:)-Fe(:))/norm(Fe(:)); clear Fn  % rel l2 error
tbl(2,:) = [en tnp tn max(MRc,MRn)];

% NFFT pre psi
flags=bitor(PRE_PSI,flags);        % PRE PSI
memorygraph('label','start NFFT pre psi plan');
ti=tic;
plan2=nfft(3,[N;N;N],M,nf,nf,nf,m,flags,fftw_flag); % use of nfft_init_guru
disp('nfft (pre psi) plan done')
plan2.x=plan.x; plan2.f=plan.f; clear plan  % to be honest about RAM
nfft_precompute_psi(plan2);                     % precomputations
tnprep=toc(ti);
fprintf('  time for nfft_precompute_psi: %g s\n',tnprep);
memorygraph('label','start NFFT pre psi');
ti=tic; nfft_adjoint(plan2);  % type 1
tnpre=toc(ti); fprintf('  time for nfft_adj(pre psi): %g s\n',tnpre);
Fnp = plan2.fhat;       % out
clear plan2
[b,t,~,c,ls,lt] = memorygraph('get'); MRnp = max(b(t>lt(5)))-ramindata;
enp = norm(Fnp(:)-Fe(:))/norm(Fe(:));  % rel l2 error
tbl(3,:) = [enp tnprep tnpre MRnp];

% show table
fprintf('\trel l2 err \tinit t (s) \trun t (s) \tRAM (GB)     bytes/NUpt\n')
for i=1:size(tbl,1)
  fprintf('%16.3g',[tbl(i,1:3) tbl(i,4)/1e9 tbl(i,4)/M]), fprintf('\n')
end

% plot
pause(0.3); [b,t,~,c,ls,lt] = memorygraph('get'); %memorygraph('plot');
memorygraph('done');
brel = b-ramindata;   % show relative to input data RAM need
figure; subplot(2,1,1); plot(t,brel/M,'.-'); xlabel('t (s)');
ylabel('RAM (bytes/NUpt)'); axis tight;
vline(lt,'r',ls); %text(0,0,'making data 40B/NUpt');
subplot(2,1,2); plot(t,c/100,'.-'); xlabel('t (s)');
ylabel('CPU usage (threads)'); axis tight; vline(lt,'r',ls); drawnow
fprintf('tot matlab mem: %.3g GB\n',total_matlab_memory/1e9) 

clear Fe Fnp; save bigprob.mat
