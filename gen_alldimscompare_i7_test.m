% run all timing expts, writing to .mat files
% barnett 3/13/18

cpu='i7'; nthreads=8;
o = []; o.memcpu = 1;

%multithreaded=1; M=1e8; NN=1e7;  % multi -------- bigger expt
multithreaded=1; M=1e7; NN=1e6;  % multi -------- bigger expt

for ty = 1 %[1 2]
  for nudist = 0 %[0 4]
    for dim = 2 %[2 3]
      N=round(NN^(1/dim)); if mod(N,2)==1, N=N+1; end   % modes per dim, eve
      o.nfftpres = 2; if dim==3 & M>1e7, o.nfftpres=1; end
      outname=sprintf('results/%s/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',cpu,dim,ty,nudist,N,M,nthreads,o.nfftpres);
      benchallcodes(ty,dim,N,M,nudist,multithreaded,outname,o);
    end
  end
end

stop

multithreaded=0; nthreads = 1; M=1e7; NN=1e6;  % single -------- smaller expt

for ty = [1 2]
  for nudist = [0 4]
    for dim = [2 3]
      N=round(NN^(1/dim)); if mod(N,2)==1, N=N+1; end   % modes per dim, even
      outname=sprintf('results/%s/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',cpu,dim,ty,nudist,N,M,nthreads,o.nfftpres);
      benchallcodes(ty,dim,N,M,nudist,multithreaded,outname,o);
    end
  end
end
