% run all timing comparison expts, writing to .mat files.
% This code is bastardized to do required expts.
% barnett 4/27/18. 5/15/18

cpu='xeon';
o = []; o.memcpu=1;   % whether to collect ram info

for multithreaded=[0 1]

  if multithreaded, nthreads = 24; M=1e8; NN=1e7;  % multi -------- bigger expt
  else, nthreads = 1; M=1e7; NN=1e6;  % single -------- smaller expt
  end

  for dim = [1 2 3]
    if dim<3, o.reps=3; else, o.reps=1; end
    if dim==1, nudists=[0]; else nudists=[0 4]; end
    for ty = [1 2]
      for nudist = nudists
        N=round(NN^(1/dim)); if mod(N,2)==1, N=N+1; end   % modes per dim, even
        o.nfftpres = 2; if dim==3 & M>1e7, o.nfftpres=1; end
        outname=sprintf('results/%s/%dd%d_nudist%d_N%d_M%d_%dthr_nfftpres%d',cpu,dim,ty,nudist,N,M,nthreads,o.nfftpres);
        benchallcodes(ty,dim,N,M,nudist,multithreaded,outname,o);
      end
    end
  end
end
