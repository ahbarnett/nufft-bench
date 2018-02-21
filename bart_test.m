% test BART nufft, conventions. This is 3D only, apparently.
% Uses: memorygraph
% Barnett 2/20/18, after new grid.c emailed to me fixing the wrapping bug.

run ~/numerics/nufft/bart-0.4.02/startup        % needed for BART
addpath ~/numerics/finufft/matlab
format long g
clear
%bart('nufft -h')  % get help
prefac = 1.00211;     % random prefactor best for BART!

if 0  % warm up

% type 1
%xyz = [1;3;-2];  % one src pt
xyz = [10;0;0];  % one src pt, at n/2
c = 1-.5i;   % strength

n=20;  % even for now
N=n^3;
cmd=sprintf('nufft -a -d %d:%d:%d',n,n,n);
F1=bart(cmd,xyz,c);

k=(-n/2:n/2-1)';   % col vec of output freq indices
%sqrt(N)*F1(:,1,1) - c*exp(pi*2i*k*xyz(1)/n)   % also sqrt(N-4) closer factor!
%sqrt(N-4)*F1(:,1,1) - c*exp(pi*2i*k*xyz(1)/n)   % also sqrt(N-4) closer factor!

sc=2*pi/n;
F2 = finufft3d1(sc*xyz(1),sc*xyz(2),sc*xyz(3),c,+1,1e-12,n,n,n);

%max(abs(sqrt(N)*F1(:,1,1) - F2(:,1,1)))
%max(abs(sqrt(N)*F1(:) - F2(:)))
F1 = F1*(sqrt(N)/prefac);
[F2(:,1,1), F1(:,1,1)]

fprintf('rel l2 err = %.3g\n',norm(F2(:)-F1(:))/norm(F2(:)))
fprintf('rel max err = %.3g\n',max(abs(F2(:)-F1(:)))/max(abs(F2(:))))

% suggests xyz coords divided by n to change to CMCL & finufft defns.
% isign = +1
end

if 1
% type 1
n=50;  % U size in each dim
M=1e6;
N = n^3;   % total output pts
c = randn(1,M)+1i*randn(1,M);
x = pi*(2*rand(1,M)-1);   % std nufft scaling for NU
y = pi*(2*rand(1,M)-1);
z = pi*(2*rand(1,M)-1);
s = 1.0; x=x*s;y=y*s;z=z*s;   % shows the prob is at the boundary of the NU box!
xyz = [x;y;z]*(n/(2*pi));  % BART relative scaling, is src box [-n/2,n/2]
disp('rand data done. starting bart...')

%o.dt=0.1;memorygraph('start',o);

bartcmd = sprintf('nufft -a -d %d:%d:%d',n,n,n);  % -a = adj = type1

tic; F1 = bart(bartcmd,xyz,c);
fprintf('bart done in %.3g s\n',toc)
%memorygraph('label','bart done');

disp('starting finufft...')
tol=1e-6;
o=[]; o.modeord=0;   % guessing: -n/2 to n/2
tic; F2 = finufft3d1(x,y,z,c,+1,tol,n,n,n,o);
fprintf('finufft done in %.3g s\n',toc)
%memorygraph('label','finufft done');
%memorygraph('plot'); memorygraph('done');

F1 = (sqrt(N)/prefac)*F1;   % rescale BART

%F1(:,1,1)./F2(:,1,1)  % look for the problem

fprintf('rel l2 err = %.3g\n',norm(F2(:)-F1(:))/norm(F2(:)))
fprintf('rel max err = %.3g\n',max(abs(F2(:)-F1(:)))/max(abs(F2(:))))
end

