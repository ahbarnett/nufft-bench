% check typ norm growth in 1d1, condition # of problem.
% Barnett 6/10/18
clear

M=1e3;
N=1e3;
x=pi*(2*rand(M,1)-1); c=randn(M,1)+1i*randn(M,1);
[f ier] = finufft1d1(x,c,+1,1e-6,N);
norm(f)/norm(c)

if M*N<=1e6   % test true cond # of (small) NUDFT prob...
  %M=1e3; N=1e3; A = exp(2i*pi*(-N/2:N/2-1)'*rand(1,M));
  A = exp(1i*(-N/2:N/2-1)'*x(:)');   % NUDFT matrix
  norm(f - A*c)/norm(f)    % NUFFT rel l2 error
  norm(A*c)/norm(c)
  [U S V] = svd(A);
  s = diag(S); [min(s), max(s), max(s)/min(s)]
  v =V(:,end);   % pick min sing vec & show it 
  [g] = finufft1d1(x,v,+1,1e-6,N);
  norm(g)/norm(v)
  norm(g - A*v)/norm(A*v)    % NUFFT rel l2 error
end

M=1e2; N=1e2; A = exp(1i*(-N/2:N/2-1)'*2*pi*rand(1,M)); min(svd(A))
c=randn(M,1)+1i*randn(M,1); norm(A*c)/norm(c)

% rand pts in half the period: note ill-cond even though M<N
M=80; N=100; A = exp(1i*(-N/2:N/2-1)'*pi*rand(1,M)); min(svd(A))
