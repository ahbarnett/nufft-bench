% Barnett trying Lederman's code for 1D classical PSWF. 6/13/18
% Also comparing approximate forms for PSWF.

% needs Roy Lederman's package: https://github.com/lederman/Prol_1D
addpath /home/alex/numerics/lederman/Prol_1D/src/matlab
clear

% adjust this:
c=10;  % freq param (beta), seems good up to 40 where roundoff hits edges

% Roy's set-up
matdim=200; minEigenvalRatio = 10^-40;
[prolate_dat, iserr] = prolate_1d_crea(c,matdim, minEigenvalRatio);

%figure; semilogy((abs(prolate_dat.nu))); ylim([10^-30,3])
%title('magnitude of eigenvalues (not scaled)')
    
xx=linspace(-1,1,1e4+1)';  h = xx(2)-xx(1);   % note one xx is at 0
index=0; [v,dv] = prolate_1d_ev(prolate_dat, index, xx);  % Roy eval PSWF
fprintf('L2 norm = %.3g, size at ends = %.3g\n',sqrt(h*sum(v.^2)),v(end))
%fprintf('1-nu_0 = %.3g\n',1-prolate_dat.nu(1))  % exp(-2c), too small to see
v = v/max(v);   % central value set to 1
%figure; plot(xx,v); title(sprintf('c=%.3g',c))
es = exp(c*(sqrt(1-xx.^2)-1));
esc = es./(1-xx.^2).^(3/8);    % -3/8 seems to be the best power!
esc(isinf(esc))=0;    % tidy endpt vals
fprintf('rel L2 err of esc(3/8) from pswf = %.g\n',norm(esc-v)/norm(v))
kb = besseli(0,c*sqrt(1-xx.^2)) / besseli(0,c);
kbasym = es./(1-xx.^2).^(1/4);
kbsl = kbasym./sqrt(1+sqrt(1-xx.^2)) * sqrt(2);    % Slepian, & match at origin
kbslh = kb./sqrt(1+sqrt(1-xx.^2)) * sqrt(2);  % hybrid tweak to Slepian
fprintf('rel L2 err of kbslh from pswf = %.g\n',norm(kbslh-v)/norm(v))
figure; semilogy(xx,v,'k-');
hold on; semilogy(xx,[es,esc,kb,kbsl,kbslh]);
legend('PSWF','ES','ES -3/8 pow', 'KB','KB Slep','KB Slep hybrid');
title(sprintf('c=%.3g',c))
figure; plot(xx,[es./v, esc./v, kb./v, kbsl./v, kbslh./v]); hline(1);
legend('ES/PSFW','ES -3/8 pow / PSWF', 'KB / PSWF','KB Slep / PSWF','KB Slep hybrid / PSWF');
title(sprintf('ratios at c=%.3g',c))
axis([.5 1 .5 1.3])

% conclusion: Slepian's asymp is better than my -3/8 power guess,
% but a slight tweak on Slepian is even better! (kbslh).

% none are exponentially close in L2 norm. Just small pointwise relative error,
% eg 1/beta sized.
