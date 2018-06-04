% plot accuracy in various norms, vs requested tol. Barnett 6/2/18
% calls gen_acc

clear     % choose params...
tols     = logspace(-1,-14,14);  % requested tols
verb = 0;

M       = 1e7;    % # of NU pts (in all dims)
N       = M;    % # of modes (approx total, used in all dims)
[er2 eth] = gen_acc(N,M,tols,verb);
save results/xeon/accNM1e6
% smaller run
M0       = 1e3;    % # of NU pts (in all dims)
N0       = M0;    % # of modes (approx total, used in all dims)
[er20 eth0] = gen_acc(N0,M0,tols,verb);

if 0  % basic plot at one M=N
figure; loglog(tols,er2(1:3,:),'+-'); hold on;loglog(tols,tols,'k-')
title(sprintf('1D: M=%.3g N=%.3g',M,N)); axis tight
legend('type 1','type 2','type 3','req tol','location','northwest')
figure; loglog(tols,er2(4:6,:),'+-'); hold on;loglog(tols,tols,'k-')
title(sprintf('2D: M=%.3g N=%.3g',M,N)); axis tight
legend('type 1','type 2','type 3','req tol','location','northwest')
figure; loglog(tols,er2(7:9,:),'+-'); hold on;loglog(tols,tols,'k-')
title(sprintf('3D: M=%.3g N=%.3g',M,N)); axis tight
legend('type 1','type 2','type 3','req tol','location','northwest')
end

for t=1:3, l{t} = sprintf('type-%d M=1e%d',t,log10(M)); end
for t=1:3, l{t+3} = sprintf('type-%d M=1e%d',t,log10(M0)); end

figure; loglog(tols,er2(1:3,:),'+-'); hold on;
loglog(tols,er20(1:3,:),'.--'); hold on;loglog(tols,tols,'k-')
xlabel('tolerance $\varepsilon$','interpreter','latex')
ylabel('$\| \tilde g - g\|_2/\|g\|_2$','interpreter','latex')
title('(a) 1D: rel. $\ell_2$ error vs req. tolerance','interpreter','latex'); axis tight
legend({l{:},'$\varepsilon$'},'location','northwest','interpreter','latex')
set(gcf,'paperposition',[0 0 3.5 3.5]);
print -depsc2 results/xeon/acc1D.eps

figure; loglog(tols,er2(4:6,:),'+-'); hold on;
loglog(tols,er20(4:6,:),'.--'); hold on;loglog(tols,tols,'k-')
xlabel('tolerance $\varepsilon$','interpreter','latex')
ylabel('$\| \tilde g - g\|_2/\|g\|_2$','interpreter','latex')
title('(b) 2D: rel. $\ell_2$ error vs req. tolerance','interpreter','latex'); axis tight
legend({l{:},'$\varepsilon$'},'location','northwest','interpreter','latex')
set(gcf,'paperposition',[0 0 3.5 3.5]);
print -depsc2 results/xeon/acc2D.eps

figure; loglog(tols,er2(7:9,:),'+-'); hold on;
loglog(tols,er20(7:9,:),'.--'); hold on;loglog(tols,tols,'k-')
xlabel('tolerance $\varepsilon$','interpreter','latex')
ylabel('$\| \tilde g - g\|_2/\|g\|_2$','interpreter','latex')
title('(c) 3D: rel. $\ell_2$ error vs req. tolerance','interpreter','latex'); axis tight
legend({l{:},'$\varepsilon$'},'location','northwest','interpreter','latex')
set(gcf,'paperposition',[0 0 3.5 3.5]);
print -depsc2 results/xeon/acc3D.eps
