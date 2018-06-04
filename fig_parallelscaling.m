% parallel scaling of finufft w/ # threads, just fig plot.
% Barnett 5/15/18

clear
dir = 'results/xeon/';  % results dir
subfigs = 'cdab';     % ordering of plots in figure row.
subfig=1;
for cc=1:2    % cases
  for type={'strong','weak'};
    ty=type{:};   % note is cell so need to convert to plain char array
    nam = [dir sprintf('%sscaling_3d_case%d',ty,cc)];  % file head
    l = load(nam);  % don't overwrite eg cc
    disp('        nthr        par efficiencies')
    disp([l.nthrs', l.ts(1,:) ./ (l.ts.*l.nthrs')])
    nbars = numel(l.nthrs);
    figure;
    Mpow10 = round(log10(l.baseM));
    Mnum = sprintf('%d\\times 10^{%d}',l.baseM/10^Mpow10,Mpow10); % sci notation
    if strcmp(ty,'strong')
      h = bar(1:nbars,l.ts(1,:) ./ (l.ts.*l.nthrs'));
      Mstr = sprintf('$M{=}%s$',Mnum);   % title text giving M (base value)
    else                                    % weak
      h = bar(1:nbars,l.ts(1,:) ./ l.ts);
      Mstr = sprintf('$M{=}(%s)p$',Mnum);
    end
    set(gca,'xticklabel',num2cellstr(l.nthrs));
    axis tight;
    title(sprintf('%s scaling: $\\epsilon{=}10^{%d}$, %s',ty,round(log10(l.eps)),Mstr),'interpreter','latex');
    xlabel('$p$ (number of threads)','interpreter','latex');
    ylabel('$\eta$ (parallel efficiency)','interpreter','latex');
    v = vline(find(l.nthrs==12),'g--'); set(v,'linewidth',2);
    text(find(l.nthrs==12)+0.1,0.98,'# phys','color',[0 0.8 0]);
    text(find(l.nthrs==12)+0.1,0.92,' cores','color',[0 0.8 0]);
    text(7,0.85,['(' subfigs(subfig) ')'],'color',[0 0 0]);
    legend('3D type-1', '3D type-2','location','southwest');
    set(gcf,'paperposition',[0 0 2.5 3.5]);
    print('-depsc2',[nam '.eps'])
    subfig=subfig+1;
  end
end
