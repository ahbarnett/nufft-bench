% parallel scaling of finufft w/ # threads, just fig plot.
% Barnett 5/15/18

clear
dir = 'results/xeon/';  % results dir
for cc=1:2    % cases
  for type={'strong','weak'};
    ty=type{:};   % note is cell so need to convert to plain char array
    nam = [dir sprintf('%sscaling_3d_case%d',ty,cc)];  % file head
    l = load(nam);  % don't overwrite eg cc
    disp('        nthr        par efficiencies')
    disp([l.nthrs', l.ts(1,:) ./ (l.ts.*l.nthrs')])
    nbars = numel(l.nthrs);
    figure;
    if strcmp(ty,'strong')
      h = bar(1:nbars,l.ts(1,:) ./ (l.ts.*l.nthrs'));
      Mstr = sprintf('$M=%.3g$',l.baseM);   % title text giving M (base value)
    else                                    % weak
      h = bar(1:nbars,l.ts(1,:) ./ l.ts);
      Mstr = sprintf('$M=(%.3g)p$',l.baseM);
    end
    set(gca,'xticklabel',num2cellstr(l.nthrs));
    axis tight;
    title(sprintf('%s scaling: $\\epsilon=10^{%d}$, %s',ty,round(log10(l.eps)),Mstr),'interpreter','latex');
    xlabel('p (number of threads)');
    ylabel('parallel efficiency');
    v = vline(find(l.nthrs==12),'g--'); set(v,'linewidth',2);
    v = text(find(l.nthrs==12)+0.1,0.98,'# phys cores','color',[0 0.8 0]);
    legend('3D type-1', '3D type-2','location','southwest');
    set(gcf,'paperposition',[0 0 3 4]);
    print('-depsc2',[nam '.eps'])
  end
end
