function [STATS] = WinBootCorFigure(STATS,condition_label_rw,Ylabel_rw)
% plot winsorized bootrapped correlations, p values and CIs across epoch.

figure; plot(STATS.xtimes,STATS.corr_results.(condition_label_rw).(Ylabel_rw).rw)
title([condition_label_rw, ' (Xdata) &',Ylabel_rw, ' (Ydata)']);
ylabel('winsorized & bootstrapped Pearsons r', 'fontsize', 15);
xlabel('time', 'fontsize', 15);

figure; plot(STATS.xtimes,STATS.corr_results.(condition_label_rw).(Ylabel_rw).p)
title([condition_label_rw, ' (Xdata) &',Ylabel_rw, ' (Ydata)']);
ylabel('p values', 'fontsize', 15);
xlabel('time', 'fontsize', 15);

figure; plot(STATS.xtimes,STATS.corr_results.(condition_label_rw).(Ylabel_rw).CI(1,:))
hold on
plot(STATS.xtimes,STATS.corr_results. Condition_label.Y_label.CI(2,:))
title([condition_label_rw, ' (Xdata) &',Ylabel_rw, ' (Ydata)']);
ylabel('CIs around rw', 'fontsize', 15);
xlabel('time', 'fontsize', 15);

end

