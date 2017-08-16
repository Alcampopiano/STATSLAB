function [STATS] = WinBootCorFigure(STATS,condition_label_rw,Ylabel_rw)

% plot winsorized bootrapped correlations, p values and CIs across epoch.
% 
% Inputs:
% 
% ***condition_label_rw***
% 
% (string) - Specify the label for the condition or contrast. This string corresponds to the string used to compute correlations in previous module. ***end***
% 
% ***Ylabel_rw***
% 
% (string) - Specify the label for the Yaxis. This string corresponds to the Ylabel string used to compute correlations in previous module.
% 
% WinbootCorFigure from the commandline:
% 
% [STATS]=WinBootCorFigure('STATS_someanalysis.mat', 'my_condition_label', 'myYlabel')
% 
% ***end***


if isempty(STATS)
    [fnamestats]=uigetfile('*.mat','Select the STATS structure for this analysis');
    load(fnamestats);
else
    load(STATS);
end

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
plot(STATS.xtimes,STATS.corr_results.(condition_label_rw).(Ylabel_rw).CI(2,:))
title([condition_label_rw, ' (Xdata) &',Ylabel_rw, ' (Ydata)']);
ylabel('CIs around rw', 'fontsize', 15);
xlabel('time', 'fontsize', 15);

end

