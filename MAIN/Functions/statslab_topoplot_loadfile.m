function statslab_topoplot_loadfile()

global tmpEEG

statslab_topoplot([],tmpEEG.chanlocs, 'style', 'blank', 'drawaxis', 'on', 'electrodes','labelpoint', 'plotrad', [], 'chaninfo', tmpEEG, 'nosedir' ,'+Y');
f=gcf;
% set(f,'units','normalized','Position',[.05 .10 .9 .8], 'Color', [.66 .76 1]);
%set(f,'Color', [.66 .76 1]); % EEGLAB default
set(f,'Color', [1 1 1]); 
end

