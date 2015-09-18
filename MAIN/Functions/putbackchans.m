function putbackchans()

global CURCLICK
global tmpEEG
global CONDFILES_SUBS
global PATHTOFILES

% load again to be safe
tmpEEG = pop_loadset('filename',CONDFILES_SUBS{1}{1},'filepath',PATHTOFILES{1});
tmpEEG = eeg_checkset(tmpEEG);
disp(['working on channel selections for file: ' CONDFILES_SUBS{1}{1}]);

statslab_topoplot([],tmpEEG.chanlocs, 'style', 'blank', 'drawaxis', 'on', 'electrodes','labelpoint', 'plotrad', [], 'chaninfo', tmpEEG, 'nosedir' ,'+Y');
f=gcf;
set(f,'Color', [1 1 1]);

%%%% not finished, needs testing
curind=cellfun(@str2num, CURCLICK);
chanlabs={EEG.chanlocs(curind).labels};
for i=1:length(curind);
    try
        findobj('String',chanlabs{i});
        set(gco, 'Color', 'green', 'FontSize',13, 'FontWeight','bold')
    catch
        findobj('String',[chanlabs{i}, ' ']);
        set(gco, 'Color', 'green', 'FontSize',13, 'FontWeight','bold')
        
    end
    
end


