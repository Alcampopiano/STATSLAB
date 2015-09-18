function [okayhit, chan_choices]=chanpick_topo(condfiles_subs,pathtofiles,numconds)

clearvars -global tmpEEG CONDFILES_SUBS PATHTOFILES CHANCHOICES CURCLICK

global tmpEEG
global CONDFILES_SUBS
global PATHTOFILES
CONDFILES_SUBS=condfiles_subs;
PATHTOFILES=pathtofiles;
global CHANCHOICES;

% temporary load
tmpEEG = pop_loadset('filename',CONDFILES_SUBS{1}{1},'filepath',PATHTOFILES{1});
tmpEEG = eeg_checkset(tmpEEG);

%size of largest number of subjects gives size of grid
tabsize=max(cell2mat(cellfun(@length,condfiles_subs,'un',0)));

% the first part of evstring1 just deletes a weird uicontrol object that supergui places
% on the figure that appears in blue. If this is not the 6th object in the
% figure, this code will not work. Better to find the object by name
% somehow.
evstring1='b=findobj(gcf); delete(b(6)); set(gcf,''Color'', [1 1 1]); statslab_topoplot_loadfile;';
evstring2='f=gcf; close(f); putbackchans;';

CHANCHOICES=cell(tabsize,numconds*2);
j=1;
for i=1:numconds;
    [rowfname colfname]=size(condfiles_subs{i});
    for q=1:rowfname
        CHANCHOICES{q,j}=condfiles_subs{i}{q};
        CHANCHOICES{q,j+1}='';
    end
    j=j+2;
end


[res, jnk, okayhit]=inputgui( ...
    'geom', ...
    {...
    {6 22 [0 0] [1 1]} ... % this is just to control the size of the GUI {6 22 [0 0] [1 1]} ...
    {6 22 [0.05 0] [2 1]} ... % load config
    {6 22 [3.95 0] [2 1]} ... %2 saveas

    }, ...
    'uilist', ...
    {...
    {'Style', 'text', 'tag','txt_ccfp','string',blanks(30)} ... %1 this is just to control the size of the GU
    {'Style', 'pushbutton', 'string', 'Load Parameters', ...
    'callback', ...
    ['[ParamName, ParamPath] = uigetfile(''*.mat'',''choose channel selection file:'',''*.mat'',''multiselect'',''off'');', ...
    'load(fullfile(ParamPath, ParamName),''-mat'');',evstring2]} ... %2
    {'Style', 'pushbutton','string','Save channel selection file', ...
    'callback',['[ParamName, ParamPath]=uiputfile(''*.*'',''channel selection file'');', ...
    'global CURCLICK;', ...
    'save(fullfile(ParamPath, ParamName),''CURCLICK'');']}, ...
    }, ...
    'title', 'STATSLAB -- statslab()',...%, ...
    'eval',evstring1 ...
    );

% this just converts char arrays to cells so that users don't have to use
% curly braces and single quotes for channel labels, but instead can simply
% delimit with spaces only when needed more than one channel.
j=2;
for i=1:numconds;
    [rowfname colfname]=size(condfiles_subs{i});
    CHANCHOICES(1:rowfname,j)=cellfun(@strsplit,CHANCHOICES(1:rowfname,j),repmat({' '},rowfname,1),'UniformOutput',0);
    j=j+2;
end
chan_choices=CHANCHOICES;
clearvars -global tmpEEG CONDFILES_SUBS PATHTOFILES CHANCHOICES CURCLICK

end

