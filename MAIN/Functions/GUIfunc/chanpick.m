function [okayhit, chan_choices]=chanpick(fnames,numconds)

%size of largest number of subjects gives size of grid
tabsize=max(cell2mat(cellfun(@length,fnames,'un',0)));

evstring='chantable;';

clear -global CHANCHOICES

global CHANCHOICES; 
CHANCHOICES=cell(tabsize,numconds*2);
j=1;
for i=1:numconds;
    [rowfname colfname]=size(fnames{i});
    for q=1:rowfname
        CHANCHOICES{q,j}=fnames{i}{q}; 
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
    'load(fullfile(ParamPath, ParamName),''-mat'');',evstring]} ... %2
    {'Style', 'pushbutton','string','Save channel selection file', ...
    'callback',['[ParamName, ParamPath]=uiputfile(''*.*'',''channel selection file'');', ...
    'global CHANCHOICES;', ...
    'save(fullfile(ParamPath, ParamName),''CHANCHOICES'');']}, ...
    }, ...
    'title', 'STATSLAB -- statslab()',...%, ...
    'eval',evstring ...
    );

% this just converts char arrays to cells so that users don't have to use
% curly braces and single quotes for channel labels, but instead can simply
% delimit with spaces only when needed more than one channel.
j=2;
for i=1:numconds;
    [rowfname colfname]=size(fnames{i});
    CHANCHOICES(1:rowfname,j)=cellfun(@strsplit,CHANCHOICES(1:rowfname,j),repmat({' '},rowfname,1),'UniformOutput',0);
    j=j+2;
end
chan_choices=CHANCHOICES;
clear -global CHANCHOICES

