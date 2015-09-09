function [okayhit, IC_choices]=ICpick(fnames,numconds)

%size of largest number of subjects gives size of grid
tabsize=max(cell2mat(cellfun(@length,fnames,'un',0)));

evstring='ICtable;';

% the ICCHOICES array is an empty array, that gets populated with fname
% elements. The adjacent rows are empty strings '' that get filled with
% components that you want to extract. The trailing rows (if you have
% unequal cell numbers), are empty []). The array gets passed globally to
% input gui in this function, and uitable in ICtable.m
global ICCHOICES; 
ICCHOICES=cell(tabsize,numconds*2);
j=1;
for i=1:numconds;
    [rowfname colfname]=size(fnames{i});
    for q=1:rowfname
        ICCHOICES{q,j}=fnames{i}{q}; 
        ICCHOICES{q,j+1}='';    
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
    ['[ParamName, ParamPath] = uigetfile(''*.mat'',''Select IC cluster file:'',''*.mat'',''multiselect'',''off'');', ...
    'load(fullfile(ParamPath, ParamName),''-mat'');']} ... %2
    {'Style', 'pushbutton','string','Save IC cluster file', ...
    'callback',['[ParamName, ParamPath]=uiputfile(''*.*'',''IC cluster file'');', ...
    'global STATSLAB_PROPERTIES; statslab_propgrid=propgrid2struct(STATSLAB_PROPERTIES);', ...
    'save(fullfile(ParamPath, ParamName),''curdata'');']}, ...
    }, ...
    'title', 'STATSLAB -- statslab()',...%, ...
    'eval',evstring ...
    );

j=2;
for i=1:numconds;
    [rowfname colfname]=size(fnames{i});
    ICCHOICES(1:rowfname,j)=cellfun(@str2num,ICCHOICES(1:rowfname,j),'UniformOutput',0)';
    j=j+2;
end
IC_choices=ICCHOICES;
clearvars -global ICCHOICES






