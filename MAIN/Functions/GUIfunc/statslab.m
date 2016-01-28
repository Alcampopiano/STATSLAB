function [STATS]=statslab(nosplash)

% version
vercheck();

% splash or no splash
if nargin==1;
    splash_fade;
end
    
% statslab_propgrid=''; % so cancel button can be hit
STATS=''; % so cancel button can be hit

PropGridStr=['global STATSLAB_PROPERTIES;' ...
    'properties=propgridbuild();' ...
    'properties = properties.GetHierarchy();' ...
    'STATSLAB_PROPERTIES = PropertyGrid(gcf,' ...
    '''Properties'', properties,' ...
    '''Position'', [.05 .10 .9 .80]);' ...
    ]; 

PropGridStr2=['global STATSLAB_PROPERTIES;' ...
    'properties=propgridbuild(statslab_propgrid);' ...
    'properties = properties.GetHierarchy();' ...
    'STATSLAB_PROPERTIES = PropertyGrid(gcf,' ...
    '''Properties'', properties,' ...
    '''Position'', [.05 .10 .9 .80]);' ...
    ];

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
    ['[ParamName, ParamPath] = uigetfile(''*.mat'',''Select Context configuration file:'',''*.mat'',''multiselect'',''off'');', ...
    'load(fullfile(ParamPath, ParamName),''-mat'');',PropGridStr2]} ... %2
    {'Style', 'pushbutton','string','Save Parameters', ...
    'callback',['[ParamName, ParamPath]=uiputfile(''*.*'',''Context configuration file'');', ...
    'global STATSLAB_PROPERTIES; statslab_propgrid=propgrid2struct(STATSLAB_PROPERTIES);', ...
    'save(fullfile(ParamPath, ParamName),''statslab_propgrid'');']}, ...
    }, ...
    'title', 'STATSLAB -- statslab()',...%, ...
    'eval',PropGridStr ...
    );

if isempty(okayhit);clearvars -global STATSLAB_PROPERTIES;return;end

global STATSLAB_PROPERTIES;
statslab_propgrid=propgrid2struct(STATSLAB_PROPERTIES);
clearvars -global STATSLAB_PROPERTIES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pass GUI inputs into stats functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if statslab_propgrid.logical_extract==1;
    
    % initiate partial function string
    funcstr=['ExtractData(statslab_propgrid.condnames,statslab_propgrid.condfiles_extract,'...
        'statslab_propgrid.levels,statslab_propgrid.design,statslab_propgrid.savestring,'];
    
    [funcstr]=varargin_spilt(funcstr,statslab_propgrid.varargin_extract);
    
    STATS=eval(funcstr);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if statslab_propgrid.logical_resample==1;
    
    % initiate partial function string
    funcstr='ResampleData([''STATS_'', statslab_propgrid.savestring, ''.mat''],statslab_propgrid.condfiles_resample,statslab_propgrid.nboot,';
    
    [funcstr]=varargin_spilt(funcstr,statslab_propgrid.varargin_resample);
    
    STATS=eval(funcstr);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if statslab_propgrid.logical_group==1;
    
    % initiate partial function string
    funcstr=['GroupStatistics([''STATS_'', statslab_propgrid.savestring, ''.mat''],statslab_propgrid.condfiles_group,'...
        'statslab_propgrid.alpha_group,statslab_propgrid.nsamp,'];
    
    [funcstr]=varargin_spilt(funcstr,statslab_propgrid.varargin_group);
    
    STATS=eval(funcstr);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if statslab_propgrid.logical_subject==1;
    
    % initiate partial function string
    funcstr=['SubjectStatistics([''STATS_'', statslab_propgrid.savestring, ''.mat''],statslab_propgrid.condfiles_subject,'...
        'statslab_propgrid.alpha_group,'];
    
    [funcstr]=varargin_spilt(funcstr,statslab_propgrid.varargin_subject);
    
    STATS=eval(funcstr);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if statslab_propgrid.logical_groupfig==1;
    
    % initiate partial function string
    funcstr='GroupFigure_sample([''STATS_'', statslab_propgrid.savestring, ''.mat''],statslab_propgrid.infodisplay_group,';
    
    [funcstr]=varargin_spilt(funcstr,statslab_propgrid.varargin_groupfig);
    
    STATS=eval(funcstr);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if statslab_propgrid.logical_subjectfig==1;
    
    % initiate partial function string
    funcstr='SubjectFigure([''STATS_'', statslab_propgrid.savestring, ''.mat''],statslab_propgrid.infodisplay_sub,';
    
    [funcstr]=varargin_spilt(funcstr,statslab_propgrid.varargin_subjectfig);
    
    STATS=eval(funcstr);
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if statslab_propgrid.logical_cor==1;
    
    % initiate partial function string
    funcstr=['WinBootCor([''STATS_'', statslab_propgrid.savestring, ''.mat''],statslab_propgrid.infodisplay_cor,'...
        'statslab_propgrid.nboot_cor, statslab_propgrid.tr, statslab_propgrid.Ylabel_cor,'];
    
    [funcstr]=varargin_spilt(funcstr,statslab_propgrid.varargin_cor);
    
    STATS=eval(funcstr);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if statslab_propgrid.logical_rw==1;
    

    % initiate partial function string
    funcstr=['WinBootCorFigure([''STATS_'', statslab_propgrid.savestring, ''.mat''],statslab_propgrid.condition_label_rw,'...
        'statslab_propgrid.Ylable_rw);'];
    
    STATS=eval(funcstr);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if statslab_propgrid.logical_slope==1;
    
    % initiate partial function string
    funcstr=['SlopeCI([''STATS_'', statslab_propgrid.savestring, ''.mat''],statslab_propgrid.infodisplay_slope, statslab_propgrid.condition_label_slope,'...
        'statslab_propgrid.Ylable_slope, statslab_propgrid.msplot_slope, statslab_propgrid.CI_color_slope, statslab_propgrid.CI_limit_slope);'];
    
    STATS=eval(funcstr);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if statslab_propgrid.logical_lowess==1;
    
    % initiate partial function string
    funcstr=['LowessCI([''STATS_'', statslab_propgrid.savestring, ''.mat''],statslab_propgrid.infodisplay_lowess, statslab_propgrid.condition_label_lowess,'...
        'statslab_propgrid.Ylable_lowess, statslab_propgrid.msplot_lowess, statslab_propgrid.nboot_lowess, statslab_propgrid.span, statslab_propgrid.nbins);'];
    
    STATS=eval(funcstr);
    
end

















