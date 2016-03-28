function [STATS]=GroupFigure_sample(STATS,infodisplay,varargin)
% 
% This function plots group-level data. 
% For ERSP and ITC measures, basic plots are created where the group-level time X frequency for the specified contrast is plotted. 
% All figures have clickable features, so click waveforms/frequencies to display more information (e.g., topographies, CIs around frequency bands).
% 
% Inputs:
% 
% ***infodisplay***
% Check box to display your condition labels and contrast matrices when plotting figures ***end***
% 
% 
% ***varargin***
% Options are specified in pairs (key -> val)
% 
% timeplot ->
% 
% 	[numeric] - the range of xaxis values to plot. Default scales to length of epoch.
% 
% topos ->
% 	 no (default) - do not calculate topographies
% 
% 	channel location file - calculate and allow plotting of topographies (only need to do once).
% 
% FactorA, FactorB, FactorAB, all (use the "all" option without a key/value pair) ->
% 
% 	[numeric] - specify which contrasts you wish to plot for each factor and the interaction
% 
% 	all (default) - plot all contrasts that were specified when statistics were calculated. 
%                   This option does not need to be paired with another option, just use it on its own 
%                   (i.e., do not pair with FactorA, FactorB, or FactorAB).
% 
% 
% For example,
% 
% FactorA
% 2
% timeplot
% -200 600
% topos
% mychanlocsfile.sfp
% 
% This set of options would plot the second Factor A comparison (otherwise just use the “all” option to plot all contrasts). 
% The type of plot would be a standard difference wave in blue for each subject bounded by a CI. 
% The channel locations file would be use to precompute topographies. 
% This allows you to click on the waveforms at any latency for any subject and see the topographies for each condition.
% 
% Using GroupFigure at the commandline:
% 
% [STATS]=GroupFigure_sample('STATS_mythesis.mat',1,'all','timeplot',[-200 600], 'topos' ,'no');
% ***end***

if isempty(STATS)
    [fnamestats]=uigetfile('*.mat','Select the STATS structure for this analysis');
    load(fnamestats);
else 
    load(STATS);
end

% set history
[hist_str]=statslab_history(['STATS_', STATS.savestring, '.mat'],infodisplay,varargin);
STATS.history.GroupFigure=hist_str;


if ~any(strcmp({'ersp' 'itc'},STATS.measure));
    [STATS]=pbgroupfig_sample(STATS,infodisplay,varargin{:});
    
elseif any(strcmp({'ersp' 'itc'},STATS.measure));
    
    [STATS]=pbgroupfigTF_sample(STATS,infodisplay,varargin{:}); 
end

end










