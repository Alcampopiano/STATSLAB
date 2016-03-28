function [STATS]=SubjectFigure(STATS,infodisplay,varargin)
% This function plots single-subject data in a variety of ways. Many options exist for non-spectral measures. 
% For ERSP and ITC measures, basic plots are created where the group-level time X frequency for the specified contrast is plotted. 
% In addition, a p-value matrix showing the proportion of subjects with statistical differences at each time and frequency is shown. 
% All figures have clickable features, so click waveforms/frequencies to display more information.
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
% plottype ->
% 	
% 	CI_MOE (default) - plot a figure showing the margin of error (one half of the CI) and display the difference between two conditions as a color bar. 
%                      Black lines show statistical differences.
% 
% 	wave – plot difference as a waveform surrounded by the CI. Red lines show statistical differences.
% 
% 	diff – plot the differences as a color bar. Black lines show statistical differences.
% 
% diffcol ->
% 
% 	[r g b] triplet -  the color of the difference waveform (default is [0 0 0]).
% 
% CIcol -> 
% 
% 	[r g b] triplet - the color of the CI (default is gray [.5 .5 .5]).
% 
% yaxis ->
% 
% 	[numeric] - the yaxis limits (for wave and diff options only). Default scales to absolute max for each subject.
% 
% zaxis -> 
% 
% 	[numeric] - the margin of error limits (for CI_MOE option only). In the range of 0-inf. Default scales to absolute max for each subject.
% caxis ->
% 
% 	 [numeric] - the color range to plot (for CI_MOE and diff options only). Default scales to absolute max for each subject.
% 
% timeplot ->
% 
% 	[numeric] - the range of xaxis values to plot. Default scales to length of epoch.
% 
% topos ->
% 	 no (default) - do not calculate topographies
%      
% 	 channel location file - calculate and allow plotting of topographies (only need to do once).
% 
% FactorA, FactorB, FactorAB, all ->
% 
% 	[numeric] - specify which contrasts you wish to plot for each factor and the interaction
% 
% 	all (default) - plot all contrasts that were specified when statistics were calculated
% 
% 
% For example,
% 
% all
% timeplot
% -200 600
% plottype
% wave
% diffcol
% [0 0 .8]
% yaxis
% -5 5
% topos
% mychanlocsfile.sfp
% 
% This set of options would plot all contrasts in your design. The type of plot would be a standard difference wave in blue for each subject bounded by a CI. The channel locations file would be use to precompute topographies. This allows you to click on the waveforms at any latency for any subject and see the topographies for each condition.
% 
% Using SubjectFigure at the commandline:
% 
% [STATS]=SubjectFigure('STATS_mythesis.mat',1,'all','timeplot',[-200 600]);
% ***end***



if isempty(STATS)
    [fnamestats]=uigetfile('*.mat','Select the STATS structure for this analysis');
    load(fnamestats);
else
    load(STATS);
end

% set history
[hist_str]=statslab_history(['STATS_', STATS.savestring, '.mat'],infodisplay,varargin);
STATS.history.SubjectFigure=hist_str;



if ~any(strcmp({'ersp' 'itc'},STATS.measure));
    [STATS]=pbsubjectfig(STATS,infodisplay,varargin{:});
    
elseif any(strcmp({'ersp' 'itc'},STATS.measure));
    
    [STATS]=pbsubjectfigTF(STATS,infodisplay,varargin{:}); 
end

end









