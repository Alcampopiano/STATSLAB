function [STATS]=GroupStatistics(STATS,condfiles,alpha,varargin)

% Calculates group-level statistics for any number of levels, and up to two-way designs of any type (between subjects, withing subjects, mixed).
%  
% Inputs:
% 
% ***condfiles***
% Leave empty and MATLAB will bring up an interface for you to load the appropriate “bootstrapped” files.
% ***end***
% 
% ***alpha***
% Statistical significance threshold. Will be adjusted if options for controlling FWE are specified.
% ***end***
% 
% ***nsamp***
% Number of resamples to take at the group level (resampling randomly from the subjects themselves). Currently this functionality is disabled so the group level analysis always contains all subjects.
% ***end***
% 
% ***varargin***
% Options are specified in pairs (key -> val)
% 
% FWE ->
% 	
% 	none  - no control for familywise error 
% 	Rom - control FWE using Rom's sequentially rejective method (Wilcox, 2012)
%   Bon - use Bonferroni method to correct for FWE
%  
% conA ->
%  
% 	[numeric] - Contrast matrix for Factor A comparisons. 
% 
% conB -> 
% 
% 	[numeric] - Contrast matrix for Factor B comparisons, if applicable. 
% 
% conAB -> 
% 
% 	[numeric] - Contrast matrix for the interaction, if applicable. 
% 
% For example,
% 
% FWE
% Rom
% conA
% 1 0 -1 0; 0 1 0 -1
% conB
% 1 -1 0 0; 0 0 1 -1
% conAB
% 1 -1 -1 1
% 
% 
% Controls for FWE using Rom's method (Wilcox, 2012; Rom, 1990). For factor A, two comparisons are made: condition 1 versus 3, and  2 versus 4. For factor B, two comparisons are made: condition 1 versus 2, and 3 versus 4. The interaction is also specified ([1 - 2] - [3 - 4]).
% 
% Using GroupStatistics at the commandline:
%  
% For a 1-way design with 2 conditions (like a t-test):
% [STATS]=GroupStatistics(STATS,.05,'FWE', 'none', 'conA', [1 -1]');
%  
% For a 2-way design with 4 conditions (2x2):
% [STATS]=GroupStatistics(STATS,.05,'FWE', 'Rom', 'conA', [1 0 -1 0; 0 1 0 -1]', 'conB', [1 -1 0 0; 0 0 1 -1]', 'conAB', [1 -1 -1 1]');
% 
% For a 1-way design with 3 conditions:
% [STATS]=GroupStatistics(STATS,.05,'FWE', 'Rom', 'conA', [1 0 -1; 0 1 -1; 1 -1 0]');
% 
% ***end***
%
% Copyright (C) <2015>  <Allan Campopiano>
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(STATS)
    [fnamestats]=uigetfile('*.mat','Select the STATS structure for this analysis');
    load(fnamestats);
else
    load(STATS);
end

% set history
[hist_str]=statslab_history(['STATS_', STATS.savestring, '.mat'],condfiles,alpha,varargin);
STATS.history.GroupStatistics=hist_str;

% update STATS structure
STATS.alpha=alpha;
%STATS.nsamp=nsamp;

% call appropriate stats functions based on input arguments
switch STATS.design
    
    case {'bw','ww','bb'};
        
        
        if any(strcmp({'ersp' 'itc'},STATS.measure));
            
            [sample_results condwaves condfiles_subs] = pbgroup2waytf(STATS, condfiles, STATS.numconds, STATS.timesout, STATS.nboot, ...
                STATS.levels(1), STATS.levels(2), STATS.alpha, STATS.design, STATS.condnames, varargin{:});
            
        elseif any(strcmp({'chanclust' 'gfa'},STATS.measure));
            
            [sample_results condwaves condfiles_subs condwaves_trim]=pbgroup2way(STATS,condfiles, STATS.numconds, STATS.numpnts, STATS.nboot, ...
                STATS.levels(1), STATS.levels(2), STATS.alpha, STATS.design, STATS.condnames, varargin{:});
                      
        end
        
        
        
        % update STATS structure
        %STATS.inferential_results=inferential_results;
        STATS.sample_results=sample_results;
        STATS.condwaves=condwaves;
        STATS.bootfiles=condfiles_subs;
        %STATS.condwaves_trim=condwaves_trim;
        
        disp('******* finished calculating statistics *******')
        disp('******* Saving STATS structure *******')
        save(['STATS_',STATS.savestring,'.mat'],'STATS');
        
        
        
    case {'b','w'};
        
        if any(strcmp({'ersp' 'itc'},STATS.measure));
            
            [sample_results condwaves condfiles_subs] = pbgroup1waytf(STATS, condfiles, STATS.numconds, STATS.timesout, STATS.nboot, ...
                STATS.levels(1), STATS.alpha, STATS.design, STATS.condnames, varargin{:});
            
        elseif any(strcmp({'chanclust' 'gfa'},STATS.measure));
            
            [sample_results condwaves condfiles_subs condwaves_trim]=pbgroup1way(STATS,condfiles, STATS.numconds, STATS.numpnts, STATS.nboot, ...
                STATS.levels(1), STATS.alpha, STATS.design, STATS.condnames, varargin{:});
            
        end
        
        
        % update STATS structure
        %STATS.inferential_results=inferential_results;
        STATS.sample_results=sample_results;
        STATS.condwaves=condwaves;
        STATS.bootfiles=condfiles_subs;
        %STATS.condwaves_trim=condwaves_trim;
        
        disp('******* finished calculating statistics *******')
        disp('******* Saving STATS structure *******')
        save(['STATS_',STATS.savestring,'.mat'],'STATS');
        
        
    otherwise
        error('design input must be ''bw'', ''ww'', ''bb'', ''w'', or ''b''')
end











end