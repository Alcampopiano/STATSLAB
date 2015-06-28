function [STATS]=GroupStatistics(STATS,condfiles,alpha,nsamp,varargin)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates group-level statistics for any number of levels, and up to
% two-way designs of any type (between subjects, withing subjects, mixed).
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input arguments:
%     STATS = structre you will be prompted to load this if the argument is left empty. Otherwise give
%             the filename to your STATS stucture in the current directory.
%
%     alpha = arbitrary threshold for significance (.05 or .01 or whatever you like). Will be adjusted with
%             Rom's method for controlling for FWE.
%
%     nsamp = number of resamples to take from the group data (1000, 5000, or what ever you like). Takes time.
%
%
%     varargin - key/val pairs, see Options for details
%
% Options:
%              FEW control method e.g., 'FWE', 'Rom', or 'FWE, 'none',
%              Rom's method is default
%
%              Contrast matrix for Factor A comparisons. For example, 'conA', [1 0 -1 0; 0 1 0 -1]'
%              You must add the transpose operator ('), as the example says.
%              See Wiki (and function usage below) for many more examples. If
%              left out, certain default contrasts will be used, but this is not
%              recommended as you should know what you want to compare.
%
%              Contrast matrix for Factor B comparisons. For example, 'conB', [1 -1 0 0; 0 0 1 -1]'
%              You must add the transpose operator ('), as the example says.
%              See Wiki (and function usage below) for many more examples. If
%              left out, certain default contrasts will be used, but this is not
%              recommended as you should know what you want to compare.
%
%              Contrast matrix for interaction comparisons. For example,
%              'conAB', [1 -1 -1 1]'. You must add the transpose operator ('), as the example says.
%              See Wiki (and function usage below) for many more examples. If
%              left out, certain default contrasts will be used, but this is not
%              recommended as you should know what you want to compare.
%
%
% Examples (these are just examples. For large designs, the number of
% possible pariwise or pooled comparisons increases exponentially):
%
% For a 1-way design with 2 conditions (like a t-test):
% [STATS]=GroupStatistics(STATS,.05,1000, [1 -1]');
%
% For a 2-way design with 4 conditions (2x2):
% [STATS]=GroupStatistics(STATS,.05,1000, 'FWE', 'Rom', [1 0 -1 0; 0 1 0 -1]', [1 0 -1 0; 0 1 0 -1]', [1 0 -1 0; 0 1 0 -1]');
%
% For a 2-way design with 8 conditions (2x4):
% [STATS]=GroupStatistics(STATS,.05,1000, 'FWE', 'Rom', [1 0 0 0 -1 0 0 0; 0 1 0 0 0 -1 0 0]', [1 -1 0 0 0 0 0 0; 0 0 0 0 1 -1 0 0]', [0 0 1 -1 0 0 -1 1]');
%
% For a 1-way design with 3 conditions:
% [STATS]=GroupStatistics(STATS,.05,1000, 'FWE', 'Rom', [1 0 -1; 0 1 -1; 1 -1 0]');
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

% update STATS structure
STATS.alpha=alpha;
STATS.nsamp=nsamp;

% call appropriate stats functions based on input arguments
switch STATS.design
    
    case {'bw','ww','bb'};
        
        
        if any(strcmp({'ersp' 'itc'},STATS.measure));
            
            [inferential_results sample_results condwaves condfiles_subs condwaves_trim] = pbgroup2waytf(STATS, condfiles, STATS.numconds, STATS.timesout, STATS.nboot, ...
                STATS.levels(1), STATS.levels(2), STATS.alpha, STATS.nsamp, STATS.design, STATS.condnames, varargin{:});
            
        elseif any(strcmp({'chanclust' 'gfa'},STATS.measure));
            
            [inferential_results sample_results condwaves condfiles_subs condwaves_trim]=pbgroup2way(condfiles, STATS.numconds, STATS.numpnts, STATS.nboot, ...
                STATS.levels(1), STATS.levels(2), STATS.alpha, STATS.nsamp, STATS.design, STATS.condnames, varargin{:});
                      
        end
        
        
        
        % update STATS structure
        STATS.inferential_results=inferential_results;
        STATS.sample_results=sample_results;
        STATS.condwaves=condwaves;
        STATS.bootfiles=condfiles_subs;
        STATS.condwaves_trim=condwaves_trim;
        
        disp('******* finished calculating statistics *******')
        disp('******* Saving STATS structure *******')
        save(['STATS_',STATS.savestring,'.mat'],'STATS');
        
        
        
    case {'b','w'};
        
        if any(strcmp({'ersp' 'itc'},STATS.measure));
            
            [inferential_results sample_results condwaves condfiles_subs condwaves_trim] = pbgroup1waytf(STATS, condfiles, STATS.numconds, STATS.timesout, STATS.nboot, ...
                STATS.levels(1), STATS.alpha, STATS.nsamp, STATS.design, STATS.condnames, varargin{:});
            
        elseif any(strcmp({'chanclust' 'gfa'},STATS.measure));
            
            [inferential_results sample_results condwaves condfiles_subs condwaves_trim]=pbgroup1way(condfiles, STATS.numconds, STATS.numpnts, STATS.nboot, ...
                STATS.levels(1), STATS.alpha, STATS.nsamp, STATS.design, STATS.condnames, varargin{:});
            
        end
        
        
        % update STATS structure
        STATS.inferential_results=inferential_results;
        STATS.sample_results=sample_results;
        STATS.condwaves=condwaves;
        STATS.bootfiles=condfiles_subs;
        STATS.condwaves_trim=condwaves_trim;
        
        disp('******* finished calculating statistics *******')
        disp('******* Saving STATS structure *******')
        save(['STATS_',STATS.savestring,'.mat'],'STATS');
        
        
    otherwise
        error('design input must be ''bw'', ''ww'', ''bb'', ''w'', or ''b''')
end











end