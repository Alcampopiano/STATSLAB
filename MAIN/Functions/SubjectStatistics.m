function [STATS]=SubjectStatistics(STATS,condfiles,alpha,varargin)

% Calculates subject-level statistics for any number of levels, and up to two-way designs.
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
% ***varargin***
% Options are specified in pairs (key -> val)
% 
% FWE ->
% 	
% 	none  - no control for familywise error 
% 	Rom - control FWE using Rom's sequentially rejective method (Wilcox, 2012)
%  
% ConA ->
%  
% 	[numeric] - Contrast matrix for Factor A comparisons. 
% 
% ConB -> 
% 
% 	[numeric] - Contrast matrix for Factor B comparisons, if applicable. 
% 
% ConAB -> 
% 
% 	[numeric] - Contrast matrix for the interaction, if applicable. 
% 
% jlables -> 
% 
% 	For between-within designs only. Name of the general within-subjects conditions. 
% 
% klabels ->
% 
% 	For between-within designs only. Name of the general between-subjects conditions.
% 
% For example,
% 
% FWE
% Rom
% ConA
% 1 0 -1 0; 0 1 0 -1
% ConB
% 1 -1 0 0; 0 0 1 -1
% ConAB
% 1 -1 -1 1
% 
% 
% Controls for FWE using Rom's method (Wilcox, 2012; Rom, 1990). For factor A, two comparisons are made: condition 1 versus 3, and  2 versus 4. For factor B, two comparisons are made: condition 1 versus 2, and 3 versus 4. The interaction is also specified ([1 - 2] - [3 – 4]).
% 
% FWE
% none
% ConA
% 1 -1 0
% klabels
% face, house, butterfly
% jlabels
% male, female
% 
% No control for FWE. In this example we have a 2x3 mixed design comparing male and females in a face processing task using three visual stimuli (face, house, butterfly). Because this is a mixed design, we can only compare the within-subjects conditions (face, house butterfly) when running single-subject statistics. 
% 
% Using SubjectStatistics at the commandline:
%  
% For a 1-way design with 2 conditions (like a t-test):
% [STATS]=SubjectStatistics(STATS,.05,'FWE', 'none', 'conA', [1 -1]');
%  
% For a 2-way design with 4 conditions (2x2):
% [STATS]=SubjectStatistics(STATS,.05,'FWE', 'Rom', 'conA', [1 0 -1 0; 0 1 0 -1]', 'conB', [1 -1 0 0; 0 0 1 -1]', 'conAB', [1 -1 -1 1]');
% 
% For a 1-way design with 3 conditions:
% [STATS]=SubjectStatistics(STATS,.05,'FWE', 'Rom', 'conA', [1 0 -1; 0 1 -1; 1 -1 0]');
% 
% For a 2-way between-within design with 3 conditions for each level of Factor A:
% [STATS]=SubjectStatistics(STATS,.05,'FWE', 'Rom', 'conA', [1 0 -1; 0 1 -1; 1 -1 0]', 'jlabels', {'male', 'female'}, 'klabels', {'face', 'house', 'butterfly'});
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
STATS.history.SubjectStatistics=hist_str;

% update STATS structure
STATS.alpha=alpha;

% call appropriate stats functions based on input arguments
switch STATS.design
    
    case 'ww';
        
        
        if any(strcmp({'ersp' 'itc'},STATS.measure));
            
            [condfiles results]=pbsubject2waytf(STATS,condfiles,STATS.numconds, STATS.nboot, ...
                STATS.levels(1), STATS.levels(2), STATS.alpha, STATS.condnames, varargin{:});
            
        elseif any(strcmp({'chanclust' 'gfa'},STATS.measure));
            
            [condfiles results]=pbsubject2way(STATS,condfiles,STATS.numconds, STATS.numpnts, STATS.nboot, ...
                STATS.levels(1), STATS.levels(2), STATS.alpha, STATS.condnames, varargin{:});
            
        end
        
        
        % update STATS structure
        STATS.subject_results=results;
        STATS.subject_bootfiles=condfiles;
        
        disp('******* finished calculating statistics *******')
        disp('******* Saving STATS structure *******')
        save(['STATS_',STATS.savestring,'.mat'],'STATS');
        
        
        
    case 'w';
        
        if any(strcmp({'ersp' 'itc'},STATS.measure));
            
            [condfiles results]=pbsubject1waytf(STATS,condfiles,STATS.numconds, STATS.nboot, ...
                STATS.levels(1), STATS.alpha, STATS.condnames, varargin{:});
            
        elseif any(strcmp({'chanclust' 'gfa'},STATS.measure));
            
            [condfiles results]=pbsubject1way(STATS,condfiles,STATS.numconds, STATS.numpnts, STATS.nboot, ...
                STATS.levels(1), STATS.alpha, STATS.condnames, varargin{:});
            
        end
        
        % update STATS structure
        STATS.subject_results=results;
        STATS.subject_bootfiles=condfiles;
        
        disp('******* finished calculating statistics *******')
        disp('******* Saving STATS structure *******')
        save(['STATS_',STATS.savestring,'.mat'],'STATS');
        
    case 'bw';
        
    
        % set default plot options
        options.jlabels={};
        options.klabels={};
        options.conA=[];
        options.FWE='Rom';

        % get field names
        optionnames = fieldnames(options);
            
        for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
            inpName = pair{1};
            
            if any(strcmp(inpName,optionnames))
                
                % overwrite default options
                options.(inpName) = pair{2};
            else
                error('%s is not a recognized parameter name',inpName)
            end
        end
        
        % move labels to STATS structure
        STATS.jlabels=options.jlabels;
        STATS.klabels=options.klabels;
        
        % each level of factor A is run as a one way, within
        % subjects design
        
        % get field names for each level of factor A
        for i=1:STATS.levels(1);
            field_name{i,1}=['Factor_A', num2str(i)];
        end
        
        % loop and gather the filenames for each level of Factor A separately
        %%%%
        if isempty(condfiles);
            
            % stop being an empty character
            clear condfiles;
            
            k=1;
            for j=1:STATS.levels(1)
                
                % load all file names subs X conditions
                for i=1:STATS.levels(2)
                    tempfname=uigetfile('*.mat',['Select all bootstrapped files in the ', STATS.condnames{k}, ' condition'], 'MultiSelect','on');
                    
                    if ~iscell(tempfname) && ~isnumeric(tempfname);
                        tempfname={tempfname};
                        condfiles(:,i)=tempfname;
                        
                    elseif isnumeric(tempfname) % if cancel was hit
                        error('no files were chosen');
                    else
                        condfiles(:,i)=tempfname;
                    end
                    k=k+1;
                end
                
                % save files, file names, and update STATS structure
                save(['condfiles_SubjectStatistics_',STATS.measure,'_',STATS.savestring, '_j',num2str(j),'.mat'],'condfiles');
                STATS.jlvlfnames{j}=['condfiles_SubjectStatistics_',STATS.measure,'_',STATS.savestring, '_j',num2str(j),'.mat'];
                STATS.subject_bootfiles{j}=condfiles;
                clear condfiles
            end
            
                % save filenames
                condfiles=STATS.jlvlfnames;
                save(['condfiles_SubjectStatistics_',STATS.measure,'_',STATS.savestring,'.mat'],'condfiles');
                % redundant
%                 condfiles=STATS.subject_bootfiles;

        else
            
            % load a file name that was given that contains the filenames X condition cell array
            condfiles_data=load(condfiles);
            condfields=fieldnames(condfiles_data);
            condfiles=condfiles_data.(condfields{1});
            
            % update STATS structure
            STATS.jlvlfnames=condfiles;
            
            % just to populate subject_bootfiles in case it needs to be updated
            for j=1:STATS.levels(1);
                tmp=load(STATS.jlvlfnames{j});
                tmpfields=fieldnames(tmp);
                tmpdat=tmp.(tmpfields{1});
                STATS.subject_bootfiles{j}=tmpdat;
            end
        end
        
        for i=1:STATS.levels(1)
            
            if any(strcmp({'ersp' 'itc'},STATS.measure));
                
                [jnk, results]=pbsubject1waytf(STATS,STATS.jlvlfnames{i},STATS.levels(2),STATS.nboot, ...
                    STATS.levels(1), STATS.alpha, STATS.condnames, varargin{:});
                
            elseif any(strcmp({'chanclust' 'gfa'},STATS.measure));
                
                % do something in a loop for factor A
                [jnk, results]=pbsubject1way(STATS,STATS.jlvlfnames{i},STATS.levels(2), STATS.numpnts, STATS.nboot, ...
                    STATS.levels(1), STATS.alpha, STATS.condnames, varargin{:});
            end
            % update STATS structure
            STATS.subject_results.(field_name{i})=results;
            disp('******* Working on next level of Factor A *******')
            
        end
        
        disp('******* finished calculating statistics *******')
        disp('******* Saving STATS structure *******')
        save(['STATS_',STATS.savestring,'.mat'],'STATS');
        
    otherwise
        error('design input must be ''bw'', ''ww'', or ''w'', in order to do single-subject statistics');
end











end