function [STATS]=SubjectStatistics(STATS,alpha,varargin)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates subject-level statistics. Keep in mind that by definition one
% can only calculate single-subject stats for within subjects designs.
% However, for a between-within (mixed) design you can use this function to
% calculate single-subject stats, this will be done along the k levels of factor B (the "within" levels).
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Input arguments:
%     STATS = structre you will be prompted to load this if the argument is left empty. Otherwise give
%             the filename to your STATS stucture in the current directory.
% 
%     alpha = arbitrary threshold for significance (.05 or .01 or whatever you like). Will be adjusted with
%             Rom's method for controlling for FWE.
% 
%     varargin - key/val pairs, see Options for details
% 
% Options:
%              Contrast matrix for Factor A comparisons. For example, [1 0 -1 0; 0 1 0 -1]'
%              You must add the transpose operator ('), as the example says.
%              See Wiki (and function usage below) for many more examples. If
%              left out, certain default contrasts will be used, but this is not
%              recommended as you should know what you want to compare.
% 
%              Contrast matrix for Factor B comparisons. For example, [1 -1 0 0; 0 0 1 -1]'
%              You must add the transpose operator ('), as the example says.
%              See Wiki (and function usage below) for many more examples. If
%              left out, certain default contrasts will be used, but this is not
%              recommended as you should know what you want to compare.
% 
%              Contrast matrix for interaction comparisons. For example, 
%              [1 -1 -1 1]'. You must add the transpose operator ('), as the example says.
%              See Wiki (and function usage below) for many more examples. If
%              left out, certain default contrasts will be used, but this is not
%              recommended as you should know what you want to compare.
% 
% 
% Examples (these are just examples. For large designs, the number of
% possible pariwise or pooled comparisons increases exponentially):
% 
% For a 1-way design with 2 conditions (like a t-test):
% [STATS]=SubjectStatistics(STATS,.05,[1 -1]');
% 
% For a 2-way, within-within design with 4 conditions (2x2):
% [STATS]=SubjectStatistics(STATS,.05,[1 0 -1 0; 0 1 0 -1]', [1 0 -1 0; 0 1 0 -1]', [1 0 -1 0; 0 1 0 -1]');
% 
% For a 2-way within-within design with 8 conditions (2x4):
% [STATS]=SubjectStatistics(STATS,.05,[1 0 0 0 -1 0 0 0; 0 1 0 0 0 -1 0 0]', [1 -1 0 0 0 0 0 0; 0 0 0 0 1 -1 0 0]', [0 0 1 -1 0 0 -1 1]');
% 
% For a 1-way within design with 3 conditions:
% [STATS]=SubjectStatistics(STATS,.05,[1 0 -1; 0 1 -1; 1 -1 0]');
% 
% For a 2-way between-within design with 4 conditions (2x2):
% [STATS]=SubjectStatistics(STATS,.05,[1 -1]',{'A','S'},{'E','H'})
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

% call appropriate stats functions based on input arguments
switch STATS.design
    
    case 'ww';
        [condfiles results]=pbsubject2way(STATS.numconds, STATS.numpnts, STATS.nboot, ...
            STATS.levels(1), STATS.levels(2), STATS.alpha, STATS.condnames, varargin{:});
        
        % update STATS structure
        STATS.subject_results=results;
        STATS.subject_bootfiles=condfiles;
        
        disp('******* finished calculating statistics *******')
        disp('******* Saving STATS structure *******')
        save(['STATS_',STATS.savestring,'.mat'],'STATS');
        
        
        
    case 'w';
        [condfiles results]=pbsubject1way('',STATS.numconds, STATS.numpnts, STATS.nboot, ...
            STATS.levels(1), STATS.alpha, STATS.condnames, varargin{:});
        
        % update STATS structure
        STATS.subject_results=results;
        STATS.subject_bootfiles=condfiles;
        
        disp('******* finished calculating statistics *******')
        disp('******* Saving STATS structure *******')
        save(['STATS_',STATS.savestring,'.mat'],'STATS');
        
    case 'bw';
        
        %         % handle varargin errors
        %         if length(varargin)~=3;
        %             error([' You must specify one contrast matrix for bw design when doing single-subject stats' ...
        %                 ' Also, you also need to specify labels for factor A and B seperatley.' ...
        %                 ' See documentation and function header for examples']);
        %
        %
        %         elseif length(varargin)==3
        %             if length(varargin{2})*length(varargin{3})~=STATS.levels(1)*STATS.levels(2)
        %                 error([' The number of labels you''ve entered is incorrect for a', num2str(STATS.levels(1)), 'X', num2str(STATS.levels(2)), 'design', ...
        %                     ' See documentation and function header for examples']);
        %
        %             else
        %                 % parse trailing inputs from varargin into labels
        %                 STATS.jlabels=varargin{2};
        %                 STATS.klabels=varargin{3};
        %
        %             end
        %         end
        
        % set default plot options
        options.jlabels='no_jlabels';
        options.klabels='no_klabels';
        
        % get field names
        optionnames = fieldnames(options);
        
        % % check to see which optional args were used and deal with accordingly
        % if isempty(varargin);
        %     warning('MATLAB:stats',['Using default contrasts matrix. You must specify one if you want a custom contrast. ' ...
        %         ' e.g., [1 -1 0; 1 0 -1]'''])
        % else
        %     % overwrite options stucture with varargin inputs if there are any
        %
        %     if ~isempty(varargin{1})
        %         options.(optionnames{1})=varargin{1};
        %     end
        % end
        
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

        % basically, each level of factor A is run as a one way, within
        % subjects design
        
        % get field names for each level of factor A
        for i=1:STATS.levels(1);
            field_name{i,1}=['Factor_A', num2str(i)];
        end

        % loop and gather the filenames for each level of Factor A separately
        k=1;
        for j=1:STATS.levels(1)
            
            % load all file names subs X conditions
            for i=1:STATS.levels(2)
                tempfname=uigetfile('*.mat',['Select all bootstrapped files in the ', STATS.condnames{k}, ' condition'], 'MultiSelect','on');
                if ~iscell(tempfname);
                    tempfname={tempfname};
                    condfiles(:,i)=tempfname;
                else
                    condfiles(:,i)=tempfname;
                end
                k=k+1;
            end
            
            % update STATS structure
            STATS.subject_bootfiles{j}=condfiles;
            clear condfiles
        end
        
        for i=1:STATS.levels(1)
            
            % do something in a loop for factor A
            [~, results]=pbsubject1way(STATS.subject_bootfiles{i},STATS.levels(2), STATS.numpnts, STATS.nboot, ...
                STATS.levels(1), STATS.alpha, STATS.condnames, varargin{:});
            
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