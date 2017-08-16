function [STATS]=WinBootCor(STATS,infodisplay,nboot,tr,Ylabel,varargin)

% Calculate winsorized bootstrapped pearson's r on EEG data at each timpoint with as many correlates (Y data) as you wish. Expects you to load Ydata (EEG correlate). This is a MATLAB varaible (.mat) with one column for each external correlate you wish to use (RTs, accuracy, personality measure etc.). The number of columns in the Y variable, must equal the length of Y label (see below). Rows correspond to number of subjects.
% 
% Inputs:
% 
% ***infodisplay***
% 
% A numerical flag (0 or 1). Set to 1 if you would like to see your contrasts, condition names, Xlabels, and Ylabels. ***end***
% 
% ***nboot***
% number of bootstrap samples to take from the paired X and Y data. ***end***
% 
% ***tr***
% percentage to Winsorize from the X and Y data. ***end***
% 
% ***Ylabel***
% A string, or cell array of strings to indicate your Y variable(s) that will be correlated one by one with your EEG waveform (X variable). For example, {'RTs','accuracy', 'perfectionism_scores'}. ***end***
% 
% ***varargin***
% 
% Options are specified in pairs (key -> val)
% 
% [string] ->
% 
%   [numeric] - String and numerical pair indicating the EEG data you would like to use in the correlation. For example,
% 
% Condition 
% 3
% 
% would correlate data from condition 3 with each Y variable. 
% 
% To use difference waves as the X variable, use 'FactorA', or 'FactorB', or 'FactorAB', a number indicating the contrast you wish to use in the correlation, and a label to identify this analysis. For example,
% 
% FactorA
% 2
% myCorrelationAnalysis
% 
% would correlate the difference score corresponding to the 2nd contrast of Factor A to each Y variable and name this analysis “my CorrelationAnalysis”
% 
% For 'bw' designs, if using contrasts (difference scores) in the correlation, you must indicate 'FactorA1', FactorA2, etc, as the strings, as contrasts were taken seperately for each level of factor A. See examples. 
% 
% WinBootCor from commandline::
% 
% [STATS]=WinBootCor('STATS_someanalysis.mat',1,1000, 20,'RT','Condition',3)
%   
% [STATS]=WinBootCor('STATS_someanalysis.mat',1,1000, 20,{'RTs', 'accuracy'},'Condition',1)
%   
% [STATS]=WinBootCor('STATS_someanalysis.mat',1,1000, 20,{'RTs', 'accuracy'},'FactorAB',1,   'interaction')
% 
% [STATS]=WinBootCor('STATS_someanalysis.mat',1,1000, 20,{'RTs', 'accuracy', 'shyness_scores'},'FactorB',2,'Easy_vs_hard_SLEEPY')
%   
% [STATS]=WinBootCor('STATS_someanalysis.mat',1,1000,20,'RTs','FactorA2',2,'Easy_vs_hard_SLEEPY')  % for 'bw' designs
% 
% ***end***
 
% To plot:
% simply use the below and after "corr_results", change to your condition/factor label (type STATS.condnames for a reminder of condition names) followed by Y label:
% figure; plot(STATS.xtimes,STATS.corr_results.SH.RT.rw) %  r values
% figure; plot(STATS.xtimes,STATS.corr_results.SH.RT.p) % p values
% 
% Or if contrasts were used in the correlation:
% figure; plot(STATS.xtimes,STATS.corr_results.interaction.RTs.rw) % r values
% figure; plot(STATS.xtimes,STATS.corr_results.interaction.RTs.p) % p values
% 
% 
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

if strcmp(STATS.design,'ww');
    
    %options = struct('Conditon', [], 'FactorA', [],'FactorB', [], 'FactorAB', []);
    
    if infodisplay
        disp('Condition names'); disp(STATS.condnames)
        %disp('Xlabel'); disp(Xlabel)
        disp('Ylabel'); disp(Ylabel)
        disp('FactorA'); disp(STATS.subject_results.subject_1.factor_A.contrasts)
        disp('FactorB'); disp(STATS.subject_results.subject_1.factor_B.contrasts)
        disp('FactorAB'); disp(STATS.subject_results.subject_1.factor_AxB.contrasts)
    end
    
elseif strcmp(STATS.design,'w');
    
    %options = struct('Conditon', [], 'FactorA', []);
    
    if infodisplay
        disp('Condition names'); disp(STATS.condnames)
        %disp('Xlabel'); disp(Xlabel)
        disp('Ylabel'); disp(Ylabel)
        disp('FactorA'); disp(STATS.subject_results.subject_1.factor_A.contrasts)
    end
    
elseif strcmp(STATS.design,'bw'); % considered not factorial in single-subject cases
    
    if infodisplay
        %disp('j level labels'); disp(STATS.jlabels)
        %disp('k level labels'); disp(STATS.klabels)
        disp('Condition names'); disp(STATS.condnames)
        %disp('Xlabel'); disp(Xlabel)
        disp('Ylabel'); disp(Ylabel)
        disp('contrasts'); disp(STATS.subject_results.Factor_A1.subject_1.factor_A.contrasts)
    end
    
elseif strcmp(STATS.design,'b');
    
    %options = struct('Conditon', [], 'FactorA', []);
    
    if infodisplay
        disp('Condition names'); disp(STATS.condnames)
        %disp('Xlabel'); disp(Xlabel)
        disp('Ylabel'); disp(Ylabel)
        %disp('FactorA'); disp(STATS.inferential_results.factor_A.contrasts)
    end
    
end

% load Ydata, can be as many cols as you have external variables
Yfname=uigetfile('*.mat','Select your EEG correlate (Ydata)');
Ystruct=load(Yfname);
Yfield=fieldnames(Ystruct);
Ydata=Ystruct.(Yfield{1});


switch STATS.design
    case {'w' 'ww'}
        
        
        if any(strcmp(varargin{1},'Condition'));
            
            % assumes that there is equal number of subjects in each cell
            [subfiles_row ~]=size(STATS.subject_bootfiles);
            Xdata=zeros(subfiles_row,STATS.numpnts);
            for i=1:subfiles_row
                tmp=load(STATS.subject_bootfiles{i,varargin{2}});
                Xdata(i,:)=mean(tmp.data,1);
                clear tmp
            end
            
            
        elseif strcmp(varargin{1},'FactorA')
            
            % get subject fields
            subnames = fieldnames(STATS.subject_results);
            
            % assumes that there is equal number of subjects in each cell
            [subfiles_row ~]=size(STATS.subject_bootfiles);
            Xdata=zeros(subfiles_row,STATS.numpnts);
            for i=1:subfiles_row
                Xdata(i,:)=STATS.subject_results.(subnames{i}).factor_A.test_stat(varargin{2},:);
            end
            
            
        elseif strcmp(varargin{1},'FactorB')
            
            % get subject fields
            subnames = fieldnames(STATS.subject_results);
            
            % assumes that there is equal number of subjects in each cell
            [subfiles_row ~]=size(STATS.subject_bootfiles);
            Xdata=zeros(subfiles_row,STATS.numpnts);
            for i=1:subfiles_row
                Xdata(i,:)=STATS.subject_results.(subnames{i}).factor_B.test_stat(varargin{2},:);
            end
            
            
            
        elseif strcmp(varargin{1},'FactorAB')
            
            % get subject fields
            subnames = fieldnames(STATS.subject_results);
            
            % assumes that there is equal number of subjects in each cell
            [subfiles_row ~]=size(STATS.subject_bootfiles);
            Xdata=zeros(subfiles_row,STATS.numpnts);
            for i=1:subfiles_row
                Xdata(i,:)=STATS.subject_results.(subnames{i}).factor_AxB.test_stat(varargin{2},:);
            end
        end
        
    case {'bw' 'b'}
        
        if any(strcmp(varargin{1},'Condition'));
            
            % load Ydata, can be as many cols as you have external variables
            Xfname=uigetfile('*.mat',['Select all subject bootstrapped files in the ' STATS.condnames{1,varargin{2}}, 'condition'], ...
                'MultiSelect' ,'On');
            Xdata=zeros(length(Xfname),STATS.numpnts);
            for i=1:length(Xfname)
                tmp=load(Xfname{i});
                subdata=tmp.data;
                Xdata(i,:)=mean(subdata,1);
                clear subdata tmp
            end
            
        else
            
            try
                % get factor fields
                factlist=fieldnames(STATS.subject_results);
                for i=1:length(factlist)
                    factnames{i}=['FactorA',num2str(i)]; % what the user inputs
                    factnames_instruct{i}=['Factor_A',num2str(i)]; % what the actual structure field is
                end
                
            catch
                disp('cannot look at correlations with difference waves, must run single subject statistics first');
            end
            
            
            
            % check if FactorA1, or FacotrA2 or FactorA3 ect was entered
            if any(strcmp(STATS.design, 'bw') && any(strcmp(varargin{1},factnames)))
                
                % find the factor A level user selected
                factind=strcmp(varargin{1},factnames);
                
                % find out how many subjects for used when stats were calculated
                subnames=fieldnames(STATS.subject_results.(factnames_instruct{factind}));
                Xdata=zeros(length(subnames),STATS.numpnts);
                
                for i=1:length(subnames)
                    Xdata(i,:)=STATS.subject_results.(factnames_instruct{factind}).(subnames{i}).factor_A.test_stat(varargin{2},:);
                end
                
            else
                error('for varargin you must use ''FactorA1'', or ''FactorA2'', ect..., followed by a number indicating the contrast');
                
            end
            
        end
        
        
end

% run correlation with X and Y data
[results]=wincor(nboot,tr,Xdata,Ydata);

% add info about which files were used for the correlations, or which factor and contrasts were used
% this makes subsequent functions easier to use
if any(strcmp(varargin{1},'Condition'));
    if strcmp(STATS.design,'ww') || strcmp(STATS.design,'w')
        results.(varargin{1})=STATS.subject_bootfiles(:,varargin{2});
    elseif strcmp(STATS.design,'bw') || strcmp(STATS.design,'b')
        results.(varargin{1})=Xfname;
    end
else
    results.(varargin{1})=varargin{2};
end

% rename fields to Ylabels
if size(Ydata,2)>1
%if ncols(Ydata)>1
    
    for i=1:size(Ydata,2)
        oldfield=['Y',num2str(i)];
        [results.(Ylabel{i})] = results.(oldfield);
        results = rmfield(results,oldfield);
        
        % add a number to tell me which y coloumn each label is associated
        % with, makes subsequent plotting functions easier to call
        results.(Ylabel{i}).ycol=i;
    end
    
else % if Ydata is only one col
    oldfield='Y1';
    [results.(Ylabel{1})] = results.(oldfield);
    results = rmfield(results,oldfield);
    
    % add a number to tell me which y coloumn each label is associated
    % with, makes subsequent plotting functions easier to call
    results.(Ylabel{1}).ycol=1;
end

% create corr_results structure in main STATS structure and give label
if ~any(strcmp(varargin{1},'Condition')); % if varargin{1} is NOT 'Condition', i.e., Factor selection of some kind
    [STATS.corr_results.(varargin{3})]=results;
else
    [STATS.corr_results.(STATS.condnames{varargin{2}})]=results;
end

disp('******* Saving STATS structure *******')
save(['STATS_',STATS.savestring,'.mat'],'STATS');


end


