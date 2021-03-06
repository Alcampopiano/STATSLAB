function [STATS]=LowessCI(STATS,infodisplay,Xlabel,Ylabel,msplot,nboot,span,nbins)

% Plot visually weighted lowess confidence intervals.
% 
% Inputs:
%  
% ***infodisplay***
%  A numerical flag (0 or 1). Set to 1 if you would like to see your contrasts, condition names, Xlabels, and Ylabels. ***end***
%  
% ***Xlabel***
% String that indicates the X variable (EEG data). Must be the same Xlabel used when running WinBootCor.m ***end***
%  
% ***Ylabel*** 
% A string indicating your Y variable(s) (correlates) Must be the same Ylabel used when running WinBootCor.m ***end***
%  
% ***msplot***
% a number indicating the ms you wish to plot ***end***
%  
% ***nboot*** 
% number of bootstrap lowess samples ***end***
% 
%  ***span***
% number in range 0 to 1, indicates size of sliding window (percent of data to use). Try .2 ***end***
%  
% ***nbins***
% number of vertical bins in which to calculate kernal density. Try a small number to begin with (50), as large numbers take more computing time and memory ***end***
%  
% LowessCI from commandline: 
%  
% LowessCI('STATS_B_analysis.mat',1,'awake','RT',100, 200, .2, 50)


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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(STATS)
    [fnamestats]=uigetfile('*.mat','Select the STATS structure for this analysis');
    load(fnamestats);
else
    load(STATS);
end

% ms to TF
x_min=min(STATS.xtimes)/1000;
x_max=max(STATS.xtimes)/1000;

% this eqauls the TF to plot at NOT the ms
msplot=round((msplot/1000-x_min)/(x_max-x_min) * (STATS.numpnts-1))+1;

if infodisplay
    disp('Winsorized/bootstrapped r value'); disp(STATS.corr_results.(Xlabel).(Ylabel).rw(msplot))
    disp('CI around rw'); disp(STATS.corr_results.(Xlabel).(Ylabel).CI(1:2,msplot));
    disp('p  value'); disp(STATS.corr_results.(Xlabel).(Ylabel).p(msplot))
    disp('Time (ms)'); disp(msplot)
end

if strcmp(STATS.design,'ww');
    
    if infodisplay
        disp('Condition names'); disp(STATS.condnames)
        disp('Xlabel'); disp(Xlabel)
        disp('Ylabel'); disp(Ylabel)
        disp('FactorA'); disp(STATS.subject_results.subject_1.factor_A.contrasts)
        disp('FactorB'); disp(STATS.subject_results.subject_1.factor_B.contrasts)
        disp('FactorAB'); disp(STATS.subject_results.subject_1.factor_AxB.contrasts)
        %disp('Winsorized/bootstrapped r value'); disp(STATS.corr_results.(Xlabel).(Ylabel).rw)
        %disp('CI around rw'); disp(STATS.corr_results.(Xlabel).(Ylabel).CI)
        %disp('p  value'); disp(STATS.corr_results.(Xlabel).(Ylabel).p)
    end
    
elseif strcmp(STATS.design,'w');
    
    %options = struct('Conditon', [], 'FactorA', []);
    
    if infodisplay
        disp('Condition names'); disp(STATS.condnames)
        disp('Xlabel'); disp(Xlabel)
        disp('Ylabel'); disp(Ylabel)
        disp('FactorA'); disp(STATS.subject_results.subject_1.factor_A.contrasts)
    end
    
elseif strcmp(STATS.design,'bw'); % considered not factorial in single-subject cases
    
    %options = struct('Conditon', []);
    
    % add fields that say FactorA1, FactorA2, etc...
    %factorAlevels=fieldnames(STATS.subject_results);
    %for i=1:length(factorAlevels);
    %    tempfieldnames=['FactorA', num2str(i)];
    %    options.(tempfieldnames)=[];
    %end
    
    if infodisplay
        %disp('j level labels'); disp(STATS.jlabels)
        %disp('k level labels'); disp(STATS.klabels)
        disp('Condition names'); disp(STATS.condnames)
        disp('Xlabel'); disp(Xlabel)
        disp('Ylabel'); disp(Ylabel)
        disp('contrasts'); disp(STATS.subject_results.Factor_A1.subject_1.factor_A.contrasts)
    end
    
elseif strcmp(STATS.design,'b');
    
    %options = struct('Conditon', [], 'FactorA', []);
    
    if infodisplay
        disp('Condition names'); disp(STATS.condnames)
        disp('Xlabel'); disp(Xlabel)
        disp('Ylabel'); disp(Ylabel)
        %disp('FactorA'); disp(STATS.inferential_results.factor_A.contrasts)
    end
    
end

% else check for name/val agreement and set defaults to empty
%nargs = length(varargin);
%if round(nargs/2)~=nargs/2
%    error('need propertyName/propertyValue pairs for optional inputs')
%end

% load Ydata, can be as many cols as you have external variables
Yfname=uigetfile('*.mat','Select your EEG correlate (Ydata)');
Ystruct=load(Yfname);
Yfield=fieldnames(Ystruct);
YdataAll=Ystruct.(Yfield{1});
Ydata=YdataAll(:,STATS.corr_results.(Xlabel).(Ylabel).ycol);

Y_fields=fieldnames(STATS.corr_results.(Xlabel));

if any(strcmp(Y_fields,'Condition'))
    
    % assumes that there is equal number of subjects in each cell
    [subfiles_row ~]=size(STATS.corr_results.(Xlabel).Condition);
    Xdata=zeros(subfiles_row,STATS.numpnts);
    for i=1:subfiles_row
        tmp=load(STATS.corr_results.(Xlabel).Condition{i});
        Xdata(i,:)=mean(tmp.data,1);
        clear tmp
    end
    
elseif any(strcmp(Y_fields,'FactorA'))
    % extract factorA waves
    % get subject fields
    subnames = fieldnames(STATS.subject_results);
    rowcon=STATS.corr_results.(Xlabel).FactorA;
    
    % assumes that there is equal number of subjects in each cell
    [subfiles_row ~]=size(subnames);
    Xdata=zeros(subfiles_row,STATS.numpnts);
    for i=1:subfiles_row
        Xdata(i,:)=STATS.subject_results.(subnames{i}).factor_A.test_stat(rowcon,:);
    end
    
    
elseif any(strcmp(Y_fields,'FactorB'))
    % extract factorB waves
    % get subject fields
    subnames = fieldnames(STATS.subject_results);
    rowcon=STATS.corr_results.(Xlabel).FactorB;
    
    % assumes that there is equal number of subjects in each cell
    [subfiles_row ~]=size(subnames);
    Xdata=zeros(subfiles_row,STATS.numpnts);
    for i=1:subfiles_row
        Xdata(i,:)=STATS.subject_results.(subnames{i}).factor_B.test_stat(rowcon,:);
    end
    
    
elseif any(strcmp(Y_fields,'FactorAB'))
    % extract factorAB waves
    % get subject fields
    subnames = fieldnames(STATS.subject_results);
    rowcon=STATS.corr_results.(Xlabel).FactorAB;
    
    % assumes that there is equal number of subjects in each cell
    [subfiles_row ~]=size(subnames);
    Xdata=zeros(subfiles_row,STATS.numpnts);
    for i=1:subfiles_row
        Xdata(i,:)=STATS.subject_results.(subnames{i}).factor_AxB.test_stat(rowcon,:);
    end
    
    
elseif strcmp(STATS.design,'bw')
    
    % levels A fields
    Alevels=fieldnames(STATS.subject_results);
    for i=1:length(Alevels)
        factnames=['FactorA',num2str(i)]; % what the user inputs
        
        if any(strcmp(factnames,Y_fields))
            % find out how many subjects for used when stats were calculated
            subnames=fieldnames(STATS.subject_results.(Alevels{i}));
            Xdata=zeros(length(subnames),STATS.numpnts);
            rowcon=STATS.corr_results.(Xlabel).(factnames);
            
            for j=1:length(subnames)
                Xdata(j,:)=STATS.subject_results.(Alevels{i}).(subnames{j}).factor_A.test_stat(rowcon,:);
            end
        end
    end
    
end

X=Xdata(:,msplot);
Y=Ydata;

% lowess kernal density
lowess_kd(X,Y,nboot,span,nbins)
h(1)=ylabel(Ylabel, 'fontsize', 18);
h(2)=xlabel(Xlabel, 'fontsize', 18);
set(h([1 2]),'Interpreter', 'none');

end









