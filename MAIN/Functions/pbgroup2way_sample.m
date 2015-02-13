function [results condwaves] = pbgroup2way_sample(numconds, numpnts, nboot, jlvls, klvls, alpha, condfiles_subs, varargin)
%{

contrasts using con2way default to pooling across one factor. However, you
may want to compare pairwise differences that are not pooled across another factor.
This difference in approach can have a practical difference in the results
and interpretation (See Wilcox, 2013).

This function should allow one to either (1) compare all pairwise
differences by pooling, (2) compare all pairwise differences without
pooling, and (3) specify a priori contrasts, which should not include ALL
pairwise differecnes, but fewer contrasts of course.

%}

% Set default contrast coefficients for 2-way
% create contrasts for 2way ANOVA (used for multi-comparisons)
[conA conB conAB] = con2way(jlvls, klvls);

% put defaults into a structure;
options=struct('conA',conA,'conB',conB,'conAB',conAB);

% get field names
optionnames = fieldnames(options);

% check to see which optional args were used and deal with accordingly
if isempty(varargin);
    warning('MATLAB:stats',['using default contrasts matrices for factor A, B, and the interaction ' ...
        'you must specify each one as separate optional input arguments if you want custom contrasts. ' ...
        'You can leave certain contrasts empty if you want the default comparisons, e.g., [], [1 1 -1 -1]'' ,[]'])
else
    % overwrite options stucture with varargin inputs if there are any
    for i=1:3;
        if ~isempty(varargin{i})
            options.(optionnames{i})=varargin{i};
        end
    end
end

% extract from options structure
conA=options.(optionnames{1});
conB=options.(optionnames{2});
conAB=options.(optionnames{3});

% used to create proper sizes in results structure
[~, conAcol]=size(conA);
[~, conBcol]=size(conB);
[~, conABcol]=size(conAB);

% load data
if isempty(condfiles_subs);
    for i=1:numconds
        tempfname=uigetfile('*.mat',['Select all subject condition ', num2str(i), ' files'], 'MultiSelect','on');       
        condfiles_subs{1,i}(:,1)=tempfname;  
    end   
    [datacell] = grandaverage(nboot,numpnts,condfiles_subs);  
else
    [datacell] = grandaverage(nboot,numpnts,condfiles_subs);
end

% get condition waveforms for plotting purposes
for i=1:length(datacell);
    condwaves(i,:)=mean(datacell{i},1);
end

%preallocate sizes
[rowconds colconds]=size(condfiles_subs);

% build results structure
results=struct('factor_A',{[]},'factor_B',{[]},'factor_AxB',{[]});
results.factor_A=struct('contrasts',{conA},'pval',{zeros(conAcol,numpnts)},'alpha',{zeros(conAcol,numpnts)},'test_stat',{zeros(conAcol,numpnts)},'CI',{cell(conAcol,1)});
%results.factor_A.z_scores=struct('z_avg',{zeros(conAcol,numpnts)},'CI_z',{cell(conAcol,1)});

for i=1:conAcol;
    results.factor_A.CI{i,1}=zeros(2,numpnts);
    %results.factor_A.z_scores.CI_z{i,1}=zeros(2,numpnts);
end

results.factor_B=struct('contrasts',{conB},'pval',{zeros(conBcol,numpnts)},'alpha',{zeros(conBcol,numpnts)},'test_stat',{zeros(conBcol,numpnts)},'CI',{cell(conBcol,1)});
%results.factor_B.z_scores=struct('z_avg',{zeros(conBcol,numpnts)},'CI_z',{cell(conBcol,1)});

for i=1:conBcol;
    results.factor_B.CI{i,1}=zeros(2,numpnts);
    %results.factor_B.z_scores.CI_z{i,1}=zeros(2,numpnts);
    
end

results.factor_AxB=struct('contrasts',{conAB},'pval',{zeros(conABcol,numpnts)},'alpha',{zeros(conABcol,numpnts)},'test_stat',{zeros(conABcol,numpnts)},'CI',{cell(conABcol,1)});
%results.factor_AxB.z_scores=struct('z_avg',{zeros(conABcol,numpnts)},'CI_z',{cell(conABcol,1)});

for i=1:conABcol;
    results.factor_AxB.CI{i,1}=zeros(2,numpnts);
    %results.factor_AxB.z_scores.CI_z{i,1}=zeros(2,numpnts);
end
%end

% load and arrange data
h2 = waitbar(0,'1','Name','analysis using all subjects','Position',[1100 486 550 40]);
childh2 = get(h2, 'Children');
set(childh2, 'Position',[5 10 538 15]);

%arrange the data for the calculations
[rowcell ~]=size(datacell{1,1});

% loop for stats at each timepoint
for timecurrent=1:numpnts;
    
    % reset data to zeros after every calculation at each timepoint
    data=zeros(rowcell,numconds);
    
    % arrange data into a matrix with subs (or single subject boot samples) X conditions
    for condcurrent=1:colconds;
        data(:,condcurrent)=datacell{1,condcurrent}(:,timecurrent);
    end

    % factor A
    con=conA;
    [psihat_stat pvalgen pcrit conflow confup]=pbstats(data, con, nboot, alpha);
    
    % passing results into results structure
    results.factor_A.pval(:,timecurrent)=pvalgen;
    results.factor_A.alpha(:,timecurrent)=pcrit;
    results.factor_A.test_stat(:,timecurrent)=psihat_stat;
    
    for i=1:conAcol;
        results.factor_A.CI{i,1}(1,timecurrent)=conflow(i);
        results.factor_A.CI{i,1}(2,timecurrent)=confup(i);
    end
    
    % factor B
    con=conB;
    [psihat_stat pvalgen pcrit conflow confup]=pbstats(data, con, nboot, alpha);
    
    % passing results into results structure
    results.factor_B.pval(:,timecurrent)=pvalgen;
    results.factor_B.alpha(:,timecurrent)=pcrit;
    results.factor_B.test_stat(:,timecurrent)=psihat_stat;
    
    for i=1:conBcol;
        results.factor_B.CI{i,1}(1,timecurrent)=conflow(i);
        results.factor_B.CI{i,1}(2,timecurrent)=confup(i);
    end
    
    % factor AxB
    con=conAB;
    [psihat_stat pvalgen pcrit conflow confup]=pbstats(data, con, nboot, alpha);
    
    % passing results into results structure
    results.factor_AxB.pval(:,timecurrent)=pvalgen;
    results.factor_AxB.alpha(:,timecurrent)=pcrit;
    results.factor_AxB.test_stat(:,timecurrent)=psihat_stat;
    
    for i=1:conABcol;
        results.factor_AxB.CI{i,1}(1,timecurrent)=conflow(i);
        results.factor_AxB.CI{i,1}(2,timecurrent)=confup(i);
    end
    
    waitbar(timecurrent/numpnts,h2,sprintf('%12s',[num2str(timecurrent),'/',num2str(numpnts)]))
end
close(h2);

%{
% the string in the sqare brackets is what it appended to the original file name
% add -mat-binary option so octave doesn't add a header or otherwise
% make it irritating to use within matlab
s=ver;
if ~isempty(find(strcmp('Octave',{s.Name})))
    save('-mat-binary', ['results_',fnames_and_identifier.identifier,'.mat'],'results');
    %elseif nargin<7
    %save('bootstrapped_ANOVA_results.mat','results');
    %junk=1;
else
    save('bootstrapped_ANOVA_results.mat','results');
    %save(['results_',fnames_and_identifier.identifier,'.mat'],'results');
end
%}

%if nargin==8
%save([savestring,'_bootstrapped_ANOVA_results.mat'],'results');
%end

end

