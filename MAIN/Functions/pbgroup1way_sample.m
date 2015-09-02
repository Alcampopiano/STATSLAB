function [results condwaves] = pbgroup1way_sample(STATS,numconds, numpnts, nboot, jlvls, alpha, condfiles_subs, varargin)
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
[conA] = con1way(jlvls);

% put defaults into a structure;
% options=struct('conA',conA);

% edited may1st/15
% set default plot options
options.conA=conA;
options.FWE='Rom';

% get field names
optionnames = fieldnames(options);

% check to see which optional args were used and deal with accordingly
% if isempty(varargin);
%     warning('MATLAB:stats',['Using default contrasts matrix. You must specify one if you want a custom contrast. ' ...
%         ' e.g., [1 -1 0; 1 0 -1]'''])
% else
%     % overwrite options stucture with varargin inputs if there are any
%     
%     if ~isempty(varargin{1})
%         options.(optionnames{1})=varargin{1};
%     end
%     
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

% extract from options structure
conA=options.conA;

% used to create proper sizes in results structure
[~, conAcol]=size(conA);

% load data
if isempty(condfiles_subs);
    for i=1:numconds
        tempfname=uigetfile('*.mat',['Select all subject condition ', num2str(i), ' files'], 'MultiSelect','on');
        condfiles_subs{1,i}(:,1)=tempfname;
    end
    [datacell] = grandaverage(STATS,nboot,numpnts,condfiles_subs);
else
    [datacell] = grandaverage(STATS,nboot,numpnts,condfiles_subs);
end

% get condition waveforms for plotting purposes
for i=1:length(datacell);
    condwaves(i,:)=mean(datacell{i},1);
end

%preallocate sizes
[rowconds colconds]=size(condfiles_subs);

% build results structure
results=struct('factor_A',{[]});
results.factor_A=struct('contrasts',{conA},'pval',{zeros(conAcol,numpnts)},'alpha',{zeros(conAcol,numpnts)},'test_stat',{zeros(conAcol,numpnts)},'CI',{cell(conAcol,1)}, 'FWE', options.FWE);
%results.factor_A.z_scores=struct('z_avg',{zeros(conAcol,numpnts)},'CI_z',{cell(conAcol,1)});

for i=1:conAcol;
    results.factor_A.CI{i,1}=zeros(2,numpnts);
    %results.factor_A.z_scores.CI_z{i,1}=zeros(2,numpnts);
end

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
    [psihat_stat, pvalgen, pcrit, conflow, confup, psihat_statz]=pbstats(data, con, nboot, alpha, options.FWE);
    
    % passing results into results structure
    results.factor_A.pval(:,timecurrent)=pvalgen;
    results.factor_A.alpha(:,timecurrent)=pcrit;
    results.factor_A.test_stat(:,timecurrent)=psihat_stat;
    
    for i=1:conAcol;
        results.factor_A.CI{i,1}(1,timecurrent)=conflow(i);
        results.factor_A.CI{i,1}(2,timecurrent)=confup(i);
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

