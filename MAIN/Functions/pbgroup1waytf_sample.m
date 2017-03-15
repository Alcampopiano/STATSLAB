function [results condwaves] = pbgroup1waytf_sample(STATS,numconds, numpnts, nboot, jlvls, alpha, condfiles_subs, varargin)
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
options.FWE='benhoch';

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

% temp
if strcmp(options.FWE, 'benhoch');
    results.factor_A.FWE=options.FWE; % temporary
end

% used to create proper sizes in results structure
[~, conAcol]=size(conA);

% load data
if isempty(condfiles_subs);
    for i=1:numconds
        tempfname=uigetfile('*.map',['Select all subject condition ', num2str(i), ' files'], 'MultiSelect','on');
        condfiles_subs{1,i}(:,1)=tempfname;
    end
    [datacell] = grandaverage(STATS,nboot,numpnts,condfiles_subs);
else
    [datacell] = grandaverage(STATS,nboot,numpnts,condfiles_subs);
end

% get condition waveforms for plotting purposes
% for TF, set condition waveforms to a file name (string), rather than data, this
% string was defined in grandaverage.m
for i=1:length(datacell);
    condwaves{i,:}=['group_TFwaves_', STATS.savestring,'_',STATS.condnames{i},'.,mat'];
end

%preallocate sizes
[rowconds colconds]=size(condfiles_subs);

% % build results structure
% results=struct('factor_A',{[]});
% results.factor_A=struct('contrasts',{conA},'pval',{zeros(conAcol,numpnts)},'alpha',{zeros(conAcol,numpnts)},'test_stat',{zeros(conAcol,numpnts)},'CI',{cell(conAcol,1)}, 'FWE', options.FWE);
%
% for i=1:conAcol;
%     results.factor_A.CI{i,1}=zeros(2,numpnts);
% end

h1 = waitbar(0,'1','Name','statistics across frequency windows','Position',[1100 549 550 40]);
childh1 = get(h1, 'Children');
set(childh1, 'Position',[5 10 538 15]);

h2 = waitbar(0,'1','Name','statistics on frequency bands','Position',[1100 486 550 40]);
childh2 = get(h2, 'Children');
set(childh2, 'Position',[5 10 538 15]);

%%%
% preallocating difference array
if strcmp(options.FWE, 'benhoch')
    results.factor_A.diffs=cell(conAcol,1);
    for i=1:conAcol;
       results.factor_A.diffs{i,1}=zeros(STATS.nboot,STATS.timesout,STATS.freqbins);
    end
end
%%%


%arrange the data for the calculations
rowcell=STATS.nboot;

% band fields
for i=1:STATS.freqbins;
    band_fields{i,1}=['band_', strrep(num2str(STATS.TF_freqs(i)),'.','_')];
end

for bandind=1:STATS.freqbins;
    
    % loop for stats at each timepoint
    for timecurrent=1:STATS.timesout;
        
        % reset data to zeros after every calculation at each timepoint
        data=zeros(rowcell,numconds);
        
        % arrange data into a matrix with subs (or single subject boot samples) X conditions
        for condcurrent=1:colconds;
            % data(:,condcurrent)=datacell{1,condcurrent}(:,timecurrent);
            data(:,condcurrent)=datacell{1,condcurrent}.Data.dat(bandind,timecurrent,:);
        end
        
        
        % factor A
        con=conA;
        [psihat, psihat_stat, pvalgen, pcrit, conflow, confup, psihat_statz]=pbstats(data, con, nboot, alpha, options.FWE);
        
        % passing results into results structure
        results.(band_fields{bandind}).factor_A.contrasts=conA;
        results.(band_fields{bandind}).factor_A.pval(:,timecurrent)=pvalgen;
        results.(band_fields{bandind}).factor_A.alpha(:,timecurrent)=pcrit;
        results.(band_fields{bandind}).factor_A.test_stat(:,timecurrent)=psihat_stat;
        
        if strcmp(options.FWE, 'benhoch')
            for i=1:conAcol;
                results.factor_A.diffs{i,1}(:,timecurrent,bandind)=psihat(:,i);
            end
        end
        
        for i=1:conAcol;            
            results.(band_fields{bandind}).factor_A.CI{i,1}(1,timecurrent)=conflow(i);
            results.(band_fields{bandind}).factor_A.CI{i,1}(2,timecurrent)=confup(i);
        end
        
        waitbar(timecurrent/numpnts,h1,sprintf('%12s',[num2str(timecurrent),'/',num2str(numpnts)]))
    end
    
    % add FWE option to structure
    results.(band_fields{bandind}).factor_A.FWE=options.FWE;
    waitbar(bandind/STATS.freqbins,h2,sprintf('%12s',[num2str(bandind),'/',num2str(STATS.freqbins)]))
      
end

close(h1,h2);

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

