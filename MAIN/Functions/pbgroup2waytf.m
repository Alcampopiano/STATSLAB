function [inferential_results sample_results condwaves condfiles_subs condwaves_trim] = pbgroup2waytf(STATS, condfiles, numconds, numpnts, nboot, jlvls, klvls, alpha, nsamp, design, condnames, varargin)
tic
%{

contrasts using con2way default to pooling across one factor. However, you
may want to compare pairwise differences that are not pooled across another factor.
This difference in approach can have a practical difference in the results
and interpretation (See Wilcox, 2013).

This function allows one to either (1) compare all pairwise
differences by pooling, (2) compare all pairwise differences without
pooling, and (3) specify a priori contrasts, which should not include ALL
pairwise differences, but fewer contrasts of course.


%}

nargs = length(varargin);
if round(nargs/2)~=nargs/2
    error('need propertyName/propertyValue pairs for optional inputs')
end

% Set default contrast coefficients for 2-way
% create contrasts for 2way ANOVA (used for multi-comparisons)
[conA conB conAB] = con2way(jlvls, klvls);

% put defaults into a structure;
% options=struct('conA',conA,'conB',conB,'conAB',conAB);

% edited may1st/15
% set default plot options
options.conA=conA;
options.conB=conB;
options.conAB=conAB;
options.FWE='Rom';

% get field names
optionnames = fieldnames(options);

% % check to see which optional args were used and deal with accordingly
% if isempty(varargin);
%     warning('MATLAB:stats',['using default contrasts matrices for factor A, B, and the interaction ' ...
%         'you must specify each one as separate optional input arguments if you want custom contrasts. ' ...
%         'You can leave certain contrasts empty if you want the default comparisons, e.g., [], [1 1 -1 -1]'' ,[]'])
% else
%     % overwrite options stucture with varargin inputs if there are any
%     for i=1:3;
%         if ~isempty(varargin{1}{i})
%             options.(optionnames{i})=varargin{1}{i};
%         end
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

% extract from options structure
conA=options.conA;
conB=options.conB;
conAB=options.conAB;

% used to create proper sizes in results structure
[~, conAcol]=size(conA);
[~, conBcol]=size(conB);
[~, conABcol]=size(conAB);

% load all file names subs X conditions
if isempty(condfiles) % allowing an input for file names
    for i=1:numconds
        tempfname=uigetfile('*.mat',['Select all bootstrapped files in the ', condnames{i}, ' condition'], 'MultiSelect','on');
        condfiles_subs{1,i}(:,1)=tempfname;
    end
    
else
    
    % load a file name that was given that contains the filenames X condition cell array
    condfiles_data=load(condfiles);
    condfields=fieldnames(condfiles_data);
    condfiles_subs=condfiles_data.(condfields{1});
    
end

% try to delete previously made mapped files
warning off
for i=1:length(STATS.condnames);
    delete(['groupboots_',STATS.savestring, '_', STATS.condnames{i},'.map']);
    delete(['mont_groupboots_',STATS.savestring, '_', STATS.condnames{i},'.map']);
    delete(['tempmont_',STATS.savestring, '_', STATS.condnames{i},'.map']);
end

for i=1:conAcol;
    delete(['zsurrogates_',STATS.savestring,'_contrastA',num2str(i),'.map']);
end

for i=1:conBcol;
    delete(['zsurrogates_',STATS.savestring,'_contrastB',num2str(i),'.map']);
end

for i=1:conABcol;
    delete(['zsurrogates_',STATS.savestring,'_contrastAB',num2str(i),'.map']);
end
warning on

%preallocate sizes
[rowconds colconds]=size(condfiles_subs);
condwaves_trim=cell(1,numconds);
%condwaves_trim=zeros(numconds,numpnts);
temparray=zeros(STATS.freqbins,numpnts);
for i=1:numconds;
    %condwaves_trim_origvals{i}=zeros(STATS.freqbins,numpnts);
    condwaves_trim{i}=zeros(STATS.freqbins,numpnts);
    %datacell=cell(1,colconds);
end
%datacell=cell(1,colconds);

% delete from disk the .map files that might have been left over from a
% previous analysia
% delete('.map');
% disp(' **** deleting stray .map files ****');

% this function builds bootstrap inds and writes them to the drive instead
% of holding them in RAM, which makes it scalable (e.g., for 100,000 resamples!)
[rowfile cond_bootvect tmpfname]=bootinds(condfiles_subs,nsamp,design,jlvls,klvls);

% preallocate cell arrays used to accumulate the nsamp CIs
CIlowbootA=cell(conAcol,1);
CIlowbootB=cell(conBcol,1);
CIlowbootAB=cell(conABcol,1);
CIupbootA=cell(conAcol,1);
CIupbootB=cell(conBcol,1);
CIupbootAB=cell(conABcol,1);

%conA, conB, conAB
% this function runs the analysis without resampling from subjects
[sample_results condwaves] = pbgroup2waytf_sample(STATS, numconds, numpnts, nboot, jlvls, klvls, alpha, condfiles_subs, 'FWE', options.FWE, 'conA', conA, 'conB', conB, 'conAB', conAB);

% % build results structure
% results=struct('factor_A',{[]},'factor_B',{[]},'factor_AxB',{[]});
% results.factor_A=struct('contrasts',{conA},'pval',{zeros(conAcol,numpnts)},'alpha',{zeros(conAcol,numpnts)},'test_stat',{zeros(conAcol,numpnts)},'CI',{cell(conAcol,1)}, 'FWE', options.FWE);
%
% for i=1:conAcol;
%     results.factor_A.CI{i,1}=zeros(2,numpnts);
% end
%
% results.factor_B=struct('contrasts',{conB},'pval',{zeros(conBcol,numpnts)},'alpha',{zeros(conBcol,numpnts)},'test_stat',{zeros(conBcol,numpnts)},'CI',{cell(conBcol,1)}, 'FWE', options.FWE);
%
%
% for i=1:conBcol;
%     results.factor_B.CI{i,1}=zeros(2,numpnts);
% end
%
% results.factor_AxB=struct('contrasts',{conAB},'pval',{zeros(conABcol,numpnts)},'alpha',{zeros(conABcol,numpnts)},'test_stat',{zeros(conABcol,numpnts)},'CI',{cell(conABcol,1)}, 'FWE', options.FWE);
%
%
% for i=1:conABcol;
%     results.factor_AxB.CI{i,1}=zeros(2,numpnts);
% end

% make identical results stucture to eventually hold inferential stats
% inferential_results=results;

h1 = waitbar(0,'1','Name','resamples from group','Position',[1100 549 550 40]);
childh1 = get(h1, 'Children');
set(childh1, 'Position',[5 10 538 15]);

h2 = waitbar(0,'1','Name','statistics on frequency bands','Position',[1100 486 550 40]);
childh2 = get(h2, 'Children');
set(childh2, 'Position',[5 10 538 15]);

% band fields
for i=1:STATS.freqbins;
    band_fields{i,1}=['band_', strrep(num2str(STATS.TF_freqs(i)),'.','_')];
end

% bootstrap loop
for bootind=1:nsamp;
    
    % this function builds datacell
    [datacell] = bootgrandaverage(STATS,condfiles_subs,numconds,nboot,numpnts,cond_bootvect,bootind,design,jlvls,klvls);
    
    % write the average of each cell in datacell to a mapped file.
    % get condition waveforms for plotting purposes
    for i=1:numconds;
        mapwrite(mean(datacell{i}.Data.dat,3),['tempmont_',STATS.savestring,'_',STATS.condnames{i},'.map'],'datsize',[STATS.freqbins STATS.timesout nsamp]);
        
        %%%%% dont worry about this, not necessary to have z score condition waves
        % also map z score effect for each monte carlo
        % mapwrite((mean(datacell{i}.Data.dat,3))./(std(datacell{i}.Data.dat,1,3)),['tempmontz_',STATS.savestring,'_',STATS.condnames{i},'.map'],'datsize',[STATS.freqbins STATS.timesout nsamp]);
    end
    
    %arrange the data for the calculations
    rowcell=STATS.nboot;
    
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
            [psihat_stat pvalgen pcrit conflow confup psihat_statz]=pbstats(data, con, nboot, alpha, options.FWE);
            
            % passing results into results structure
            results.(band_fields{bandind}).factor_A.contrasts=conA;
            results.(band_fields{bandind}).factor_A.pval(:,timecurrent)=pvalgen;
            results.(band_fields{bandind}).factor_A.alpha(:,timecurrent)=pcrit;
            results.(band_fields{bandind}).factor_A.test_stat(:,timecurrent)=psihat_stat;
            results.(band_fields{bandind}).factor_A.test_statz(:,timecurrent)=psihat_statz;
            
            for i=1:conAcol;
                results.(band_fields{bandind}).factor_A.CI{i,1}(1,timecurrent)=conflow(i);
                results.(band_fields{bandind}).factor_A.CI{i,1}(2,timecurrent)=confup(i);
            end
            
            % factor B
            con=conB;
            [psihat_stat pvalgen pcrit conflow confup psihat_statz]=pbstats(data, con, nboot, alpha, options.FWE);
            
            % passing results into results structure
            results.(band_fields{bandind}).factor_B.contrasts=conB;
            results.(band_fields{bandind}).factor_B.pval(:,timecurrent)=pvalgen;
            results.(band_fields{bandind}).factor_B.alpha(:,timecurrent)=pcrit;
            results.(band_fields{bandind}).factor_B.test_stat(:,timecurrent)=psihat_stat;
            results.(band_fields{bandind}).factor_B.test_statz(:,timecurrent)=psihat_statz;
            
            for i=1:conBcol;
                results.(band_fields{bandind}).factor_B.CI{i,1}(1,timecurrent)=conflow(i);
                results.(band_fields{bandind}).factor_B.CI{i,1}(2,timecurrent)=confup(i);
            end
            
            % factor AxB
            con=conAB;
            [psihat_stat pvalgen pcrit conflow confup psihat_statz]=pbstats(data, con, nboot, alpha, options.FWE);
            
            % passing results into results structure
            results.(band_fields{bandind}).factor_AxB.contrasts=conAB;
            results.(band_fields{bandind}).factor_AxB.pval(:,timecurrent)=pvalgen;
            results.(band_fields{bandind}).factor_AxB.alpha(:,timecurrent)=pcrit;
            results.(band_fields{bandind}).factor_AxB.test_stat(:,timecurrent)=psihat_stat;
            results.(band_fields{bandind}).factor_AxB.test_statz(:,timecurrent)=psihat_statz;
            
            
            for i=1:conABcol;
                results.(band_fields{bandind}).factor_AxB.CI{i,1}(1,timecurrent)=conflow(i);
                results.(band_fields{bandind}).factor_AxB.CI{i,1}(2,timecurrent)=confup(i);
            end
            
            
        end
        waitbar(bandind/STATS.freqbins,h2,sprintf('%12s',[num2str(bandind),'/',num2str(STATS.freqbins)]))
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%This is where we extract only what we need from each bootstrap%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ext=1:conAcol
        for bandind=1:STATS.freqbins;
            
            %[~,tmpfname_tmp]=fileparts(tempname);
            %tmpfname_diff{ext}=tmpfname_tmp;
            temparray(bandind,:)=results.(band_fields{bandind}).factor_A.test_statz(ext,:);
            %mapwrite(results.(band_fields{bandind}).factor_A.test_statz(ext,:),[tmpfname_diff{ext},'.map'],'datsize',[STATS.freqbins numpnts nsamp]);
            
            
        end
        
        % here are the z surrogates that we do stats on, the arrray holds
        % all effect sizes for every band and timepoint for each monte carlo
        mapwrite(temparray,['zsurrogates_',STATS.savestring,'_contrastA',num2str(ext),'.map'],'datsize',[STATS.freqbins numpnts nsamp]);
        %mapwrite(temparray,[tmpfname_diff{ext},'.map'],'datsize',[STATS.freqbins numpnts nsamp]);
    end
    
    for ext=1:conBcol
        for bandind=1:STATS.freqbins;
            
            %[~,tmpfname_tmp]=fileparts(tempname);
            %tmpfname_diff{ext}=tmpfname_tmp;
            temparray(bandind,:)=results.(band_fields{bandind}).factor_B.test_statz(ext,:);
            %mapwrite(results.(band_fields{bandind}).factor_A.test_statz(ext,:),[tmpfname_diff{ext},'.map'],'datsize',[STATS.freqbins numpnts nsamp]);
            
            
        end
        
        % here are the z surrogates that we do stats on, the arrray holds
        % all effect sizes for every band and timepoint for each monte carlo
        mapwrite(temparray,['zsurrogates_',STATS.savestring,'_contrastB',num2str(ext),'.map'],'datsize',[STATS.freqbins numpnts nsamp]);
        %mapwrite(temparray,[tmpfname_diff{ext},'.map'],'datsize',[STATS.freqbins numpnts nsamp]);
    end
    
    for ext=1:conABcol
        for bandind=1:STATS.freqbins;
            
            %[~,tmpfname_tmp]=fileparts(tempname);
            %tmpfname_diff{ext}=tmpfname_tmp;
            temparray(bandind,:)=results.(band_fields{bandind}).factor_AxB.test_statz(ext,:);
            %mapwrite(results.(band_fields{bandind}).factor_A.test_statz(ext,:),[tmpfname_diff{ext},'.map'],'datsize',[STATS.freqbins numpnts nsamp]);
            
            
        end
        
        % here are the z surrogates that we do stats on, the arrray holds
        % all effect sizes for every band and timepoint for each monte carlo
        mapwrite(temparray,['zsurrogates_',STATS.savestring,'_contrastAB',num2str(ext),'.map'],'datsize',[STATS.freqbins numpnts nsamp]);
        %mapwrite(temparray,[tmpfname_diff{ext},'.map'],'datsize',[STATS.freqbins numpnts nsamp]);
    end
    
    %%% diffwaves should be iteratively written to drive - mem map
    %%% stealing differnce waves in order to calculate "real" CIs
%     for ext=1:conAcol
%         diffwaveA{ext,1}(bootind,:)=results.factor_A.test_stat(ext,:);
%     end
%     
%     for ext=1:conBcol
%         diffwaveB{ext,1}(bootind,:)=results.factor_B.test_stat(ext,:);
%     end
%     
%     for ext=1:conABcol
%         diffwaveAB{ext,1}(bootind,:)=results.factor_AxB.test_stat(ext,:);
%     end
    
    waitbar(bootind/nsamp,h1,sprintf('%12s',[num2str(bootind),'/',num2str(nsamp)]))
end


% put lower and upper bounds into a cell, lowers 1st, uppers second
%CIA={CIlowbootA,CIupbootA};
%CIB={CIlowbootB,CIupbootB};
%CIAB={CIlowbootAB,CIupbootAB};


%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
% Mother CIs, the way they are here, do NOT need to be cacultaed, instead,
% just find the CIs from the bootstrapped average differecnce wave
% (psihat_stat). You should do that here so that you dont have to keep
% calculating them after running this function, which is getting annoying
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
% calculate the 95% "mother" CIs based on what you have accumulated
%[CIlowA CIupA]=grand_pb_bootCI(CIA, conAcol, numpnts, nsamp, alpha);
%[CIlowB CIupB]=grand_pb_bootCI(CIB, conBcol, numpnts, nsamp, alpha);
%[CIlowAB CIupAB]=grand_pb_bootCI(CIAB, conABcol, numpnts, nsamp, alpha);


% map write the CI arrays, might be cool to see them at some point
% for ext=1:conAcol
%     fidm=mapwrite(CIlowbootA{ext,1},['CIlowbootA',num2str(ext),'.map'],'datsize',[nsamp numpnts]);
%     fidm=mapwrite(CIupbootA{ext,1},['CIupbootA',num2str(ext),'.map'],'datsize',[nsamp numpnts]);
%     fidm=mapwrite(diffwaveA{ext,1},['diffwaveA',num2str(ext),'.map'],'datsize',[nsamp numpnts]);
%
% end
%
% for ext=1:conBcol
%     fidm=mapwrite(CIlowbootB{ext,1},['CIlowbootB',num2str(ext),'.map'],'datsize',[nsamp numpnts]);
%     fidm=mapwrite(CIupbootB{ext,1},['CIupbootB',num2str(ext),'.map'],'datsize',[nsamp numpnts]);
%     fidm=mapwrite(diffwaveB{ext,1},['diffwaveB',num2str(ext),'.map'],'datsize',[nsamp numpnts]);
% end
%
% for ext=1:conABcol
%     fidm=mapwrite(CIlowbootAB{ext,1},['CIlowbootAB',num2str(ext),'.map'],'datsize',[nsamp numpnts]);
%     fidm=mapwrite(CIupbootAB{ext,1},['CIupbootAB',num2str(ext),'.map'],'datsize',[nsamp numpnts]);
%     fidm=mapwrite(diffwaveAB{ext,1},['diffwaveAB',num2str(ext),'.map'],'datsize',[nsamp numpnts]);
% end

close(h1,h2);
%%%%%%%%%%%%%%%%%% inferential statistics %%%%%%%%%%%%%%%%%

% get condition waveforms to plot
for i=1:numconds;
    
    %     %%%%%% dont worry about this one, make var equal to the regular unit waveforms
    %     % also map z score effect for each monte carlo
    %     dat_tempmontz=mapread(['tempmontz_',STATS.savestring, '_', STATS.condnames{i},'.map'], 'dat','datsize',[STATS.freqbins STATS.timesout nsamp]);
    %     condwaves_trim{i}=mean(dat_tempmontz.Data.dat,3);
    
    dat_tempmont=mapread(['tempmont_',STATS.savestring, '_', STATS.condnames{i},'.map'], 'dat');
    condwaves_trim{i}=mean(dat_tempmont.Data.dat,3);
end

%preallocate samll data arrays
data_A=zeros(nsamp,conAcol);
data_B=zeros(nsamp,conBcol);
data_AB=zeros(nsamp,conABcol);


% access the big difference wave arrays
for i=1:conAcol
    
    diffdata.(['A',num2str(i)])=mapread(['zsurrogates_',STATS.savestring,'_contrastA',num2str(ext),'.map'],'dat');
    
end

for i=1:conBcol
    diffdata.(['B',num2str(i)])=mapread(['zsurrogates_',STATS.savestring,'_contrastB',num2str(ext),'.map'],'dat');
end

for i=1:conABcol
    diffdata.(['AB',num2str(i)])=mapread(['zsurrogates_',STATS.savestring,'_contrastAB',num2str(ext),'.map'],'dat');
end

% waitbar for final stages, doing inferential stats
h3 = waitbar(0,'1','Name','inferential statistics on frequency bands','Position',[1100 486 550 40]);
childh3 = get(h3, 'Children');
set(childh3, 'Position',[5 10 538 15]);

for bandind=1:STATS.freqbins;
    
    % loop for stats at each timepoint
    for timecurrent=1:numpnts;
        
        % factor A
        con=conA;
        for i=1:conAcol;
            data_A(:,i)=diffdata.(['A',num2str(i)]).Data.dat(bandind,timecurrent,:);
        end
        
        [psihat_stat pvalgen pcrit conflow confup]=pbstats_diff(data_A, con, nsamp, alpha, options.FWE);
        
        % passing results into results structure
        inferential_results.(band_fields{bandind}).factor_A.contrasts=conA;
        inferential_results.(band_fields{bandind}).factor_A.pval(:,timecurrent)=pvalgen;
        inferential_results.(band_fields{bandind}).factor_A.alpha(:,timecurrent)=pcrit;
        inferential_results.(band_fields{bandind}).factor_A.test_stat(:,timecurrent)=psihat_stat;
        
        for i=1:conAcol;
            inferential_results.(band_fields{bandind}).factor_A.CI{i,1}(1,timecurrent)=conflow(i);
            inferential_results.(band_fields{bandind}).factor_A.CI{i,1}(2,timecurrent)=confup(i);
        end
        
        % factor B
        con=conB;
        for i=1:conBcol;
            data_B(:,i)=diffdata.(['B',num2str(i)]).Data.dat(bandind,timecurrent,:);
        end
        
        [psihat_stat pvalgen pcrit conflow confup]=pbstats_diff(data_B, con, nsamp, alpha, options.FWE);
        
        % passing results into results structure
        inferential_results.(band_fields{bandind}).factor_B.contrasts=conB;
        inferential_results.(band_fields{bandind}).factor_B.pval(:,timecurrent)=pvalgen;
        inferential_results.(band_fields{bandind}).factor_B.alpha(:,timecurrent)=pcrit;
        inferential_results.(band_fields{bandind}).factor_B.test_stat(:,timecurrent)=psihat_stat;
        
        for i=1:conBcol;
            inferential_results.(band_fields{bandind}).factor_B.CI{i,1}(1,timecurrent)=conflow(i);
            inferential_results.(band_fields{bandind}).factor_B.CI{i,1}(2,timecurrent)=confup(i);
        end
        
        % factor A
        con=conAB;
        for i=1:conABcol;
            data_AB(:,i)=diffdata.(['AB',num2str(i)]).Data.dat(bandind,timecurrent,:);
        end
        
        [psihat_stat pvalgen pcrit conflow confup]=pbstats_diff(data_AB, con, nsamp, alpha, options.FWE);
        
        % passing results into results structure
        inferential_results.(band_fields{bandind}).factor_AxB.contrasts=conAB;
        inferential_results.(band_fields{bandind}).factor_AxB.pval(:,timecurrent)=pvalgen;
        inferential_results.(band_fields{bandind}).factor_AxB.alpha(:,timecurrent)=pcrit;
        inferential_results.(band_fields{bandind}).factor_AxB.test_stat(:,timecurrent)=psihat_stat;
        
        
        for i=1:conABcol;
            inferential_results.(band_fields{bandind}).factor_AxB.CI{i,1}(1,timecurrent)=conflow(i);
            inferential_results.(band_fields{bandind}).factor_AxB.CI{i,1}(2,timecurrent)=confup(i);
        end
        
        %waitbar(timecurrent/numpnts,h3,sprintf('%12s',[num2str(timecurrent),'/',num2str(numpnts)]))
    end
    waitbar(bandind/STATS.freqbins,h3,sprintf('%12s',[num2str(bandind),'/',num2str(STATS.freqbins)]))
end

% edit may 8th/15
% clean temporary mapped files
if iscell(tmpfname)
    for i=1:length(tempfname);
        delete(tempfname{i});
    end
    
else
    delete(tempfname);
end

close(h3)
end


