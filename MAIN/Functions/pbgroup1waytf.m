function [sample_results condwaves condfiles_subs] = pbgroup1waytf(STATS, condfiles, numconds, numpnts, nboot, jlvls, alpha, nsamp, design, condnames, varargin)
tic


nargs = length(varargin);
if round(nargs/2)~=nargs/2
    error('need propertyName/propertyValue pairs for optional inputs')
end

% Set default contrast coefficients for 2-way
% create contrasts for 1way ANOVA (used for multi-comparisons)
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
%         if ~isempty(varargin{1})
%         options.(optionnames{1})=varargin{1}{1};
%         end
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

% load all file names subs X conditions
if isempty(condfiles) % allowing an input for file names
    for i=1:numconds
        tempfname=uigetfile('*.map',['Select all bootstrapped files in the ', condnames{i}, ' condition'], 'MultiSelect','on');
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
warning on

%preallocate sizes
[rowconds colconds]=size(condfiles_subs);
temparray=zeros(STATS.freqbins,numpnts);
condwaves_trim=cell(1,numconds);
for i=1:numconds;
    %condwaves_trim_origvals{i}=zeros(STATS.freqbins,numpnts);
    condwaves_trim{i}=zeros(STATS.freqbins,numpnts);
    %datacell=cell(1,colconds);
end

% delete from disk the .map files that might have been left over from a
% previous analysia
% delete('*boot*.map','*wave*.map');
% disp(' **** deleting stray .map files ****');

% this function builds bootstrap inds and writes them to the drive instead
% of holding them in RAM, which makes it scalable (e.g., for 100,000 resamples!)
[rowfile cond_bootvect tmpfname]=bootinds(condfiles_subs,nsamp,design,jlvls);

% preallocate cell arrays used to accumulate the nsamp CIs
CIlowbootA=cell(conAcol,1);
CIupbootA=cell(conAcol,1);

% this function runs the analysis without resampling from subjects
[sample_results condwaves] = pbgroup1waytf_sample(STATS, numconds, STATS.timesout, nboot, jlvls, alpha, condfiles_subs, 'FWE', options.FWE, 'conA', conA);

% % % build results structure
% % results=struct('factor_A',{[]});
% % results.factor_A=struct('contrasts',{conA},'pval',{zeros(conAcol,numpnts)},'alpha',{zeros(conAcol,numpnts)},'test_stat',{zeros(conAcol,numpnts)},'CI',{cell(conAcol,1)}, 'FWE', options.FWE);
% %
% % for i=1:conAcol;
% %     results.factor_A.CI{i,1}=zeros(2,numpnts);
% % end
%
% % make identical results stucture to eventually hold inferential stats
% %inferential_results=results;
%
% % load and arrange data
% h1 = waitbar(0,'1','Name','resamples from group','Position',[1100 549 550 40]);
% childh1 = get(h1, 'Children');
% set(childh1, 'Position',[5 10 538 15]);
%
% h2 = waitbar(0,'1','Name','statistics on frequency bands','Position',[1100 486 550 40]);
% childh2 = get(h2, 'Children');
% set(childh2, 'Position',[5 10 538 15]);
%
% % band fields
% for i=1:STATS.freqbins;
%     band_fields{i,1}=['band_', strrep(num2str(STATS.TF_freqs(i)),'.','_')];
% end
%
% % build temporary file names
% % for ext=1:conAcol
% %     [~,tmpfname_tmp]=fileparts(tempname);
% %     tmpfname_diff{ext}=tmpfname_tmp;
% % end
%
% % bootstrap loop
% for bootind=1:nsamp;
%
%     % this function builds datacell
%     [datacell] = bootgrandaverage(STATS,condfiles_subs,numconds,nboot,numpnts,cond_bootvect,bootind,design,jlvls);
%
%     % write the average of each cell in datacell to a mapped file.
%     % get condition waveforms for plotting purposes
%     for i=1:numconds;
%         mapwrite(mean(datacell{i}.Data.dat,3),['tempmont_',STATS.savestring,'_',STATS.condnames{i},'.map'],'datsize',[STATS.freqbins STATS.timesout nsamp]);
%
%         %%%%% dont worry about this, not necessary to have z score condition waves
%         % also map z score effect for each monte carlo
%         % mapwrite((mean(datacell{i}.Data.dat,3))./(std(datacell{i}.Data.dat,1,3)),['tempmontz_',STATS.savestring,'_',STATS.condnames{i},'.map'],'datsize',[STATS.freqbins STATS.timesout nsamp]);
%     end
%
%     %arrange the data for the calculations
%     rowcell=STATS.nboot;
%
%     for bandind=1:STATS.freqbins;
%
%         % loop for stats at each timepoint
%         for timecurrent=1:STATS.timesout;
%
%             % reset data to zeros after every calculation at each timepoint
%             data=zeros(rowcell,numconds);
%
%             % arrange data into a matrix with subs (or single subject boot samples) X conditions
%             for condcurrent=1:colconds;
%                 % data(:,condcurrent)=datacell{1,condcurrent}(:,timecurrent);
%                 data(:,condcurrent)=datacell{1,condcurrent}.Data.dat(bandind,timecurrent,:);
%             end
%
%
%             % factor A
%             con=conA;
%             [psihat_stat pvalgen pcrit conflow confup psihat_statz]=pbstats(data, con, nboot, alpha, options.FWE);
%
%             % passing results into results structure
%             results.(band_fields{bandind}).factor_A.contrasts=conA;
%             results.(band_fields{bandind}).factor_A.pval(:,timecurrent)=pvalgen;
%             results.(band_fields{bandind}).factor_A.alpha(:,timecurrent)=pcrit;
%             results.(band_fields{bandind}).factor_A.test_stat(:,timecurrent)=psihat_stat;
%             results.(band_fields{bandind}).factor_A.test_statz(:,timecurrent)=psihat_statz;
%
%             for i=1:conAcol;
%                 results.(band_fields{bandind}).factor_A.CI{i,1}(1,timecurrent)=conflow(i);
%                 results.(band_fields{bandind}).factor_A.CI{i,1}(2,timecurrent)=confup(i);
%             end
%
%             %waitbar(timecurrent/numpnts,h2,sprintf('%12s',[num2str(timecurrent),'/',num2str(numpnts)]))
%         end
%         waitbar(bandind/STATS.freqbins,h2,sprintf('%12s',[num2str(bandind),'/',num2str(STATS.freqbins)]))
%     end
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%This is where we extract only what we need from each bootstrap%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     %     for ext=1:conAcol
%     %         CIlowbootA{ext,1}(bootind,:)=results.factor_A.CI{ext,1}(1,:);
%     %         CIupbootA{ext,1}(bootind,:)=results.factor_A.CI{ext,1}(2,:);
%     %     end
%
%     %%% get effect sizes for each bootstrap test
%     %     for ext=1:conAcol
%     %
%     %         % avg difference
%     %         %diffwaveA{ext,1}(bootind,:)=results.factor_A.test_stat(ext,:);
%     %
%     %         % z effect
%     %         diffwaveA{ext,1}(bootind,:)=results.factor_A.test_statz(ext,:);
%     %     end
%     for ext=1:conAcol
%         for bandind=1:STATS.freqbins;
%
%             %[~,tmpfname_tmp]=fileparts(tempname);
%             %tmpfname_diff{ext}=tmpfname_tmp;
%             temparray(bandind,:)=results.(band_fields{bandind}).factor_A.test_statz(ext,:);
%             %mapwrite(results.(band_fields{bandind}).factor_A.test_statz(ext,:),[tmpfname_diff{ext},'.map'],'datsize',[STATS.freqbins numpnts nsamp]);
%
%
%         end
%
%         % here are the z surrogates that we do stats on, the arrray holds
%         % all effect sizes for every band and timepoint for each monte carlo
%         mapwrite(temparray,['zsurrogates_',STATS.savestring,'_contrastA',num2str(ext),'.map'],'datsize',[STATS.freqbins numpnts nsamp]);
%         %mapwrite(temparray,[tmpfname_diff{ext},'.map'],'datsize',[STATS.freqbins numpnts nsamp]);
%     end
%
%     waitbar(bootind/nsamp,h1,sprintf('%12s',[num2str(bootind),'/',num2str(nsamp)]))
% end
%
%
% % put lower and upper bounds into a cell, lowers 1st, uppers second
% %CIA={CIlowbootA,CIupbootA};
% %CIB={CIlowbootB,CIupbootB};
% %CIAB={CIlowbootAB,CIupbootAB};
%
% % map write the CI arrays, might be cool to see them at some point
% % for ext=1:conAcol
% %     [~,tmpfname_tmp]=fileparts(tempname);
% %     tmpfname_diff{ext}=tmpfname_tmp;
% %     mapwrite(diffwaveA{ext,1},[tmpfname_diff{ext},'.map'],'datsize',[nsamp numpnts]);
% %
% % end
%
% close(h1,h2);
% %%%%%%%%%%%%%%%%%% inferential statistics %%%%%%%%%%%%%%%%%
%
%
% % get condition waveforms to plot
% for i=1:numconds;
%
%     %     %%%%%% dont worry about this one, make var equal to the regular unit waveforms
%     %     % also map z score effect for each monte carlo
%     %     dat_tempmontz=mapread(['tempmontz_',STATS.savestring, '_', STATS.condnames{i},'.map'], 'dat','datsize',[STATS.freqbins STATS.timesout nsamp]);
%     %     condwaves_trim{i}=mean(dat_tempmontz.Data.dat,3);
%
%     dat_tempmont=mapread(['tempmont_',STATS.savestring, '_', STATS.condnames{i},'.map'], 'dat');
%     condwaves_trim{i}=mean(dat_tempmont.Data.dat,3);
% end
%
% %preallocate samll data arrays
% data_A=zeros(nsamp,conAcol);
%
% % access the big difference wave arrays
% for i=1:conAcol
%     %diffdata.(['A',num2str(i)])=mapread([tmpfname_diff{i},'.map'],'dat');
%     diffdata.(['A',num2str(i)])=mapread(['zsurrogates_',STATS.savestring,'_contrastA',num2str(ext),'.map'],'dat');
% end
%
% % waitbar for final stages, doing inferential stats
% h3 = waitbar(0,'1','Name','inferential statistics on frequency bands','Position',[1100 486 550 40]);
% childh3 = get(h3, 'Children');
% set(childh3, 'Position',[5 10 538 15]);
%
%
% for bandind=1:STATS.freqbins;
%
%     % loop for stats at each timepoint
%     for timecurrent=1:numpnts;
%
%         % factor A
%         con=conA;
%         for i=1:conAcol;
%             data_A(:,i)=diffdata.(['A',num2str(i)]).Data.dat(bandind,timecurrent,:);
%         end
%
%         [psihat_stat pvalgen pcrit conflow confup]=pbstats_diff(data_A, con, nsamp, alpha, options.FWE);
%
%         % passing results into results structure
%         inferential_results.(band_fields{bandind}).factor_A.contrasts=conA;
%         inferential_results.(band_fields{bandind}).factor_A.pval(:,timecurrent)=pvalgen;
%         inferential_results.(band_fields{bandind}).factor_A.alpha(:,timecurrent)=pcrit;
%         inferential_results.(band_fields{bandind}).factor_A.test_stat(:,timecurrent)=psihat_stat;
%
%         for i=1:conAcol;
%             inferential_results.(band_fields{bandind}).factor_A.CI{i,1}(1,timecurrent)=conflow(i);
%             inferential_results.(band_fields{bandind}).factor_A.CI{i,1}(2,timecurrent)=confup(i);
%         end
%
%         %waitbar(timecurrent/numpnts,h3,sprintf('%12s',[num2str(timecurrent),'/',num2str(numpnts)]))
%     end
%     waitbar(bandind/STATS.freqbins,h3,sprintf('%12s',[num2str(bandind),'/',num2str(STATS.freqbins)]))
% end

% edit may 8th/15
% clean temporary mapped files
try
    if iscell(tmpfname)
        for i=1:length(tmpfname);
            delete([tmpfname{i}, '.map']);
        end
        
    else
        delete([tmpfname, '.map']);
    end
catch
end

% if iscell(tmpfname_diff)
%     for i=1:length(tmpfname_diff);
%         % delete([tmpfname_CIup{i}, '.map']);
%         % delete([tmpfname_CIlow{i}, '.map']);
%         delete([tmpfname_diff{i}, '.map']);
%     end
%
% else
%     %    delete([tmpfname_CIup, '.map']);
%     %    delete([tmpfname_CIlow, '.map']);
%     delete([tmpfname_diff, '.map']);
% end

%close(h3)
end



