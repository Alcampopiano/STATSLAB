function [inferential_results sample_results condwaves condfiles_subs condwaves_trim] = pbgroup1way(STATS,condfiles, numconds, numpnts, nboot, jlvls, alpha, nsamp, design, condnames, varargin)
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
        tempfname=uigetfile('*.mat',['Select all bootstrapped files in the ', condnames{i}, ' condition'], 'MultiSelect','on');
        condfiles_subs{1,i}(:,1)=tempfname;
    end
    
else
    
    % load a file name that was given that contains the filenames X condition cell array
    condfiles_data=load(condfiles);
    condfields=fieldnames(condfiles_data);
    condfiles_subs=condfiles_data.(condfields{1});
    
end


%load('condfiles_subs.mat') 

%preallocate sizes
[rowconds colconds]=size(condfiles_subs);
condwaves_trim_gather=cell(1,numconds);
condwaves_trim=zeros(numconds,numpnts);
%datacell=cell(1,colconds);

% delete from disk the .map files that might have been left over from a
% previous analysia
 delete('*boot*.map','*wave*.map');
% disp(' **** deleting stray .map files ****');

% this function builds bootstrap inds and writes them to the drive instead
% of holding them in RAM, which makes it scalable (e.g., for 100,000 resamples!)
[rowfile cond_bootvect tmpfname]=bootinds(condfiles_subs,nsamp,design,jlvls);
    
% preallocate cell arrays used to accumulate the nsamp CIs
CIlowbootA=cell(conAcol,1);
CIupbootA=cell(conAcol,1);

% this function runs the analysis without resampling from subjects
[sample_results condwaves] = pbgroup1way_sample(STATS, numconds, numpnts, nboot, jlvls, alpha, condfiles_subs, 'FWE', options.FWE, 'conA', conA);

% build results structure
results=struct('factor_A',{[]});
results.factor_A=struct('contrasts',{conA},'pval',{zeros(conAcol,numpnts)},'alpha',{zeros(conAcol,numpnts)},'test_stat',{zeros(conAcol,numpnts)},'CI',{cell(conAcol,1)}, 'FWE', options.FWE);

for i=1:conAcol;
    results.factor_A.CI{i,1}=zeros(2,numpnts);
end

% make identical results stucture to eventually hold inferential stats
inferential_results=results;

% load and arrange data
h1 = waitbar(0,'1','Name','resamples from group','Position',[1100 549 550 40]);
childh1 = get(h1, 'Children');
set(childh1, 'Position',[5 10 538 15]);

h2 = waitbar(0,'1','Name','stats across time','Position',[1100 486 550 40]);
childh2 = get(h2, 'Children');
set(childh2, 'Position',[5 10 538 15]);

% bootstrap loop
for bootind=1:nsamp;
   
    % this function builds datacell
    [datacell] = bootgrandaverage(STATS,condfiles_subs,numconds,nboot,numpnts,cond_bootvect,bootind,design,jlvls);
                                    
    % get condition waveforms for plotting purposes
    %%%%%%%%%%%%%%%
    % why bother doing this gathering? Maybe just take mean from the
    % analyses with al subjects, or maybe its fine.
    %%%%%%%%%%%%%%
    for i=1:numconds;
        condwaves_trim_gather{i}(bootind,:)=mean(datacell{i},1);
    end
    
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
        [psihat_stat pvalgen pcrit conflow confup]=pbstats(data, con, nboot, alpha, options.FWE);
        
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%This is where we extract only what we need from each bootstrap%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ext=1:conAcol
        CIlowbootA{ext,1}(bootind,:)=results.factor_A.CI{ext,1}(1,:);
        CIupbootA{ext,1}(bootind,:)=results.factor_A.CI{ext,1}(2,:);
    end

    %%% diffwaves should be iteratively written to drive - mem map
    %%% stealing differnce waves in order to calculate "real" CIs
    for ext=1:conAcol
        diffwaveA{ext,1}(bootind,:)=results.factor_A.test_stat(ext,:);
    end
   
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
for ext=1:conAcol
    
    [~,tmpfname_tmp]=fileparts(tempname);
    tmpfname_CIup{ext}=tmpfname_tmp;
    fidm=mapwrite(CIlowbootA{ext,1},[tmpfname_CIup{ext},'.map'],'datsize',[nsamp numpnts]);
    
    [~,tmpfname_tmp]=fileparts(tempname);
    tmpfname_CIlow{ext}=tmpfname_tmp;
    fidm=mapwrite(CIupbootA{ext,1},[tmpfname_CIlow{ext},'.map'],'datsize',[nsamp numpnts]);
    
    [~,tmpfname_tmp]=fileparts(tempname);
    tmpfname_diff{ext}=tmpfname_tmp;
    fidm=mapwrite(diffwaveA{ext,1},[tmpfname_diff{ext},'.map'],'datsize',[nsamp numpnts]);
    
end

close(h1,h2);
%%%%%%%%%%%%%%%%%% inferential statistics %%%%%%%%%%%%%%%%%

% get condition waveforms to plot
for i=1:numconds;
    condwaves_trim(i,:)=mean(condwaves_trim_gather{i},1);
end

%preallocate samll data arrays
data_A=zeros(nsamp,conAcol);

% access the big difference wave arrays
for i=1:conAcol
    diffdata.(['A',num2str(i)])=mapread([tmpfname_diff{i},'.map'],'dat','datsize',[nsamp numpnts]);
end

% waitbar for final stages, doing inferential stats
h3 = waitbar(0,'1','Name','inferential statistics','Position',[1100 486 550 40]);
childh3 = get(h3, 'Children');
set(childh3, 'Position',[5 10 538 15]);

% loop for stats at each timepoint
for timecurrent=1:numpnts;
    
    % factor A
    con=conA;
    for i=1:conAcol;
        data_A(:,i)=diffdata.(['A',num2str(i)]).Data.dat(:,timecurrent);
    end
    
    [psihat_stat pvalgen pcrit conflow confup]=pbstats_diff(data_A, con, nsamp, alpha, options.FWE);
    
    % passing results into results structure
    inferential_results.factor_A.pval(:,timecurrent)=pvalgen;
    inferential_results.factor_A.alpha(:,timecurrent)=pcrit;
    inferential_results.factor_A.test_stat(:,timecurrent)=psihat_stat;
    
    for i=1:conAcol;
        inferential_results.factor_A.CI{i,1}(1,timecurrent)=conflow(i);
        inferential_results.factor_A.CI{i,1}(2,timecurrent)=confup(i);
    end

    waitbar(timecurrent/numpnts,h3,sprintf('%12s',[num2str(timecurrent),'/',num2str(numpnts)]))
end

% edit may 8th/15
% clean temporary mapped files
if iscell(tmpfname)
    for i=1:length(tmpfname);
        delete([tmpfname{i}, '.map']);
    end
    
else
    delete([tmpfname, '.map']);
end


%if iscell(tmpfname_CIup)
    for i=1:length(tmpfname_CIup);
        delete([tmpfname_CIup{i}, '.map']);
        delete([tmpfname_CIlow{i}, '.map']);
        delete([tmpfname_diff{i}, '.map']);
    end
    
%else
%    delete([tmpfname_CIup, '.map']);
%    delete([tmpfname_CIlow, '.map']);
%    delete([tmpfname_diff, '.map']);
%end


close(h3)
end


