function [condfiles results] = pbsubject1waytf(STATS, condfiles, numconds, nboot, jlvls, alpha, condnames, varargin)
tic
% bla bla bla

nargs = length(varargin);
if round(nargs/2)~=nargs/2
    error('need propertyName/propertyValue pairs for optional inputs')
end

% Set default contrast coefficients for 1-way
[conA] = con1way(jlvls);

% put defaults into a structure;
% options=struct('conA',conA);

% edited may1st/15
% set default plot options
options.conA=conA;
options.FWE='Rom';

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

% extract from options structure
conA=options.conA;

% used to create proper sizes in results structure
[~, conAcol]=size(conA);

% if condfiles is empty, bring up browser, if not then you likely created
% them in some parent function, perhaps to run single-subject stats in a
% 'bw' design (multiple one-way tests...)

if isempty(condfiles);
    
    % get rid of condfiles as an empty char string ''
    clear condfiles
    
    % load all file names subs X conditions
    for i=1:numconds
        tempfname=uigetfile('*.map',['Select all bootstrapped files in the ', condnames{i}, ' condition'], 'MultiSelect','on');
        if ~iscell(tempfname);
            tempfname={tempfname};
            condfiles(:,i)=tempfname;
        else
            condfiles(:,i)=tempfname;
        end
    end
    
    % save so one can load without gui
    save(['condfiles_SubjectStatistics_',STATS.measure,'_',STATS.savestring,'.mat'],'condfiles');
    
else
    
    % load a file name that was given that contains the filenames X condition cell array
    condfiles_data=load(condfiles);
    condfields=fieldnames(condfiles_data);
    condfiles=condfiles_data.(condfields{1});
    
    %%%% there is probably a better way of doing this, but I was tired at
    %%%% that moment
    % fix cell within a cell
%     condfiles=cell(1,numconds);
%     for i=1:numconds
%         [rowcondcell colcondcell]=size(condfiles_cellcell{i});
%         for j=1:rowcondcell
%             condfiles{j,i}=condfiles_cellcell{i}{j};
%         end
%     end
    
end


%preallocate sizes
[rowconds colconds]=size(condfiles);
datacell=cell(1,colconds);

for i=1:rowconds;
    field_name{i,1}=['subject_', num2str(i)];
end

% band fields
for i=1:STATS.freqbins;
    band_fields{i,1}=['band_', strrep(num2str(STATS.TF_freqs(i)),'.','_')];
end

% for testcurrent=1:rowconds;
%     
%     results.(field_name{testcurrent})=struct('factor_A',{[]});
%     results.(field_name{testcurrent}).factor_A=struct('contrasts',{conA},'pval',{zeros(conAcol,numpnts)},'alpha',{zeros(conAcol,numpnts)},'test_stat',{zeros(conAcol,numpnts)},'CI',{cell(conAcol,1)}, 'FWE', options.FWE);
%     
%     
%     for i=1:conAcol;
%         results.(field_name{testcurrent}).factor_A.CI{i,1}=zeros(2,numpnts);
%     end
%     
% end


%%%%%%%%% finding near val for bin reduction
% frange=[5.2 7.5];
% 
% freqs=[3 3.5 5 5.5 7 7.2 8 8.9 10];
% [lowval lowind] = min(abs(freqs-frange(1)))
% low_near = freqs(lowind) 
% 
% [hival hiind] = min(abs(freqs-frange(2)))
% hi_near = freqs(hiind) 

h1 = waitbar(0,'1','Name','subject progress','Position',[1100 549 550 40]);
childh1 = get(h1, 'Children');
set(childh1, 'Position',[5 10 538 15]);

h2 = waitbar(0,'1','Name','statistics on frequency bin','Position',[1100 486 550 40]);
childh2 = get(h2, 'Children');
set(childh2, 'Position',[5 10 538 15]);

for filecurrent=1:rowconds;
    
    for condcurrent=1:colconds;
        %conds=load(condfiles{filecurrent,condcurrent});
        
        % memory map load
        datamap=mapread(condfiles{filecurrent,condcurrent},'dat','datsize',[STATS.freqbins,STATS.timesout,STATS.nboot]);    
        
%         if ~strcmp(STATS.tfbsline,'none');
%             
%             %%% Remove baseline for each subject
%             %%%%%%%%%%%%%%%%%%%%%%
%             meanrem=zeros(STATS.freqbins,STATS.timesout, STATS.nboot);
%             for pg=1:STATS.nboot;
%                 
%                 % need flexible inputs
%                 meangather=mean(datamap.Data.dat(:,STATS.tfbsline(1):STATS.tfbsline(2),pg),2); % baseline period
%                 meanrep=repmat(meangather,1,STATS.timesout);
%                 meanrem(:,:,pg)=datamap.Data.dat(:,:,pg)-meanrep;
%             end
%             
%             [~,tmpmeanrem{condcurrent}]=fileparts(tempname);
%             mapwrite(meanrem,[tmpmeanrem{condcurrent},'.map'],'datsize',[STATS.freqbins STATS.timesout,STATS.nboot]);
%             datamap=mapread([tmpmeanrem{condcurrent},'.map'],'dat','datsize',[STATS.freqbins,STATS.timesout,STATS.nboot]);
%             
%         end
        datacell{1,condcurrent}=datamap;
        clear datamap
        
        % write the condition waveforms for each subject to disk
        condition=mean(datacell{1,condcurrent}.Data.dat,3);
        save([condfiles{filecurrent,condcurrent}(1:end-4),'_TF_waves.mat'],'condition');
        clear condition
        
    end
    
    %arrange the data for the calculations
    %[rowcell ~]=size(datacell{1,1});  
    
    rowcell=STATS.nboot;
    
    for bandind=1:STATS.freqbins;
        
        % stats at each timepoint
        %for timecurrent=1:numpnts;
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
            results.(field_name{filecurrent}).(band_fields{bandind}).factor_A.contrasts=conA;
            results.(field_name{filecurrent}).(band_fields{bandind}).factor_A.pval(:,timecurrent)=pvalgen;
            results.(field_name{filecurrent}).(band_fields{bandind}).factor_A.alpha(:,timecurrent)=pcrit;
            results.(field_name{filecurrent}).(band_fields{bandind}).factor_A.test_stat(:,timecurrent)=psihat_stat;
          
            
            for i=1:conAcol;
                results.(field_name{filecurrent}).(band_fields{bandind}).factor_A.CI{i,1}(1,timecurrent)=conflow(i);
                results.(field_name{filecurrent}).(band_fields{bandind}).factor_A.CI{i,1}(2,timecurrent)=confup(i);
            end
            
            %waitbar(timecurrent/numpnts,h2,sprintf('%12s',[num2str(timecurrent),'/',num2str(numpnts)]))
        end
        
        waitbar(bandind/STATS.freqbins,h2,sprintf('%12s',[num2str(bandind),'/',num2str(STATS.freqbins)]))
    end
    
    waitbar(filecurrent/rowconds,h1,sprintf('%12s',[num2str(filecurrent),'/',num2str(rowconds)]))
end
close(h1,h2);
toc






end

