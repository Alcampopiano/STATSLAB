function [condfiles results] = pbsubject1way(STATS, condfiles, numconds, numpnts, nboot, jlvls, alpha, condnames, varargin)
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
%options.conB=[]; % here as a place holder for bw designs (which are done as seperate 1-ways)
options.conA=conA;
options.FWE='Rom';

% for bw designs
options.jlabels={};
options.klabels={};

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

% move to STATS struct
STATS.jlabels=options.jlabels;
STATS.klabels=options.klabels;

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
        tempfname=uigetfile('*.mat',['Select all bootstrapped files in the ', condnames{i}, ' condition'], 'MultiSelect','on');
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
%         condfiles{j,i}=condfiles_cellcell{i}{j};
%         end  
%     end
    
end


%preallocate sizes
[rowconds colconds]=size(condfiles);
datacell=cell(1,colconds);

for i=1:rowconds;
    field_name{i,1}=['subject_', num2str(i)];
end

for testcurrent=1:rowconds;
    
    results.(field_name{testcurrent})=struct('factor_A',{[]});
    results.(field_name{testcurrent}).factor_A=struct('contrasts',{conA},'pval',{zeros(conAcol,numpnts)},'alpha',{zeros(conAcol,numpnts)},'test_stat',{zeros(conAcol,numpnts)},'CI',{cell(conAcol,1)}, 'FWE', options.FWE);
    
    
    for i=1:conAcol;
        results.(field_name{testcurrent}).factor_A.CI{i,1}=zeros(2,numpnts);
    end
    
end


h1 = waitbar(0,'1','Name','subject progress','Position',[1100 549 550 40]);
childh1 = get(h1, 'Children');
set(childh1, 'Position',[5 10 538 15]);

h2 = waitbar(0,'1','Name','stats across time','Position',[1100 486 550 40]);
childh2 = get(h2, 'Children');
set(childh2, 'Position',[5 10 538 15]);

for filecurrent=1:rowconds;
    
    for condcurrent=1:colconds;
        conds=load(condfiles{filecurrent,condcurrent});
        %conds=load(condfiles{1,condcurrent}{filecurrent});
        datacell{1,condcurrent}=conds.data;
    end
    
    %arrange the data for the calculations
    [rowcell ~]=size(datacell{1,1});
    
    
    % stats at each timepoint
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
        results.(field_name{filecurrent}).factor_A.pval(:,timecurrent)=pvalgen;
        results.(field_name{filecurrent}).factor_A.alpha(:,timecurrent)=pcrit;
        results.(field_name{filecurrent}).factor_A.test_stat(:,timecurrent)=psihat_stat;
        
        
        for i=1:conAcol;
            results.(field_name{filecurrent}).factor_A.CI{i,1}(1,timecurrent)=conflow(i);
            results.(field_name{filecurrent}).factor_A.CI{i,1}(2,timecurrent)=confup(i);
        end
        
        waitbar(timecurrent/numpnts,h2,sprintf('%12s',[num2str(timecurrent),'/',num2str(numpnts)]))
    end
    
    waitbar(filecurrent/rowconds,h1,sprintf('%12s',[num2str(filecurrent),'/',num2str(rowconds)]))
end
close(h1,h2);
toc
end

