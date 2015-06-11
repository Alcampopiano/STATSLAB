function [alldatacell] = grandaverage(nboot,numpnts,condfiles_subs)

% creates grand averages that can be used for robust group stats
% like in Rousselet 2008, and Desjardins 2013 (percentile bootstrap tests).

% It will stack subjects along the 3rd dimension (e.g., 1000 surrogates X time X subjects)
% and then average across subjects so the resulting array is for example 1000 X TFs.

%{
    Inputs:
    nboot = number of bootstrapps that are in the incoming subject files
    numpnts = TFs
    numsubs = Number of subjects
    condfiles_subs=cell array of strings (indicating a subjects surrogates X TF) for each condition separated by cell
%}

% preallocate
[rowfile colfile]=size(condfiles_subs);
alldatacell=cell(1,colfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datacellind=randi(rowfile,1000,rowfile); % jun4th/15, remove after testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:colfile;
    [subrow subcol]=size(condfiles_subs{i});
    subdata=zeros(nboot,numpnts,subrow);
    
    for j=1:subrow;
        tempload=load(condfiles_subs{1,i}{j,:});
        subdata(:,:,j)=tempload.data;
        clear tempload
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TESTING resampling procedure, remove or incorporate this after testing jun4th/15
    % resample in a different way from Datacell, the result of
    % GroupFigure_sample & GroupStatistics_sample will now be a resampling, NOT
    % simply including all subjects as this function should normally do
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for q=1:1000;
        alldatacell{1,i}(q,:)=trimmean(subdata(q,:,datacellind(q,:)),40,3);
    end
    % this line SHOULD be here but was commented out due to testing on
    % jun4th/15
    % create grand avergae surrogates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %alldatacell{1,i}=mean(subdata,3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end









