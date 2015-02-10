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

for i=1:colfile;
    [subrow subcol]=size(condfiles_subs{i});
    subdata=zeros(nboot,numpnts,subrow);
    
    for j=1:subrow;
        tempload=load(condfiles_subs{1,i}{j,:});
        subdata(:,:,j)=tempload.data;
        clear tempload
    end
    
    % create grand avergae surrogates
    alldatacell{1,i}=mean(subdata,3);
end
