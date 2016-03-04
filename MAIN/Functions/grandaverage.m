function [alldatacell] = grandaverage(STATS,nboot,numpnts,condfiles_subs)

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

disp('deleting temp map files first');
delete('tp*.map')

% preallocate
[rowfile colfile]=size(condfiles_subs);
alldatacell=cell(1,colfile);


if any(strcmp({'chanclust' 'gfa'},STATS.measure));
    
    for i=1:colfile;
        [subrow subcol]=size(condfiles_subs{i});
        subdata=zeros(nboot,numpnts,subrow);
        
        for j=1:subrow;
            tempload=load(condfiles_subs{1,i}{j,:});
            subdata(:,:,j)=tempload.data;
            clear tempload
        end
        
        alldatacell{1,i}=mean(subdata,3);
        
    end
    
elseif any(strcmp({'ersp' 'itc'},STATS.measure));
    
    for i=1:colfile;
        [subrow subcol]=size(condfiles_subs{i});
        
        for j=1:subrow;
            
            % memory map load
            datamap=mapread(condfiles_subs{1,i}{j,:},'dat','datsize',[STATS.freqbins,STATS.timesout,STATS.nboot]);
            
            % write the condition waveforms for each subject to disk
            condition=mean(datamap.Data.dat,3);
            
            save([condfiles_subs{1,i}{j,:}(1:end-4),'_TF_waves.mat'],'condition');
            clear condition
            
            if j==1;
                dat_avg=datamap;
            else
                
                % create a temp file name to store large freq surrogates
                [~,tmpfname{j}]=fileparts(tempname);
                mapwrite((dat_avg.Data.dat+datamap.Data.dat),[tmpfname{j},'.map'],'datsize',[STATS.freqbins STATS.timesout,STATS.nboot]);
                dat_avg=mapread([tmpfname{j},'.map'],'dat','datsize',[STATS.freqbins,STATS.timesout,STATS.nboot]);
                              
                if j>2;
                    delete([tmpfname{j-1},'.map']);
                end
                
            end
            clear datamap
        end
       
        % take average based on subrow subjects in each condition
        [~,tmpfname_mean]=fileparts(tempname);
        mapwrite((dat_avg.Data.dat./subrow),[tmpfname_mean,'.map'],'datsize',[STATS.freqbins STATS.timesout,STATS.nboot]);
        clear dat_avg
        dat_avg=mapread([tmpfname_mean,'.map'],'dat','datsize',[STATS.freqbins,STATS.timesout,STATS.nboot]);
        
        % try to delete previous mapped files;
        warning off
        delete(['groupboots_',STATS.savestring, '_', STATS.condnames{i},'.map']);
        warning on
        
        % this should be here instead of meanrem when no baseline removal option is chosen
        mapwrite(dat_avg.Data.dat,['groupboots_',STATS.savestring, '_', STATS.condnames{i},'.map'],'datsize',[STATS.freqbins STATS.timesout STATS.nboot]);
        
        % this should be here instead of meanrem when no baseline removal option is chosen
        condition=mean(dat_avg.Data.dat,3);

        save(['group_TFwaves_',STATS.savestring, '_', STATS.condnames{i},'.mat'],'condition');
        clear condition
        
        % fill up datacell with the full TF surrogates to be used in statistics.
        alldatacell{1,i}=mapread(['groupboots_',STATS.savestring, '_', STATS.condnames{i},'.map'], 'dat','datsize',[STATS.freqbins,STATS.timesout,STATS.nboot]);
        
        clear dat_avg
        
    end
    
    % clean up
    delete('tp*.map')
end





