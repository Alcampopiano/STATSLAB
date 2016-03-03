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

% preallocate
[rowfile colfile]=size(condfiles_subs);
alldatacell=cell(1,colfile);


if any(strcmp({'chanclust' 'gfa'},STATS.measure));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %datacellind=randi(rowfile,1000,rowfile); % jun4th/15, remove after testing
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
        alldatacell{1,i}=mean(subdata,3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
elseif any(strcmp({'ersp' 'itc'},STATS.measure));
    
    % convert ms to TFs for tfbsline
%     if ~strcmp(STATS.tfbsline,'none');
%         MStoTF(1)=round((STATS.tfbsline(1)/1000-STATS.xmin)/(STATS.xmax-STATS.xmin) * (STATS.numpnts-1))+1;
%         MStoTF(2)=round((STATS.tfbsline(2)/1000-STATS.xmin)/(STATS.xmax-STATS.xmin) * (STATS.numpnts-1))+1;  
%     end
%     TFmin=min(STATS.TF_times);
%     TFmax=max(STATS.TF_times);
%     
%     
%     

% standard baseline
[val ind]=min(abs(STATS.TF_times));
MStoTF(1)=1;
MStoTF(2)=ind;
    
    for i=1:colfile;
        [subrow subcol]=size(condfiles_subs{i});
        
        for j=1:subrow;
            
            % memory map load
            datamap=mapread(condfiles_subs{1,i}{j,:},'dat','datsize',[STATS.freqbins,STATS.timesout,STATS.nboot]);
            
            
            
            if ~strcmp(STATS.tfbsline,'none');
                
                %%% Remove baseline for each subject
                %%%%%%%%%%%%%%%%%%%%%%
                meanrem=zeros(STATS.freqbins,STATS.timesout, STATS.nboot);
                for pg=1:STATS.nboot;
                   
                    meangather=mean(datamap.Data.dat(:,MStoTF(1):MStoTF(2),pg),2);
                    meanrep=repmat(meangather,1,STATS.timesout);
                    meanrem(:,:,pg)=datamap.Data.dat(:,:,pg)-meanrep;
                end
                
                [~,tmpmeanrem1]=fileparts(tempname);
                mapwrite(meanrem,[tmpmeanrem1,'.map'],'datsize',[STATS.freqbins STATS.timesout,STATS.nboot]);
                datamap=mapread([tmpmeanrem1,'.map'],'dat','datsize',[STATS.freqbins,STATS.timesout,STATS.nboot]);
                
                
            end
            
            
            %clear datamap
            %%%%%%%%%%%%%%%%%%%%%%%
            % edit
            % write the condition waveforms for each subject to disk
            condition=mean(datamap.Data.dat,3);
            
            save([condfiles_subs{1,i}{j,:}(1:end-4),'_TF_waves.mat'],'condition');
            clear condition
            
            if j==1;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                if strcmp(STATS.tfbsline,'none');
                    
                    % should be here when no baseline rem
                    dat_avg=datamap;
                    
                else
                    % writing and reading same array to disk to keep cumulative
                    % averaging code consistent when doing baseline rem
                    [~,tmpmeanrem2]=fileparts(tempname);
                    mapwrite(datamap.Data.dat,[tmpmeanrem2,'.map'],'datsize',[STATS.freqbins STATS.timesout,STATS.nboot]);
                    dat_avg=mapread([tmpmeanrem2,'.map'],'dat','datsize',[STATS.freqbins,STATS.timesout,STATS.nboot]);
                    
                end
                % create a temp file name to store large freq surrogates
                %[~,tmpfname]=fileparts(tempname);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            elseif j>1;
                
                if j>2;
                    oldtmpfname=tmpfname;
                end
                % create a temp file name to store large freq surrogates
                [~,tmpfname]=fileparts(tempname);
                %mapwrite((dat_avg.Data.dat+datamap.Data.dat)/2,[tmpfname,'.map'],'datsize',[STATS.freqbins STATS.timesout,STATS.nboot]);
                mapwrite((dat_avg.Data.dat+datamap.Data.dat),[tmpfname,'.map'],'datsize',[STATS.freqbins STATS.timesout,STATS.nboot]);
                dat_avg=mapread([tmpfname,'.map'],'dat','datsize',[STATS.freqbins,STATS.timesout,STATS.nboot]);
                
                try
                    warning off
                    delete([oldtmpfname,'.map']);
                    warning on
                catch
                end
                
                
            end
            clear datamap
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % take average based on subrow subjects in each condition
        [~,tmpfname_mean]=fileparts(tempname);
        mapwrite((dat_avg.Data.dat./subrow),[tmpfname_mean,'.map'],'datsize',[STATS.freqbins STATS.timesout,STATS.nboot]);
        clear dat_avg
        dat_avg=mapread([tmpfname_mean,'.map'],'dat','datsize',[STATS.freqbins,STATS.timesout,STATS.nboot]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % try to delete previous mapped files;
        warning off
        delete(['groupboots_',STATS.savestring, '_', STATS.condnames{i},'.map']);
        warning on
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~strcmp(STATS.tfbsline,'none');
            % removal of baseline, edit
            meanrem=zeros(STATS.freqbins,STATS.timesout, STATS.nboot);
            for pg=1:STATS.nboot;
                
                % need flexible inputs
                meangather=mean(dat_avg.Data.dat(:,MStoTF(1):MStoTF(2),pg),2); % baseline period
                meanrep=repmat(meangather,1,STATS.timesout);
                meanrem(:,:,pg)=dat_avg.Data.dat(:,:,pg)-meanrep;
            end
            
            % save full surrogate arrays X condition
            mapwrite(meanrem,['groupboots_',STATS.savestring, '_', STATS.condnames{i},'.map'],'datsize',[STATS.freqbins STATS.timesout STATS.nboot]);
            
            
            % save mean TF waveforms
            condition=mean(meanrem,3);
            
        else
            
            % this should be here instead of meanrem when no baseline removal option is chosen
            mapwrite(dat_avg.Data.dat,['groupboots_',STATS.savestring, '_', STATS.condnames{i},'.map'],'datsize',[STATS.freqbins STATS.timesout STATS.nboot]);
            
            % this should be here instead of meanrem when no baseline removal option is chosen
            condition=mean(dat_avg.Data.dat,3); %%%% keep
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        save(['group_TFwaves_',STATS.savestring, '_', STATS.condnames{i},'.mat'],'condition');
        clear condition
        %mapwrite(mean(dat_avg.Data.dat,3),['group_TFwaves_',STATS.savestring, '_', STATS.condnames{i},'.map'],'datsize',[STATS.freqbins STATS.timesout]);
        
        % fill up datacell with the full TF surrogates to be used in statistics.
        alldatacell{1,i}=mapread(['groupboots_',STATS.savestring, '_', STATS.condnames{i},'.map'], 'dat','datsize',[STATS.freqbins,STATS.timesout,STATS.nboot]);
        
        
        % try to delete previous mapped files;
        try
            warning off
            delete([tmpmeanrem1, '.map']);
            delete([tmpmeanrem2, '.map']);
            delete([tmpfname,'.map']);
            delete([tmpfname_mean,'.map']);
            warning on
        catch
        end
        clear dat_avg
        
    end
    
end





