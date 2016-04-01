function [STATS]=topobuild(STATS, locsfile, figtype)

% coming from subject figure, or group figure
if strcmp(figtype, 'group')
    
    % determine the set files that were used in the statistics
    topocell=STATS.bootfiles;
    bootfiles=STATS.bootfiles;
    
elseif strcmp(figtype, 'subject') && ~strcmp(STATS.design, 'bw')
    % determine the set files that were used in the statistics
    %topocell=STATS.subject_bootfiles;
    %bootfiles=STATS.subject_bootfiles;
    
    topocell=cell(1,STATS.numconds);
    for k=1:STATS.levels(1);
        topocell{k}(:,1)=STATS.subject_bootfiles(:,k)';
        
    end
    bootfiles=topocell;
    
elseif strcmp(figtype, 'subject') && strcmp(STATS.design, 'bw')
    % determine the set files that were used in the statistics
    
    % unnest subject_bootfiles
    topocell=cell(1,STATS.numconds);
    q=1;
    for j=1:STATS.levels(1);
        for k=1:STATS.levels(2);
            topocell{q}=STATS.subject_bootfiles{j}(:,k);
            q=q+1;
        end  
    end
    bootfiles=topocell;
end

STATS.subtopofiles=topocell;
STATS.grouptopofiles=cell(1,length(STATS.condnames));

for i=1:STATS.numconds;
    
    [rowbt colbt]=size(bootfiles{i});
    %[rowbt colbt]=size(STATS.bootfiles{i});
    for k=1:rowbt;
        
        ind=strfind(bootfiles{i}{k},STATS.datatype);
        topocell{i}{k}=[topocell{i}{k}(1:ind-2), '.set'];
    end
    
end

res='';
for i=1:STATS.numconds;
    
    [row col]=size(topocell{i});
    
    for j=1:row;
        
        EEG=pop_loadset('filename', topocell{i}{j});
        EEG = eeg_checkset( EEG );
        
        % is this ICA data
        if isfield(STATS, 'ICretain') && any(strfind(STATS.datatype, 'ica'))
            disp('assuming this is an analysis on IC data');
            
            % get orignal IC inds
            icorig=[1:min(size(EEG.icawinv))];

            % just the numbers
            ICinds=STATS.ICretain(:,2:2:end);
            
            % just names
            ICfnames=STATS.ICretain(:,1:2:end);
            
            % find correct IC inds
            ind=find(strcmp(topocell{i}{j}, ICfnames(:,i)));
            compswant=ICinds(ind,i);
            compswant=compswant{1};

            % is it residual
            if ~any(strfind(bootfiles{1}{1},[STATS.datatype, '_', 'RESextracted']));
                
                
                % remove inds
                icorig(compswant)=[];
                EEG = pop_subcomp( EEG, icorig, 0);
                EEG = eeg_checkset(EEG);
            else
                res='RES';
                % remove inds
                EEG = pop_subcomp( EEG, compswant, 0);
                EEG = eeg_checkset(EEG);
            end
            
        end
        
        % do all the nasty stuff
        EEG.data=mean(EEG.data,3); % making interpolating much faster
        EEG.trials=1; % so interpolation will work
        EEG=interpmont(EEG,locsfile);
        EEG = eeg_checkset( EEG );
        EEG=pop_chanedit(EEG, 'load',{locsfile 'filetype' 'autodetect'});
        EEG = eeg_checkset( EEG );
        
        % take EEG fields you need for group topos
        if j==1 && i==1;
            tmpEEG.xmin=min(EEG.times/1000);
            tmpEEG.xmax=max(EEG.times/1000);
            tmpEEG.chanlocs=EEG.chanlocs;
            tmpEEG.pnts=EEG.pnts;
            tmpEEG.nbchan=EEG.nbchan;
            tmpEEG.trials=1;
        end
        
        if j==1;
            
            EEGsubs=EEG.data;
            
            % save sub
            save([topocell{i}{j}(1:end-4), '_subtopo_', res, STATS.savestring, '.mat'], 'EEG');
            
            %populate array for stats struct
            STATS.subtopofiles{i}{j}=[topocell{i}{j}(1:end-4), '_subtopo_', res, STATS.savestring, '.mat'];
            clear EEG
            
        else
            % accumulate until averaging
            EEGsubs=EEGsubs+EEG.data;
            
            % save sub
            save([topocell{i}{j}(1:end-4), '_subtopo_', res, STATS.savestring, '.mat'], 'EEG');
            
            %populate array for stats struct
            STATS.subtopofiles{i}{j}=[topocell{i}{j}(1:end-4), '_subtopo_', res, STATS.savestring, '.mat'];
            clear EEG
        end
        
    end
    
    % mean
    EEGsubs=EEGsubs/row;
    tmpEEG.data=EEGsubs;
    save(['grouptopo_', res, STATS.savestring, '_',STATS.condnames{i}, '.mat'], 'tmpEEG');
    
    % populate stats
    STATS.grouptopofiles{i}=['grouptopo_', res, STATS.savestring, '_',STATS.condnames{i}, '.mat'];
    
end




