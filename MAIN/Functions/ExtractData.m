%                                                                                          
% This function extracts various dependent measures (e.g., ICs, gfa, microvolts) from EEGLAB's .set files. The output files are then used in the subsequent module ResampleData.m
%                                                                                          
% 
% Inputs:
% 
% 
% ***condfiles***
% 
% Leave emtpy and MATLAB will bring up an interface for you to load the appropriate SET files. After this is done the first time, a file is saved (condfiles*.mat). For subsequent calls type this filename in to avoid having to manually choose files again. ***end***
% 
% ***condnames***
% 
% Type in your condition labels on separate lines. For example,
% 
% face
% house 
% object ***end***
% 
% ***levels***
% 
% Levels for each factor in your design. Leave a space between factors.
% 
% For example, 
% 
% 2 3
% 
% denotes a 2x3 design.
% 
% 2 
% 
% denotes a two-condition design. Three-way designs not allowed. ***end***
% 
% ***design***
% 
% A string indicating your design. Choose from the drop down menu.
% 
% w
% 
% denotes a within-subjects design with two conditions. 
% 
% ww 
% 
% denotes a two-factor, fully within-subjects design. 
% 
% bw 
% denotes a mixed design  ***end***
% 
% ***savestring***
% 
% A keyword that will be appended to important files as you work through STATSLAB modules. It should identify your study and/or analysis. The STATS structure will have savestring appended to it. 
% 
% For example,
% 
% Fcz_ICA_analysis ***end***
% 
% ***varargin***
% 
% Options are specified in pairs (key -> val)
%
% Measure ->
% 		
% 	icamax		project IC to channel with max weight
% 		
% 	icagfa		project IC to scalp and measure scalp GFA
% 		
% 	icaitc		project IC to selected channel and calculate inter-trial coherence
% 		
% 	icaersp		project IC to selected channel and calculate event-related spectral perturbation
% 		
% 	icascalp	project IC to selected channel and measure microvolt
% 		
% 	scalpitc	inter-trial coherence for selected channel
% 		
% 	scalpersp	event-related spectral perturbation for selected channel
% 		
% 	scalpgfa	global field power for scalp
% 		
% 	scalpchan	measure microvolts for selected scalp channel
% 		
% ICs ->		
% 	persubject		bring up GUI to enter IC indexes for each subject and condition
% 
% 	your_IC_file.mat	name of file that hold IC indexes (you can make file using the “persubject” option first
% 		
% 		
% Chans ->		
% 
% 	persubject	bring up montage GUI to select channels for each subject and condition
% 
% 	chanlabels 	The channel labels for electrodes you are analyzing (applied to all subjects)
% 
% 		
% tfcycles,freqs,nfreqs -> 		
% 		
%       To be used with the itc and ersp options. See EEGLAB's newtimef.m for information on these key & val options	
% 		
% For example,
% 		
% measure		
% scalpchan		
% chans		
% Fcz Cz	
% 
% will extract standard scalp microvolts for Fcz and Cz (channels are averaged together)	
% 		
% For example,		
% 		
% measure		
% icascalp		
% ICs		
% persubject		
% chans		
% persubject		
% 
% A GUI will appear to select IC indexes. Those ICs will then be projected back to the channels selected in the montage GUI. 
% 
% Using ExtractData at the commandline:
% 
% ExtractData({'face' 'house'}, [], 2, 'w', 'my_Oz_analysis', 'measure', 'scalpchan', 'chans', 'persubject');
% ExtractData({'old_cong' 'old_incong' 'young_cong' 'young_incong'}, [], [2 2], 'bw', 'my_Oz_ICAanalysis', 'measure', 'icascalp', 'ICs', 'my_IC_file.mat', 'chans', 'Oz');
% 
% ***end***
%
%
% Copyright (C) <2015>  <Allan Campopiano>
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [STATS]=ExtractData(condnames,condfiles,levels,design,savestring,varargin)

% set history
[hist_str]=statslab_history(condnames,condfiles,levels,design,savestring,varargin);
STATS.history.ExtractData=hist_str;


%%%%% options %%%%%
options.measure=[];
options.chans='persubject'; % key word default, can be chan label -> 'Cz' or -> {'Fcz' 'Cz'}
options.ICs=NaN;
options.tfcycles=[3 .5]; % spectral opts
options.freqs=[3 30]; % spectral opts
options.nfreqs=27;
%options.timesout=600; % spectral opts
%options.tfbsline='none'; % TF baseline default


% get field names
optionnames = fieldnames(options);

nargs=length(varargin);
if round(nargs/2)~=nargs/2
    error('need NAME/VALUE pairs for optional inputs')
end

% handle varargin
for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
    inpName = pair{1};
    
    if any(strcmp(inpName,optionnames))
        
        % overwrite default options
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name',inpName)
    end
end

% turn char to cell array if chan was specified as a char
if ischar(options.chans) && ~strcmp(options.chans, 'persubject')
    
    % then its a char string and needs to be a cell
    options.chans={options.chans};
end


% there should be a label for each conditon.
numconds=length(condnames);

% if user selected file names once, they do not have to go through that again, provided they give file name for the condfiles argument.
if isempty(condfiles);
    
    %load all file names subs X conditions
    for i=1:numconds
        [tempfname temppath]=uigetfile('*.set',['Select all .set files in the ', condnames{i}, ' condition'], 'MultiSelect','on');
        if ~iscell(tempfname)
            tempfname={tempfname};
            %warning('MATLAB:questionableactivity','You only selected one subject for this condition. Is that what you want?');
        end
        
        % see if there are any zeros, if so give error msg. This would happen if the user clicked cancel in the uigetfile pane.
        % if this error occurs a useless array of zeros is not saved later down the line, which is good.
        if ~any(tempfname{1});
            error(['You did not select any files for some condition. ' ...
                'Once you select the condition files properly the first time, they are saved and you do not need to select them again. ' ...
                'Just provide condfiles.mat as an input argument on subsequent calls'])
        else
            condfiles_subs{1,i}(:,1)=tempfname;
            pathtofiles{1,i}(1,:)=temppath;
        end
    end
    
    % so, if I didnt find any zeros, save the array. Assuming that care was
    % put into its creation by the user.
    condfiles{1}=condfiles_subs;
    condfiles{2}=pathtofiles;
    save(['condfiles_ExtractData_',options.measure,'_',savestring,'.mat'],'condfiles');
    
else
    load(condfiles);
    condfiles_subs=condfiles{1};
    pathtofiles=condfiles{2};
end


% parse j and k lvls, or just jlvls depending on length of levels.
% parse j and k levels
if length(levels)==2;
    jlvls=levels(1);
    klvls=levels(2);
    
elseif length(levels)==1;
    jlvls=levels(1);
    
else
    error('levels can only have length of 1 (for 1-way design) or 2 (for 2-way design)');
end


%%
% GUI options for ICs
if strcmp(options.ICs, 'persubject');
    
    % GUI IC picker
    [okayhit, ICretain]=ICpick(condfiles_subs,numconds);
    
    if isempty(okayhit)
        STATS=[]; % to exit cleanly
        return
    end
    
elseif ~any(isnan(options.ICs)) %&& strcmp(options.ICs(end-3:end), '.mat')
    
    % load a file
    tmp=load(options.ICs);
    ICretain=tmp.ICCHOICES;
    
    % turn IC selections into numbers
    j=2;
    for i=1:numconds;
        [rowfname colfname]=size(condfiles_subs{i});
        ICretain(1:rowfname,j)=cellfun(@str2num,ICretain(1:rowfname,j),'UniformOutput',0)';
        j=j+2;
    end
    
end

%%
% GUI for channel selections

persubject=0;
if strcmp(options.chans,'persubject')
    persubject=1;
    
    % GUI channel picker not written yet
    [cancel_hit chanfile]=chanpick_topo(condfiles_subs, pathtofiles, numconds);
    
    if cancel_hit
        STATS=[]; % to exit cleanly
        return
    end
    
elseif ~isempty(strfind(options.chans{1},'.mat'))
    
    % load a file
    tmp=load(options.chans{1});
    chanfile=tmp.data.chanarray;
end
%% cases

switch options.measure
    
    case {'icascalp', 'icaersp', 'icaitc'}
        disp(' ***** projecting selected ICs to scalp and extracting channels ***** ')
        
        for i=1:numconds
            [rowcond colcond]=size(condfiles_subs{i});
            
            for j=1:rowcond;
                
                % load file
                EEG = pop_loadset('filename',condfiles_subs{i}{j},'filepath',pathtofiles{i});
                EEG = eeg_checkset(EEG);
                
                % get orignal IC inds
                icorig=[1:min(size(EEG.icawinv))];
                
                compswant=ICretain(j,i+i);
                compswant=compswant{1};
                
                icorig(compswant)=[];
                
                % pick ICs you want to retain
                EEGretain = pop_subcomp( EEG, icorig, 0);
                EEGretain = eeg_checkset(EEGretain);
                
                if persubject==0;
                    
                    % when some option to select a channel is given
                    try
                        % scroll through chans the user wants and collect relavent indices
                        for ii=1:length(options.chans)
                            chanind(ii)=find(strcmp({EEG.chanlocs.labels},options.chans{ii}));
                        end
                        
                    catch
                        error(['Channel ' options.chans{ii}, 'does not exist for subject ', condfiles_subs{i}{j}]);
                    end
                    
                elseif persubject==1;
                    
                    % when some option to select a channel is given
                    try
                        % scroll through chans the user wants and collect relavent indices
                        for ii=1:length(chanfile{j,i+i})
                            chanind(ii)=find(strcmp({EEG.chanlocs.labels},chanfile{j,i+i}{ii}));
                        end
                        
                    catch
                        error(['Channel ' chanfile{j,i+i}{ii}, 'does not exist for subject ', condfiles_subs{i}{j}]);
                    end
                    
                end
                
                data=EEGretain.data(chanind,:,:);
                
                % save it with original filename but get rid of original
                % extention (hence the 1:end-4)
                save([condfiles_subs{i}{j}(1:end-4),'_',options.measure,'_extracted.mat'],'data');
                clear data
                
                % get residual
                EEGres = pop_subcomp( EEG, compswant, 0);
                EEGres = eeg_checkset(EEGres);
                
                % put EEG.data into the variable data
                data=EEGres.data(chanind,:,:);
                
                % save it with original filename but get rid of original
                % extention (hence the 1:end-4)
                save([condfiles_subs{i}{j}(1:end-4),'_',options.measure,'_RESextracted.mat'],'data');
                clear data
                
            end
            
        end
        
        xtimes=EEG.times;
        disp(' ***** finished projecting selected ICs to specified scalp channels ***** ')
        
    case 'icagfa'
        disp(' ***** projecting selected ICs to scalp and extracting data array for later GFA calculation ***** ')
        
        for i=1:numconds
            [rowcond colcond]=size(condfiles_subs{i});
            
            for j=1:rowcond;
                
                % load file
                EEG = pop_loadset('filename',condfiles_subs{i}{j},'filepath',pathtofiles{i});
                EEG = eeg_checkset(EEG);
                
                % get orignal IC inds
                icorig=[1:min(size(EEG.icawinv))];
                
                compswant=ICretain(j,i+i);
                
                icorig(compswant)=[];
                
                % pick ICs you want to retain
                EEGretain = pop_subcomp( EEG, icorig, 0);
                EEGretain = eeg_checkset(EEGretain);
                
                data=EEGretain.data;
                
                % save it with original filename but get rid of original
                % extention (hence the 1:end-4)
                save([condfiles_subs{i}{j}(1:end-4),'_',options.measure,'_extracted.mat'],'data');
                clear data
                
                % get residual
                EEGres = pop_subcomp( EEG, compswant, 0);
                EEGres = eeg_checkset(EEGres);
                
                % put EEG.data into the variable data
                data=EEGres.data;
                
                % save it with original filename but get rid of original
                % extention (hence the 1:end-4)
                save([condfiles_subs{i}{j}(1:end-4),'_',options.measure,'_RESextracted.mat'],'data');
                clear data
                
            end
            
        end
        
        xtimes=EEG.times;
        disp(' ***** finished projecting selected ICs to scalp for later GFA calculations ***** ')
        
    case 'scalpgfa'
        disp(' ***** extracting the data array for later GFA calculations ***** ')
        
        for k=1:numconds;
            
            [rowcond colcond]=size(condfiles_subs{k});
            
            for s=1:rowcond; % scroll through subjects
                
                % load file
                EEG = pop_loadset('filename',condfiles_subs{k}{s},'filepath',pathtofiles{k});
                EEG = eeg_checkset(EEG);
                
                % put EEG.data into the variable data
                data=EEG.data;
                
                % save it with original filename but get rid of original
                % extention (hence the 1:end-4)
                save([condfiles_subs{k}{s}(1:end-4),'_',options.measure,'_extracted.mat'],'data');
                clear data
                
            end
            
        end
        xtimes=EEG.times;
        disp(' ***** finished extracting the data array for later GFA calculations *****')
        
        % likely a temporary case option for sysc14, unless
        % I can keep everything internally consistent without too much trouble
%     case 'bipolar'
%         disp(' ***** extracting the data array for bipolar ***** ')
%         
%         %for k=1:numconds; % always one
%         
%         [rowcond colcond]=size(condfiles_subs{1});
%         
%         miscinfo=cell(rowcond+1,10);
%         miscinfo(1,:)={'sub', 'ind', 'chan', 'val', 'ntrials', 'sub', 'ind', 'chan', 'val', 'ntrials'};
%         miscinfo(2:end,1)=condfiles_subs{1}(:);
%         miscinfo(2:end,6)=condfiles_subs{2}(:);
%         
%         for s=1:rowcond; % scroll through subjects
%             
%             % load file from cond 1
%             EEG = pop_loadset('filename',condfiles_subs{1}{s},'filepath',pathtofiles{1});
%             EEG = eeg_checkset(EEG);
%             miscinfo{s+1,5}=EEG.trials;
%             
%             %create cond 1 chan erps
%             erps_cond1=mean(EEG.data,3);
%             
%             % load file from cond 2
%             EEG = pop_loadset('filename',condfiles_subs{2}{s},'filepath',pathtofiles{1});
%             EEG = eeg_checkset(EEG);
%             miscinfo{s+1,10}=EEG.trials;
%             
%             %create cond 2 chan erps
%             erps_cond2=mean(EEG.data,3);
%             
%             % create the difference ERP waves for each channel
%             erp_diff=erps_cond1-erps_cond2;
%             
%             clear erps_cond1 erps_cond2
%             
%             % copy EEG.data
%             tmpEEG=EEG;
%             tmpEEG.data=erp_diff;
%             
%             % set trials ==1 to trick pop_timtopo into plotting
%             tmpEEG.trials=1;
%             h=figure; pop_timtopo(tmpEEG, [-200  500], [NaN], 'ERP data and scalp maps of left_check');
%             
%             timewin=input('type in the window(ms) to calculate min and max vals\n');
%             close(h);
%             
%             if isempty(timewin)
%                 
%                 if strcmp(condnames{1},'face')
%                     timewin=[150 180]; % for face house
%                     
%                 elseif strcmp(condnames{1},'Left')
%                     timewin=[80 120]; % for left right checker
%                 end
%                 
%             end
%             
%             MStoTF_min=round((timewin(1)/1000-EEG.xmin)/(EEG.xmax-EEG.xmin) * (EEG.pnts-1))+1;
%             MStoTF_max=round((timewin(2)/1000-EEG.xmin)/(EEG.xmax-EEG.xmin) * (EEG.pnts-1))+1;
%             
%             % for each channel get max value within a window
%             [max_val_time max_ind_time]=max(tmpEEG.data(:,MStoTF_min:MStoTF_max,:),[],2);
%             
%             % then, get the index for the max channel
%             [max_val_chan max_ind_chan]=max(max_val_time);
%             
%             % for each channel get min value within a window
%             [min_val_time min_ind_time]=min(tmpEEG.data(:,MStoTF_min:MStoTF_max,:),[],2);
%             
%             % then, get the index for the min channel
%             [min_val_chan min_ind_chan]=min(min_val_time);
%             
%             clear tmpEEG
%             
%             % now we have max and min channel index, steal it from original
%             % data, cond 2 is already loaded so take it from there.
%             datamax=EEG.data(max_ind_chan,:,:);
%             datamin=EEG.data(min_ind_chan,:,:);
%             
%             % create bipolar
%             data=datamax-datamin;
%             miscinfo{s+1,7}=max_ind_chan;
%             miscinfo{s+1,8}=EEG.chanlocs(max_ind_chan).labels;
%             miscinfo{s+1,9}=max_val_chan;
%             
%             % save it with original filename but get rid of original
%             % extention (hence the 1:end-4)
%             save([condfiles_subs{2}{s}(1:end-4),'_',options.measure,'_extracted.mat'],'data');
%             clear data datamax datamin
%             
%             % load cond 1 file
%             EEG = pop_loadset('filename',condfiles_subs{1}{s},'filepath',pathtofiles{1});
%             EEG = eeg_checkset(EEG);
%             miscinfo{s+1,2}=min_ind_chan;
%             miscinfo{s+1,3}=EEG.chanlocs(min_ind_chan).labels;
%             miscinfo{s+1,4}=min_val_chan;
%             
%             % now we have max and min channel index, steal it from original
%             % data, from cond1.
%             datamax=EEG.data(max_ind_chan,:,:);
%             datamin=EEG.data(min_ind_chan,:,:);
%             
%             % create bipolar
%             data=datamax-datamin;
%             
%             % save it with original filename but get rid of original
%             % extention (hence the 1:end-4)
%             save([condfiles_subs{1}{s}(1:end-4),'_',options.measure,'_extracted.mat'],'data');
%             clear data datamax datamin
%             
%         end
%         
%         STATS.miscinfo=miscinfo;
%         xtimes=EEG.times;
%         disp(' ***** finished extracting the data array for bipolar analysis *****')
        
    case 'scalpchan'
        disp('***** extracting selected channel(s). Multiple channels will be averaged together in the next step ***** ')
        
        for k=1:numconds;
            
            [rowcond colcond]=size(condfiles_subs{k});
            
            for s=1:rowcond; % scroll through subjects
                
                % load file
                EEG = pop_loadset('filename',condfiles_subs{k}{s},'filepath',pathtofiles{k});
                EEG = eeg_checkset(EEG);
                
                if persubject==0;
                    
                    try
                        % scroll through chans the user wants and collect relavent indices
                        for i=1:length(options.chans)
                            chanind(i)=find(strcmp({EEG.chanlocs.labels},options.chans{i}));
                        end
                        
                    catch
                        error(['Channel ' options.chans{i}, 'does not exist for subject ', condfiles_subs{k}{s}]);
                    end
                    
                    
                elseif persubject==1;
                    
                    % when some option to select a channel is given
                    try
                        % scroll through chans the user wants and collect relavent indices
                        for ii=1:length(chanfile{s,k+k})
                            chanind(ii)=find(strcmp({EEG.chanlocs.labels},chanfile{s,k+k}{ii}));
                        end
                        
                    catch
                        error(['Channel ' chanfile{s,k+k}{ii}, 'does not exist for subject ', condfiles_subs{k}{s}]);
                    end
                    
                end
                
                
                % steal data from the channels you are interested in
                data=EEG.data(chanind,:,:);
                
                % save it with original filename but get rid of original
                % extention (hence the 1:end-4)
                save([condfiles_subs{k}{s}(1:end-4),'_',options.measure,'_extracted.mat'],'data');
                clear data
                
            end
        end
        
        xtimes=EEG.times;
        disp('***** finished extracting selected channel(s) ***** ')
        
%     case 'scalpchansub'
%         disp('***** extracting selected channel(s). Multiple channels will be averaged together in the next step ***** ')
% 
%         for k=1:numconds;
%             
%             [rowcond colcond]=size(condfiles_subs{k});
%             
%             if k==1;
%                 
%                 miscinfo=cell(rowcond+1,4);
%                 miscinfo(1,:)={'sub', 'ind', 'chan', 'val'};
%                 miscinfo(2:end,1)=condfiles_subs{1}(:);
%                 
%                 for s=1:rowcond; % scroll through subjects
%                     
%                     % load file from first condition
%                     EEG = pop_loadset('filename',condfiles_subs{k}{s},'filepath',pathtofiles{k});
%                     EEG = eeg_checkset(EEG);
%                     
%                     % determin the largest neg channel during window of interest
%                     h=figure; title(['condition as, subject ', num2str(s)]);
%                     pop_timtopo(EEG, [-100 400], [NaN], ['condition as, subject ', num2str(s)],'electrodes','on');
%                     
%                     timewin=input('type in timewindow around N170\n');
%                     close(h);
%                     
%                     MStoTF_min=round((timewin(1)/1000-EEG.xmin)/(EEG.xmax-EEG.xmin) * (EEG.pnts-1))+1;
%                     MStoTF_max=round((timewin(2)/1000-EEG.xmin)/(EEG.xmax-EEG.xmin) * (EEG.pnts-1))+1;
%                     
%                     % for each channel get min value within a window
%                     [val_time ind_time]=min(mean(EEG.data(:,MStoTF_min:MStoTF_max,:),3),[],2);
%                     
%                     % then, get the index for the min channel
%                     [val_chan(s) ind_chan(s)]=min(val_time);
%                     
%                     % add chanlabel to speadsheet
%                     miscinfo{s+1,2}=ind_chan(s);
%                     miscinfo{s+1,4}=val_chan(s);
%                     
%                     chanlab{s}=EEG.chanlocs(ind_chan(s)).labels;
%                     miscinfo{s+1,3}=chanlab{s};
%                     
%                     % steal data from the channels you are interested in
%                     data=EEG.data(ind_chan(s),:,:);
%                     
%                     % save it with original filename but get rid of original
%                     % extention (hence the 1:end-4)
%                     save([condfiles_subs{k}{s}(1:end-4),'_',savestring,'_',options.measure,'_extracted.mat'],'data');
%                     clear data
%                     
%                 end
%                 
%             end
%             
%             if k>1;
%                 
%                 for s=1:rowcond; % scroll through subjects
%                     
%                     % load file from first condition
%                     EEG = pop_loadset('filename',condfiles_subs{k}{s},'filepath',pathtofiles{k});
%                     EEG = eeg_checkset(EEG);
%                     
%                     % determin the largest neg channel during window of interest
%                     %figure; h=pop_timtopo(EEG, [-100 400], [NaN], 'ERP data and scalp maps of as','electrodes','on');
%                     
%                     %chanlab{s}=input('type in chanlabel you want');
%                     %close(h);
%                     
%                     % add chanlabel to speadsheet
%                     %miscinfo(s,3)=chanlab;
%                     
%                     %chanind(s)=find(strcmp({EEG.chanlocs.labels},chanlab));
%                     %miscinfo(s,2)=chanind;
%                     
%                     
%                     % steal data from the channels you are interested in
%                     data=EEG.data(ind_chan(s),:,:);
%                     
%                     % save it with original filename but get rid of original
%                     % extention (hence the 1:end-4)
%                     save([condfiles_subs{k}{s}(1:end-4),'_',savestring,'_',options.measure,'_extracted.mat'],'data');
%                     clear data
%                     
%                 end
%                 
%             end
%             
%         end
%         STATS.miscinfo=miscinfo;
%         xtimes=EEG.times;
%         disp('***** finished extracting selected channel(s) ***** ')
        
        
    case 'scalpersp'
        disp('***** extracting selected channel(s). Multiple channels will be averaged together in the next step ***** ')
        
        for k=1:numconds;
            
            [rowcond colcond]=size(condfiles_subs{k});
            
            for s=1:rowcond; % scroll through subjects
                
                % load file
                EEG = pop_loadset('filename',condfiles_subs{k}{s},'filepath',pathtofiles{k});
                EEG = eeg_checkset(EEG);
                
                
                 if persubject==0;
                    
                    try
                        % scroll through chans the user wants and collect relavent indices
                        for i=1:length(options.chans)
                            chanind(i)=find(strcmp({EEG.chanlocs.labels},options.chans{i}));
                        end
                        
                    catch
                        error(['Channel ' options.chans{i}, 'does not exist for subject ', condfiles_subs{k}{s}]);
                    end
                    
                    
                elseif persubject==1;
                    
                    % when some option to select a channel is given
                    try
                        % scroll through chans the user wants and collect relavent indices
                        for ii=1:length(chanfile{s,k+k})
                            chanind(ii)=find(strcmp({EEG.chanlocs.labels},chanfile{s,k+k}{ii}));
                        end
                        
                    catch
                        error(['Channel ' chanfile{s,k+k}{ii}, 'does not exist for subject ', condfiles_subs{k}{s}]);
                    end
                    
                end
                
                % steal data from the channels you are interested in
                data=EEG.data(chanind,:,:);
                
                % save it with original filename but get rid of original
                % extention (hence the 1:end-4)
                save([condfiles_subs{k}{s}(1:end-4),'_',options.measure,'_extracted.mat'],'data');
                clear data
                
            end
            
        end
        xtimes=EEG.times;
        disp('***** finished extracting selected channel(s) ***** ')
        
        
    case 'scalpitc'
        disp('***** extracting selected channel(s). Multiple channels will be averaged together in the next step ***** ')
        
        for k=1:numconds;
            
            [rowcond colcond]=size(condfiles_subs{k});
            
            for s=1:rowcond; % scroll through subjects
                
                % load file
                EEG = pop_loadset('filename',condfiles_subs{k}{s},'filepath',pathtofiles{k});
                EEG = eeg_checkset(EEG);
                
                if persubject==0;
                    
                    try
                        % scroll through chans the user wants and collect relavent indices
                        for i=1:length(options.chans)
                            chanind(i)=find(strcmp({EEG.chanlocs.labels},options.chans{i}));
                        end
                        
                    catch
                        error(['Channel ' options.chans{i}, 'does not exist for subject ', condfiles_subs{k}{s}]);
                    end
                    
                    
                elseif persubject==1;
                    
                    % when some option to select a channel is given
                    try
                        % scroll through chans the user wants and collect relavent indices
                        for ii=1:length(chanfile{s,k+k})
                            chanind(ii)=find(strcmp({EEG.chanlocs.labels},chanfile{s,k+k}{ii}));
                        end
                        
                    catch
                        error(['Channel ' chanfile{s,k+k}{ii}, 'does not exist for subject ', condfiles_subs{k}{s}]);
                    end
                    
                end
                
                % steal data from the channels you are interested in
                data=EEG.data(chanind,:,:);
                
                % save it with original filename but get rid of original
                % extention (hence the 1:end-4)
                save([condfiles_subs{k}{s}(1:end-4),'_',options.measure,'_extracted.mat'],'data');
                clear data
                
            end
            
        end
        xtimes=EEG.times;
        disp('***** finished extracting selected channel(s) ***** ')
end


% populate the all important STATS structure
STATS.orig_setfiles=condfiles_subs;
STATS.pathtofiles=pathtofiles;
STATS.levels=levels;
STATS.design=design;
STATS.savestring=savestring;
STATS.xtimes=xtimes;
STATS.numpnts=size(STATS.xtimes,2);
STATS.condnames=condnames;
STATS.datatype=options.measure;
STATS.numconds=numconds;
STATS.srate=EEG.srate;
STATS.xmin=EEG.xmin;
STATS.xmax=EEG.xmax;
STATS.tfcycles=options.tfcycles;
STATS.freqs=options.freqs;
STATS.nfreqs=options.nfreqs;
%STATS.timesout=options.timesout;
STATS.tfbsline=options.tfbsline;
STATS.chanlabels=options.chans;


% set measures for following resampling procedure
if any(strcmp({'icascalp','scalpchan','scalpchansub'},options.measure))
    STATS.measure='chanclust';
elseif any(strcmp({'icagfa','scalpgfa'},options.measure))
    STATS.measure='gfa';
elseif any(strcmp({'icaersp','scalpersp'},options.measure))
    STATS.measure='ersp';
elseif any(strcmp({'scalpitc','icaitc'},options.measure))
    STATS.measure='itc';
end

try STATS.ICretain=ICretain; catch, STATS.ICretain=[]; end
try STATS.chaninfo=chanfile; catch, STATS.chaninfo=[]; end

disp('******* Saving STATS structure *******')
save(['STATS_',savestring,'.mat'],'STATS');

end



