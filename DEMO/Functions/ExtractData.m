%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This function extracts data in various forms (e.g., ICs, gfa, scalp, see below)
from EEGLAB's .set files. The output files are then used in ResampleData.m with the
ultimate goal of calculating percentile bootstrap statistics at every timepoint
for single-subjects as well as the group (SubjectStatistics.m, GroupStatistics.m),
given any design (up to 2-way).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Inputs:


condfiles - cell array of filenames for each subject and condition.
            Leave emtpy (ie., []) and MATLAB will bring up an interface for you to
            load the appropriate subject condition files. After this is done, the
            file is saved (e.g., confiles_scalpchan.mat) and can be entered for subsequent
            calls instead of using [] and having to manually load.

condnames - cell array of condition labels. For example, {'face' 'house' 'object'}

varargin  - key/val pairs. See Options.


Options:

icamax    - Project selected ICs to scalp channel with maximum weight determined
            by the weight matrix. All ICA options must be followed with
            a filename of text file (csv) indicating
            components to retain for each subject. See demo data to
            understand the construction of this file.

icagfa    - Project selected ICs to the scalp and extract entire data array
            for later GFA calculations. All ICA options must be followed with
            a filename of text file (csv) indicating
            components to retain for each subject. See demo data to
            understand the construction of this file.

scalpgfa   - Extract the full scalp data array for later GFA
             calculations. Does not require any following input arguments.

scalpchan  - Extract specified channel(s). Must be followed with cell array
             of channel labels. E.g., {'C11' 'C15' 'C12'}. If asking for
             more than one channel, the channel group will be averaged
             together during the resampling stages (ResampleData.m) giving
             you a channel cluster. This function does not handle doing
             stats on multiple channels independently in the same call.

Examples:
[STATS]=ExtractData({'AE' 'AH' 'SE' 'SH'},[],[2 2],'ww','Occipital_Analysis','scalpchan',{'A23', 'A24' })
[STATS]=ExtractData({'AE' 'AH' 'SE? ?SH?}, ?condfiles_icamax_analysis.mat?, [2 2],'ww','icamax_analysis','icamax',?wwICfile.txt?)


Copyright (C) <2015>  <Allan Campopiano>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [STATS]=ExtractData(condnames,condfiles,levels,design,savestring,varargin)

% basic error handling for varargin options
nargs=length(varargin);
if nargs==1 && strcmp(varargin,'scalpgfa')
    varargin{2}=[];
    
elseif round(nargs/2)~=nargs/2
    error('need NAME/VALUE pairs for optional inputs')
    
elseif nargs>2
    error('you are asking to extract too many things. See Options for acceptable choices')
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
    save(['condfiles_',varargin{1},'.mat'],'condfiles');
    
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

% if varargin{1} is this, do that, etc
switch varargin{1}
    
    case 'icamax'
        disp(' ***** projecting selected ICs to scalp and extracting channel with maximum weight ***** ')
        
        % load txt file with ICs to retain
        % trailing unequal positions in ICretain will be filled with zeros,
        % normally this is problematic but since we only loop through any
        % given column rowcond (number of subjects in a condition) times,
        % it should be okay (i.e., the zero is never accessed).
        ICretain=csvread(varargin{2});
        
        %if strcmp(design,'bw');
        if strcmp(design,'bw');
            
            kk=0;
            jj=1;
            for j=1:jlvls;
                
                [rowcond colcond]=size(condfiles_subs{jj});
                
                for s=1:rowcond; % scroll through subjects
                    
                    for k=1:klvls;
                        
                        % load file
                        EEG = pop_loadset('filename',condfiles_subs{k+kk}{s},'filepath',pathtofiles{k+kk});
                        EEG = eeg_checkset(EEG);
                        
                        % get orignal IC inds
                        icorig=[1:min(size(EEG.icawinv))];
                        
                        compswant=ICretain(s,j);
                        
                        icorig(compswant)=[];
                        
                        % pick ICs you want to retain
                        EEGretain = pop_subcomp( EEG, icorig, 0);
                        EEGretain = eeg_checkset(EEGretain);
                        
                        % figure out which channel has max weight
                        for i=1:size(EEGretain.icawinv,2);
                            [maxweight(i) maxind(i)]=max(abs(EEGretain.icaweights(i,:)));
                        end
                        
                        % put EEG.data into the variable data
                        data=EEGretain.data(maxind,:,:);
                        
                        % save it with original filename but get rid of original
                        % extention (hence the 1:end-4)
                        save([condfiles_subs{k+kk}{s}(1:end-4),'_',varargin{1},'_extracted.mat'],'data');
                        clear data
                        
                        % get residual
                        EEGres = pop_subcomp( EEG, compswant, 0);
                        EEGres = eeg_checkset(EEGres);
                        
                        % put EEG.data into the variable data
                        data=EEGres.data(maxind,:,:);
                        
                        % save it with original filename but get rid of original
                        % extention (hence the 1:end-4)
                        save([condfiles_subs{k+kk}{s}(1:end-4),'_',varargin{1},'_RESextracted.mat'],'data');
                        clear data
                        
                    end
                    
                end
                
                kk=kk+klvls;
                jj=jj+klvls;
            end
            
            %elseif strcmp(design,'ww')||strcmp(design,'w');
        elseif strcmp(design,'ww')||strcmp(design,'w');
            
            % always assuming that within subjects factors have equal
            % number of subjects. So I only query the size of condfiles_subs{1}
            [rowcond colcond]=size(condfiles_subs{1});
            
            for s=1:rowcond; % scroll through subjects
                
                for k=1:numconds;
                    
                    % load file
                    EEG = pop_loadset('filename',condfiles_subs{k}{s},'filepath',pathtofiles{k});
                    EEG = eeg_checkset(EEG);
                    
                    % get orignal IC inds
                    icorig=[1:min(size(EEG.icawinv))];
                    
                    compswant=ICretain(s);
                    
                    icorig(compswant)=[];
                    
                    % pick ICs you want to retain
                    EEGretain = pop_subcomp( EEG, icorig, 0);
                    EEGretain = eeg_checkset(EEGretain);
                    
                    % figure out which channel has max weight
                    for i=1:size(EEGretain.icawinv,2);
                        [maxweight(i) maxind(i)]=max(abs(EEGretain.icaweights(i,:)));
                    end
                    
                    % put EEG.data into the variable data
                    data=EEGretain.data(maxind,:,:);
                    
                    % save it with original filename but get rid of original
                    % extention (hence the 1:end-4)
                    save([condfiles_subs{k}{s}(1:end-4),'_',varargin{1},'_extracted.mat'],'data');
                    clear data
                    
                    % get residual
                    EEGres = pop_subcomp( EEG, compswant, 0);
                    EEGres = eeg_checkset(EEGres);
                    
                    % put EEG.data into the variable data
                    data=EEGres.data(maxind,:,:);
                    
                    % save it with original filename but get rid of original
                    % extention (hence the 1:end-4)
                    save([condfiles_subs{k}{s}(1:end-4),'_',varargin{1},'_RESextracted.mat'],'data');
                    clear data
                    
                end
                
            end
            
            %elseif strcmp(design,'bb')||strcmp(design,'b');
        elseif strcmp(design,'bb')||strcmp(design,'b');
            
            for k=1:numconds;
                
                [rowcond colcond]=size(condfiles_subs{k});
                
                for s=1:rowcond; % scroll through subjects
                    
                    % load file
                    EEG = pop_loadset('filename',condfiles_subs{k}{s},'filepath',pathtofiles{k});
                    EEG = eeg_checkset(EEG);
                    
                    % get orignal IC inds
                    icorig=[1:min(size(EEG.icawinv))];
                    
                    compswant=ICretain(s,k);
                    
                    icorig(compswant)=[];
                    
                    % pick ICs you want to retain
                    EEGretain = pop_subcomp( EEG, icorig, 0);
                    EEGretain = eeg_checkset(EEGretain);
                    
                    % figure out which channel has max weight
                    for i=1:size(EEGretain.icawinv,2);
                        [maxweight(i) maxind(i)]=max(abs(EEGretain.icaweights(i,:)));
                    end
                    
                    % put EEG.data into the variable data
                    data=EEGretain.data(maxind,:,:);
                    
                    % save it with original filename but get rid of original
                    % extention (hence the 1:end-4)
                    save([condfiles_subs{k}{s}(1:end-4),'_',varargin{1},'_extracted.mat'],'data');
                    clear data
                    
                    % get residual
                    EEGres = pop_subcomp( EEG, compswant, 0);
                    EEGres = eeg_checkset(EEGres);
                    
                    % put EEG.data into the variable data
                    data=EEGres.data(maxind,:,:);
                    
                    % save it with original filename but get rid of original
                    % extention (hence the 1:end-4)
                    save([condfiles_subs{k}{s}(1:end-4),'_',varargin{1},'_RESextracted.mat'],'data');
                    clear data
                    
                end
                
            end
            
        end
        xtimes=EEG.times;
        disp(' ***** finished projecting selected ICs to scalp and extracting channel with maximum weight ***** ')
        
    case 'icagfa'
        disp(' ***** projecting selected ICs to scalp and extracting data array for later GFA calculation ***** ')
        
        % load txt file with ICs to retain
        % trailing unequal positions in ICretain will be filled with zeros,
        % normally this is problematic but since we only loop through any
        % given column rowcond (number of subjects in a condition) times,
        % it should be okay (i.e., the zero is never accessed).
        ICretain=csvread(varargin{2});
        
        % if strcmp(design,'bw');
        if strcmp(design,'bw');
            
            kk=0;
            jj=1;
            for j=1:jlvls;
                
                [rowcond colcond]=size(condfiles_subs{jj});
                
                for s=1:rowcond; % scroll through subjects
                    
                    for k=1:klvls;
                        
                        % load file
                        EEG = pop_loadset('filename',condfiles_subs{k+kk}{s},'filepath',pathtofiles{k+kk});
                        EEG = eeg_checkset(EEG);
                        
                        % get orignal IC inds
                        icorig=[1:min(size(EEG.icawinv))];
                        
                        compswant=ICretain(s,j);
                        
                        icorig(compswant)=[];
                        
                        % pick ICs you want to retain
                        EEGretain = pop_subcomp( EEG, icorig, 0);
                        EEGretain = eeg_checkset(EEGretain);
                        
                        % put EEG.data into the variable data
                        data=EEGretain.data;
                        
                        % save it with original filename but get rid of original
                        % extention (hence the 1:end-4)
                        save([condfiles_subs{k+kk}{s}(1:end-4),'_',varargin{1},'_extracted.mat'],'data');
                        clear data
                        
                        % get residual
                        EEGres = pop_subcomp( EEG, compswant, 0);
                        EEGres = eeg_checkset(EEGres);
                        
                        % put EEG.data into the variable data
                        data=EEGres.data;
                        
                        % save it with original filename but get rid of original
                        % extention (hence the 1:end-4)
                        save([condfiles_subs{k+kk}{s}(1:end-4),'_',varargin{1},'_RESextracted.mat'],'data');
                        clear data
                        
                    end
                    
                end
                
                kk=kk+klvls;
                jj=jj+klvls;
            end
            
            %elseif strcmp(design,'ww')||strcmp(design,'w');
        elseif strcmp(design,'ww')||strcmp(design,'w');
            
            % always assuming that within subjects factors have equal
            % number of subjects. So I only query the size of condfiles_subs{1}
            [rowcond colcond]=size(condfiles_subs{1});
            
            for s=1:rowcond; % scroll through subjects
                
                for k=1:numconds;
                    
                    % load file
                    EEG = pop_loadset('filename',condfiles_subs{k}{s},'filepath',pathtofiles{k});
                    EEG = eeg_checkset(EEG);
                    
                    % get orignal IC inds
                    icorig=[1:min(size(EEG.icawinv))];
                    
                    compswant=ICretain(s);
                    
                    icorig(compswant)=[];
                    
                    % pick ICs you want to retain
                    EEGretain = pop_subcomp( EEG, icorig, 0);
                    EEGretain = eeg_checkset(EEGretain);
                    
                    % put EEG.data into the variable data
                    data=EEGretain.data;
                    
                    % save it with original filename but get rid of original
                    % extention (hence the 1:end-4)
                    save([condfiles_subs{k}{s}(1:end-4),'_',varargin{1},'_extracted.mat'],'data');
                    clear data
                    
                    % get residual
                    EEGres = pop_subcomp( EEG, compswant, 0);
                    EEGres = eeg_checkset(EEGres);
                    
                    % put EEG.data into the variable data
                    data=EEGres.data;
                    
                    % save it with original filename but get rid of original
                    % extention (hence the 1:end-4)
                    save([condfiles_subs{k}{s}(1:end-4),'_',varargin{1},'_RESextracted.mat'],'data');
                    clear data
                    
                end
                
            end
            
            %elseif strcmp(design,'bb')||strcmp(design,'b');
        elseif strcmp(design,'bb')||strcmp(design,'b');
            
            for k=1:numconds;
                
                [rowcond colcond]=size(condfiles_subs{k});
                
                for s=1:rowcond; % scroll through subjects
                    
                    % load file
                    EEG = pop_loadset('filename',condfiles_subs{k}{s},'filepath',pathtofiles{k});
                    EEG = eeg_checkset(EEG);
                    
                    % get orignal IC inds
                    icorig=[1:min(size(EEG.icawinv))];
                    
                    compswant=ICretain(s,k);
                    
                    icorig(compswant)=[];
                    
                    % pick ICs you want to retain
                    EEGretain = pop_subcomp( EEG, icorig, 0);
                    EEGretain = eeg_checkset(EEGretain);
                    
                    % put EEG.data into the variable data
                    data=EEGretain.data;
                    
                    % save it with original filename but get rid of original
                    % extention (hence the 1:end-4)
                    save([condfiles_subs{k}{s}(1:end-4),'_',varargin{1},'_extracted.mat'],'data');
                    clear data
                    
                    % get residual
                    EEGres = pop_subcomp( EEG, compswant, 0);
                    EEGres = eeg_checkset(EEGres);
                    
                    % put EEG.data into the variable data
                    data=EEGres.data;
                    
                    % save it with original filename but get rid of original
                    % extention (hence the 1:end-4)
                    save([condfiles_subs{k}{s}(1:end-4),'_',varargin{1},'_RESextracted.mat'],'data');
                    clear data
                    
                end
                
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
                save([condfiles_subs{k}{s}(1:end-4),'_',varargin{1},'_extracted.mat'],'data');
                clear data
                
            end
            
        end
        xtimes=EEG.times;
        disp(' ***** finished extracting the data array for later GFA calculations *****')
        
    case 'scalpchan'
        disp('***** extracting selected channel(s). Multiple channels will be averaged together in the next step ***** ')
        
        for k=1:numconds;
            
            [rowcond colcond]=size(condfiles_subs{k});
            
            for s=1:rowcond; % scroll through subjects
                
                % load file
                EEG = pop_loadset('filename',condfiles_subs{k}{s},'filepath',pathtofiles{k});
                EEG = eeg_checkset(EEG);
                
                try
                    % scroll through chans the user wants and collect relavent indices
                    for i=1:length(varargin{2})
                        chanind(i)=find(strcmp({EEG.chanlocs.labels},varargin{2}{i}));
                    end
                    
                catch
                    error(['One of the channels you want does not exist for some subjects. ' ...
                        'If the files are interpolated, make sure they are all interpolated to the same montage. ' ...
                        'If they are not interpolated, the channels you are looking for must exist for every subject. ' ...
                        'You must have channel information loaded into the EEGLABs edit channel locations pane. ' ...
                        'If you have not loaded channel information (see the blue EEGLAB pane for info), google how to do that.']);
                end
                
                % steal data from the channels you are interested in
                data=EEG.data(chanind,:,:);
                
                % save it with original filename but get rid of original
                % extention (hence the 1:end-4)
                save([condfiles_subs{k}{s}(1:end-4),'_',varargin{1},'_extracted.mat'],'data');
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
STATS.datatype=varargin{1};
STATS.numconds=numconds;
STATS.srate=EEG.srate;

if strcmp('scalpchan',varargin{1})
    STATS.chanlabels=varargin{2};
else
    STATS.chanlabels=[];
end

if any(strcmp({'icamax','scalpchan'},varargin{1}))
    STATS.measure='chanclust';
elseif any(strcmp({'icagfa','scalpgfa'},varargin{1}))
    STATS.measure='gfa';
end

try STATS.ICretain=ICretain; catch, STATS.ICretain=[]; end

disp('******* Saving STATS structure *******')
save(['STATS_',savestring,'.mat'],'STATS');

end



