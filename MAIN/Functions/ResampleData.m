% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function bootstraps single trial EEG data for each subject. The 20% trimmed mean is taken
% across trials (i.e., trimmed ERPS) at each timepoint. Simply load in as many NxMxP .mat files
% as you wish when the browser comes up. Each file is saved with the original file name with
% 'bootstrapped' appended to it. This function can take some time to run depending on the number
% of bootstrapp samples as well as the size and number of files. However, your computer will
% likely not be rendered usless while it runs as not much RAM is required, just time.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input arguments:
%     STATS stucture - you will be prompted to load this
%     nboot = number of resamples you wish to take from the single-trials
%
% varargin - key/val pairs, see Options for details
%
% Options:
%       trialcap - a string paired with the number of the cap (e.g., 50).
%       The trialcap option limits the number of trials used in each resample
%       to the cap. This helps to equate noise levels due to unequal number of trials in
%       different conditions. Usually this only a problem when using GFA
%       measures IF there are highly dissproportionate numbers of trials
%       across conditions, as one condition with more noise (less trials)
%       will offset the GFA waveform compared to a condition with less noise (more trials).
%       Typical ERPs occilate around zero, so trialcap does not likely need to be set for
%       ERPs. Omit this option alltogether if you wish to resample from the
%       max number of trials (default).
%
%
% Examples:
%     ResampleData(STATS,1000,'trialcap', 50);
%     ResampleData(STATS,50000); % will take a while
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
%%
function [STATS]=ResampleData(STATS,condfiles,nboot,varargin)

% bring in STATS stucture
if isempty(STATS)
    [fnamestats]=uigetfile('*.mat','Select the STATS structure for this analysis');
    load(fnamestats);
else
    load(STATS);
end

% add field to STATS
STATS.nboot=nboot;

% set defaults, with only one option this makes no sense, but if more get
% added it is exandable in this for
options.trialcap='none';

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

if strcmp(options.trialcap,'none') || isempty(varargin);
    disp('***** not using trial cap *****');
    capflag=0;
else
    STATS.trialcap=options.trialcap;
    disp(['***** using trial cap of ', num2str(STATS.trialcap), ' *****']);
    capflag=1;
end

% this is a '1 X number of files' cell array of strings (fnames)
if isempty(condfiles)
    fnames=uigetfile('*.mat','Select all files you wish to bootstrap', 'MultiSelect','on');
    
else
    condfiles=load(condfiles);
    condfields=fieldnames(condfiles);
    fnames=condfiles.(condfields{1});
end


% set a function handle to whichever measure the user specified ('gfa' or 'chanclust')
if strcmp(STATS.measure,'gfa')
    measure_han=@(x) std(x,1);
elseif strcmp(STATS.measure,'chanclust');
    measure_han=@(x) mean(x,1);
else
    error('must choose either gfa or chanclust as input arguments');
end

[rowfile colfile]=size(fnames);
datacell=cell(1,colfile);

h2 = waitbar(0,'1','Name','Processing Files','Position',[1100 549 550 40]);
childh2 = get(h2, 'Children');
set(childh2, 'Position',[5 10 538 15]);

h3 = waitbar(0,'1','Name','Bootstrapping','Position',[1100 486 550 40]);
childh3 = get(h3, 'Children');
set(childh3, 'Position',[5 10 538 15]);

% loading data
for filecurrent=1:colfile;
    subdata=load(fnames{1,filecurrent});
    
    %datacell{1,filecurrent}=subdata.data;
    datacell{1}=subdata.data;
    
    
    % data is set to zero and the size is taken at the start of each new file
    
    %[jnk numpnts pageEEG]=size(datacell{filecurrent});
    [jnk numpnts pageEEG]=size(datacell{1});
    trialEEG=pageEEG;
    
    
    % deal with trial cap if there is one
    if capflag
        if trialEEG>STATS.trialcap
            trialEEG=STATS.trialcap;
        end
    end
    
    switch STATS.measure
        
        case {'gfa', 'chanclust'}
            
            data=zeros(nboot,numpnts);
            
            % boot loop
            for bootcurrent=1:nboot;
                
                % resample with replacement from datacell, creating bootcell
                bootvect=randi(pageEEG,1,trialEEG);
                
                %bootcell=datacell{filecurrent}(:,:,bootvect);
                bootcell=datacell{1}(:,:,bootvect);
                
                % trimmed channel ERPs
                chan_erps=trimmean(bootcell,40,3);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % this statement makes use of function handles.
                user_measure=measure_han(chan_erps);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %store bootstrap surrogates
                data(bootcurrent,:)=user_measure;
                
                waitbar(bootcurrent/nboot,h3,sprintf('%12s',[num2str(bootcurrent),'/',num2str(nboot)]))
            end
            
        case {'ersp', 'itc'}
            
            if length(datacell{1},1)>1;
                chanlength='multi';
            else
                chanlength='single';
            end
            
            % boot loop
            for bootcurrent=1:nboot;
                
                % resample with replacement from datacell, creating bootcell
                bootvect=randi(pageEEG,1,trialEEG);
                
                switch chanlength
                    case 'multi'
                                               
                        for i=1:length(datacell{1},1) % num rows/channels
                            
                            % compute TF info to get coeficients for all channels
                            figure('visible','off');
                            [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = newtimef(datacell{1}(i,:,:), ...
                                STATS.xtimes, [STATS.xmin STATS.xmax]*1000, ...
                                STATS.srate, STATS.tfcycles,'freqs',STATS.freqs,'timesout',STATS.timesout,'plotersp','off','plotitc','off');
                            % STATS.srate, [3 .5],'freqs',[3 30],'timesout',600);
                            clear ersp itc powbase times freqs erspboot itcboot
                            
                            if strcmp(STATS.measure, 'ersp')
                                clear itc powbase times freqs erspboot itcboot
                                
                                mapwrite(trimmean(abs(tfdata(:,:,bootvect)),40,3),[fnames{1,filecurrent}, '_tempbootstrapped.map'],'datsize',[STATS.TF_freqs,STATS.TF_times,length(datacell{1})]);
                                
                            elseif strcmp(STATS.measure, 'itc')
                                clear itc powbase times freqs erspboot itcboot
                                
                                % no trimmed mean?
                                mapwrite(abs(mean(((tfdata(:,:,bootvect,:))./abs(tfdata(:,:,bootvect,:))),4)),[fnames{1,filecurrent}, '_tempbootstrapped.map'],'datsize',[STATS.TF_freqs,STATS.TF_times,length(datacell{1})]);
                                
                            end
                            
                        end
                        
                        % this holds multiple channel spectral information
                        datamap=mapread([fnames{1,filecurrent}, '_tempbootstrapped.map'],'dat','datsize',[STATS.TF_freqs,STATS.TF_times,length(datacell{1})]);
                        
                        % this is the avreage of the multi channel spectral information, like a spectral ROI
                        mapwrite(mean(datamap.Data.dat,3),[fnames{1,filecurrent}, '_bootstrapped.mat'],'datsize',[STATS.TF_freqs,STATS.TF_times,STATS.nboot]);
                        
                        % get rid or gathering arrays
                        delete([fnames{1,filecurrent}, '_tempbootstrapped.map']);
                        
                    case 'single'
                        
                        % compute TF info to get coeficients
                        figure('visible','off');
                        [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = newtimef(datacell{1}(1,:,:), ...
                            STATS.xtimes, [STATS.xmin STATS.xmax]*1000, ...
                            STATS.srate, STATS.tfcycles,'freqs',STATS.freqs,'timesout',STATS.timesout,'plotersp','off','plotitc','off');
                        % STATS.srate, [3 .5],'freqs',[3 30],'timesout',600);
                        clear ersp itc powbase times freqs erspboot itcboot
                        
                        if strcmp(STATS.measure, 'ersp')
                            clear itc powbase times freqs erspboot itcboot
                            
                            mapwrite(trimmean(abs(tfdata(:,:,bootvect)),40,3),[fnames{1,filecurrent}, '_bootstrapped.mat'],'datsize',[STATS.TF_freqs,STATS.TF_times,STATS.nboot]);
                            
                        elseif strcmp(STATS.measure, 'itc')
                            clear itc powbase times freqs erspboot itcboot
                            
                            % no trimmed mean
                            mapwrite(abs(mean(((tfdata(:,:,bootvect,:))./abs(tfdata(:,:,bootvect,:))),4)),[fnames{1,filecurrent}, '_bootstrapped.mat'],'datsize',[STATS.TF_freqs,STATS.TF_times,STATS.nboot]);
                            
                        end
                        
                end
                
                waitbar(bootcurrent/nboot,h3,sprintf('%12s',[num2str(bootcurrent),'/',num2str(nboot)]))
            end
    end
    
    
    if capflag % this can probably be taken out after system comparison
        
        % the string in the sqare brackets is what it appended to the original file name
        save([fnames{filecurrent}(1,1:end-4),'_bootstrapped_', num2str(STATS.trialcap),'.mat'],'data');
    else
        
        save([fnames{filecurrent}(1,1:end-4),'_bootstrapped.mat'],'data');
    end
    
    waitbar(filecurrent/colfile,h2,sprintf('%12s',[num2str(filecurrent),'/',num2str(colfile)]))
    
end
close(h2,h3);

disp('******* Finished resmapling from single trials *******')
disp('******* Saving STATS structure *******')
save(['STATS_',STATS.savestring,'.mat'],'STATS');
end

