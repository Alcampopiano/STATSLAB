% This function resamples single trial EEG data for each input file. The 20% trimmed mean is taken
% across trials at each timepoint. Each output file is saved with the original file name with
% “bootstrapped” appended to it.
% 
% Inputs:
% 
% ***condfiles***
% Leave empty and MATLAB will bring up an interface for you to load the appropriate “extracted” files.
% ***end***
% 
% ***nboot*** 
% number of resamples you wish to take from the single-trials ***end***
% 
% ***varargin***
% Options are specified in pairs (key -> val)
% 
% trialcap ->
% 	
% 	[numeric] - the trialcap option limits the # of trials used in each resample to the specified value
% 
% 	none - sample a new set of trials that is equal in size to the original set 
% 	
% For example,
% 
% trialcap
% 100
% 
% will resample from the single-trial data, randomly choosing with replacement 100 trials.
% 
% Using ResampleData at the commandline:
% 
% ResampleData(STATS, [], 1000, 'trialcap', 'none');
% ***end***
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
    
    % save so one can load without gui
    save(['condfiles_ResampleData_',STATS.measure,'_',STATS.savestring,'.mat'],'fnames');
    
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
    %else
    %    error('must choose either gfa or chanclust as input arguments');
end

[rowfile colfile]=size(fnames);
datacell=cell(1,colfile);

h2 = waitbar(0,'1','Name','Processing Files','Position',[1100 549 550 40]);
childh2 = get(h2, 'Children');
set(childh2, 'Position',[5 10 538 15]);

h3 = waitbar(0,'1','Name','Bootstrapping','Position',[1100 486 550 40]);
childh3 = get(h3, 'Children');
set(childh3, 'Position',[5 10 538 15]);


% get rid of previously mapped bootstrapped files
for filecurrent=1:colfile;
    warning off
    delete([fnames{1,filecurrent}(1,1:end-4), '_bootstrapped.map']);
    warning on
end

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
        
%         %%% trialcap edit when doing trial cap on a per subject basis
%         %%%%%%%%%%%%%%%%%%%%
%         load('trialmins.mat')
%         load('trialpre.mat')
%         
%         for cp=1:length(trialmins)
%             if strfind(fnames{1,filecurrent},trialpre{cp})
%                 trialEEG=trialmins(cp);
%             end
%         end       
        
        %%%% this should be here when setting cap that applies to
        %%%% everyone, although the above edit should be made as an option
        %%%% so that trial caps can be used on a per subject basis
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
            
            if capflag % this can probably be taken out after system comparison
                
                % the string in the sqare brackets is what it appended to the original file name
                save([fnames{filecurrent}(1,1:end-4),'_bootstrapped_', num2str(STATS.trialcap),'.mat'],'data');
            else
                
                save([fnames{filecurrent}(1,1:end-4),'_bootstrapped.mat'],'data');
            end
            
            waitbar(filecurrent/colfile,h2,sprintf('%12s',[num2str(filecurrent),'/',num2str(colfile)]))
            
        case {'ersp', 'itc'}
            
            %%%% get freq bin length
            %%%% probably a better way of doing this
            %             junk=zeros(1,STATS.numpnts,1);
            %             [ersp] = newtimef(junk, ...
            %                 STATS.numpnts, [STATS.xmin STATS.xmax]*1000, ...
            %                 STATS.srate, STATS.tfcycles,'freqs',STATS.freqs,'timesout',1,'plotersp','off','plotitc','off');
            %             freqbins=length(ersp);
            
            clear ersp
            
            if size(datacell{1},1)>1;
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
                        
                        for i=1:size(datacell{1},1) % num rows/channels
                            
                            
                            if bootcurrent==1
                                % compute TF info to get coeficients for all channels
                                figure('visible','off');
                                [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = newtimef(datacell{1}(i,:,:), ...
                                    STATS.numpnts, [STATS.xmin STATS.xmax]*1000, ...
                                    STATS.srate, STATS.tfcycles,'freqs',STATS.freqs,'nfreqs',STATS.nfreqs,'timesout',-1,'plotersp','off','plotitc','off');
                                    %STATS.srate, STATS.tfcycles,'freqs',STATS.freqs,'timesout',STATS.timesout,'plotersp','off','plotitc','off');
                                freqbins=size(ersp,1);
                                STATS.freqbins=freqbins;
                                STATS.TF_times=times;
                                STATS.TF_freqs=freqs;
                                STATS.timesout=size(tfdata,2);
                                clear ersp itc powbase times freqs erspboot itcboot
                            end
                            
                            
                            if strcmp(STATS.measure, 'ersp')
                                %clear itc powbase times freqs erspboot itcboot
                                
                                mapwrite(trimmean(abs(tfdata(:,:,bootvect)),40,3),[fnames{1,filecurrent}(1,1:end-4), '_', STATS.savestring,'_tempbootstrapped.map'],'datsize',[freqbins,STATS.timesout,size(datacell{1},1)]);
                                
                            elseif strcmp(STATS.measure, 'itc')
                                %clear itc powbase times freqs erspboot itcboot
                                
                                % no trimmed mean?
                                mapwrite(abs(mean(((tfdata(:,:,bootvect,:))./abs(tfdata(:,:,bootvect,:))),4)),[fnames{1,filecurrent}(1,1:end-4), '_', STATS.savestring,'_tempbootstrapped.map'],'datsize',[freqbins,STATS.timesout,size(datacell{1},1)]);
                                
                            end
                            
                        end
                        
                        % this holds multiple channel spectral information
                        datamap=mapread([fnames{1,filecurrent}(1,1:end-4), '_', STATS.savestring, '_tempbootstrapped.map'],'dat');
                        
                        % get rid of previously mapped files
                        %warning off
                        %delete([fnames{1,filecurrent}(1,1:end-4), '_bootstrapped.map']);
                        %warning on
                        
                        % this is the avreage of the multi channel spectral information, like a spectral ROI
                        mapwrite(mean(datamap.Data.dat,3),[fnames{1,filecurrent}(1,1:end-4), '_', STATS.savestring, '_bootstrapped.map'],'datsize',[freqbins,STATS.timesout,STATS.nboot]);
                        
                        % get rid or gathering arrays
                        delete([fnames{1,filecurrent}(1,1:end-4), '_', STATS.savestring, '_tempbootstrapped.map']);
                        
                    case 'single'
                        
                        if bootcurrent==1;
                            
                            % compute TF info to get coeficients
                            figure('visible','off');
                            [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = newtimef(datacell{1}(1,:,:), ...
                                STATS.numpnts, [STATS.xmin STATS.xmax]*1000, ...
                                STATS.srate, STATS.tfcycles,'freqs',STATS.freqs,'nfreqs',STATS.nfreqs,'timesout',-1,'plotersp','off','plotitc','off');
                                    %STATS.srate, STATS.tfcycles,'freqs',STATS.freqs,'timesout',STATS.timesout,'plotersp','off','plotitc','off');
                            freqbins=size(ersp,1);
                            STATS.freqbins=freqbins;
                            STATS.TF_times=times;
                            STATS.TF_freqs=freqs;
                            STATS.timesout=size(tfdata,2);
                            clear ersp itc powbase times freqs erspboot itcboot
                            
                        end
                        
                        
                        if strcmp(STATS.measure, 'ersp')
                            %clear itc powbase times freqs erspboot itcboot
                            
                            % get rid of previously mapped files
                            %warning off
                            %delete([fnames{1,filecurrent}(1,1:end-4), '_bootstrapped.map']);
                            %warning on
                            
                            mapwrite(trimmean(abs(tfdata(:,:,bootvect)),40,3),[fnames{1,filecurrent}(1,1:end-4), '_', STATS.savestring, '_bootstrapped.map'],'datsize',[freqbins,STATS.timesout,STATS.nboot]);
                            
                        elseif strcmp(STATS.measure, 'itc')
                            %clear itc powbase times freqs erspboot itcboot
                            
                            % get rid of previously mapped files
                            %warning off
                            %delete([fnames{1,filecurrent}(1,1:end-4), '_bootstrapped.map']);
                            %warning on
                            
                            % no trimmed mean
                            mapwrite(abs(mean(((tfdata(:,:,bootvect,:))./abs(tfdata(:,:,bootvect,:))),4)),[fnames{1,filecurrent}(1,1:end-4), '_', STATS.savestring, '_bootstrapped.map'],'datsize',[freqbins,STATS.timesout,STATS.nboot]);
                            
                        end
                        
                end
                
                waitbar(bootcurrent/nboot,h3,sprintf('%12s',[num2str(bootcurrent),'/',num2str(nboot)]))
            end
    end
    
    
    %     if capflag % this can probably be taken out after system comparison
    %
    %         % the string in the sqare brackets is what it appended to the original file name
    %         save([fnames{filecurrent}(1,1:end-4),'_bootstrapped_', num2str(STATS.trialcap),'.mat'],'data');
    %     else
    %
    %         save([fnames{filecurrent}(1,1:end-4),'_bootstrapped.mat'],'data');
    %     end
    
    waitbar(filecurrent/colfile,h2,sprintf('%12s',[num2str(filecurrent),'/',num2str(colfile)]))
    
end
close(h2,h3);

disp('******* Finished resmapling from single trials *******')
disp('******* Saving STATS structure *******')
save(['STATS_',STATS.savestring,'.mat'],'STATS');
end

