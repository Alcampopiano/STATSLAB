function [okayhit, chan_choices]=chanpick_topo(condfiles_subs,pathtofiles,numconds,varargin)

clearvars -global tmpEEG CHANCHOICES CURCLICK SAVEHIT LOADHIT

global CURCLICK tmpEEG CHANCHOICES SAVEHIT LOADHIT



% the first part of evstring1 just deletes a weird uicontrol object that supergui places
% on the figure that appears in blue. If this is not the 6th object in the
% figure, this code will not work. Better to find the object by name
% somehow.
evstring1='b=findobj(gcf); delete(b(6)); set(gcf,''Color'', [1 1 1]); statslab_topoplot_loadfile;';


if isempty(varargin{1})
    CHANCHOICES=cell(tabsize,numconds*2);
    j=1;
    for i=1:numconds;
        [rowfname colfname]=size(condfiles_subs{i});
        for q=1:rowfname
            CHANCHOICES{q,j}=condfiles_subs{i}{q};
            CHANCHOICES{q,j+1}='';
        end
        j=j+2;
    end
    
    
else
    
    % load CHANCHOICES
    [ParamName, ParamPath] = uigetfile('*.mat','choose channel selection file:','*.mat','multiselect','off');
    load(fullfile(ParamPath, ParamName), '-mat');
    
end


for jj=1:numconds;
    for ii=1:length(condfiles_subs{jj}); % number of subs in each condition,
        
        if ~isempty(varargin{1})
            tmpEEG = pop_loadset('filename',condfiles_subs{jj}{ii},'filepath',pathtofiles{1});
            tmpEEG = eeg_checkset(tmpEEG);
            
            
            if ~isempty(CHANCHOICES{ii,jj*2});
                for q=1:length(CHANCHOICES{ii,jj*2});
                    labs=tmpEEG.chanlocs(CHANCHOICES{ii,jj*2}{q}).labels;
                    
                    oh=findobj('String',[labs, '_']);
                    
                    if isempty(oh)
                        oh=findobj('String',labs);
                    end
                    set(oh, 'Color', 'green', 'FontSize',13, 'FontWeight','bold');
                end
                
            end
            LOADHIT=0;
            
        else
            tmpEEG = pop_loadset('filename',condfiles_subs{jj}{ii},'filepath',pathtofiles{1});
            tmpEEG = eeg_checkset(tmpEEG);
            disp(['working on channel selections for file: ' condfiles_subs{jj}{ii}]);
            
            [res, jnk, okayhit]=inputgui( ...
                'geom', ...
                {...
                {6 22 [0 0] [1 1]} ... % this is just to control the size of the GUI {6 22 [0 0] [1 1]} ...
                {6 22 [0.05 0] [2 1]} ... % load config
                {6 22 [3.95 0] [2 1]} ... %2 saveas
                
                }, ...
                'uilist', ...
                {...
                {'Style', 'text', 'tag','txt_ccfp','string',blanks(30)} ... %1 this is just to control the size of the GU
                {'Style', 'pushbutton', 'string', 'Load Parameters', ...
                'callback', 'global LOADHIT; LOADHIT=1; close(gcf)'} ... %2
                {'Style', 'pushbutton','string','Save channel selection file', ...
                'callback','global SAVEHIT; SAVEHIT=1;'}, ...
                }, ...
                'title', 'STATSLAB -- statslab()',...%, ...
                'eval',evstring1 ...
                );
            
            
            
            
            
            
            
        end
        
        if SAVEHIT==1;
            CHANCHOICES{ii,jj*2}=CURCLICK;
            [ParamName, ParamPath]=uiputfile('*.*','channel selection file');
            save(fullfile(ParamPath, ParamName),'CHANCHOICES');
            SAVEHIT=0;
            okayhit='hello';
            
        elseif LOADHIT==1;
            
            [okayhit, chan_choices]=chanpick_topo(condfiles_subs,pathtofiles,numconds,1);
            
            % load CHANCHOICES
            [ParamName, ParamPath] = uigetfile('*.mat','choose channel selection file:','*.mat','multiselect','off');
            load(fullfile(ParamPath, ParamName), '-mat');
            %LOADHIT=0;
            okayhit='hello';
            
            
        else
            CHANCHOICES{ii,jj*2}=CURCLICK;
        end
        
        % is cancel is hit or fig is closed
        if isempty(okayhit)
            chan_choices=CHANCHOICES;
            %clearvars -global tmpEEG CURCLICK CHANCHOICES SAVEHIT
            return
        end
    end
end

chan_choices=CHANCHOICES;
clearvars -global tmpEEG CHANCHOICES CURCLICK SAVEHIT LOADHIT
end