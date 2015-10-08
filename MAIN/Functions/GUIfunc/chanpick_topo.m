function [okayhit, chan_choices]=chanpick_topo(condfiles_subs,pathtofiles,numconds)

%clearvars -global CURCLICK tmpEEG CHANCHOICES SAVEHIT CURSUB LOADHIT CONDFILES_SUBS PATHTOFILES NUMCONDS NEXT BACK

global CURCLICK tmpEEG CHANCHOICES SAVEHIT CURSUB LOADHIT CONDFILES_SUBS PATHTOFILES NUMCONDS NEXT BACK INLOAD OKAYPUSH CURPOS

% make global so eval can use them in inputgui
CONDFILES_SUBS=condfiles_subs;
PATHTOFILES=pathtofiles;
NUMCONDS=numconds;

% deletes a weird uicontrol object that supergui places on the figure that appears in blue.
evstring1='global CURPOS; CURPOS=get(gcf, ''Position''); statslab_topoplot_loadfile; b=findobj(gcf); delete(b(9)); set(gcf,''Color'', [1 1 1]);';

% recall function
evstring2=['global CURPOS LOADHIT CONDFILES_SUBS PATHTOFILES NUMCONDS; CURPOS=get(gcf, ''Position''); LOADHIT=1;' ...
    'close(gcf); chanpick_topo(CONDFILES_SUBS,PATHTOFILES,NUMCONDS);'];

if LOADHIT==1;
    [ParamName ParamPath]=uigetfile('*.mat','choose channel selection file:','*.mat','multiselect','off');
    tmp=load(fullfile(ParamPath, ParamName), '-mat');
    field=fieldnames(tmp);
    CHANCHOICES=tmp.(field{1});
    INLOAD=1;
    OKAYPUSH=[];
else
    % grid size
    tabsize=max(cell2mat(cellfun(@length,condfiles_subs,'un',0)));
    
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
end

jj=1;
ii=1;
sel=true;

screenpos='cent';
try
    while sel
        
        CURSUB=condfiles_subs{jj}{ii};
        
        tmpEEG = pop_loadset('filename',CURSUB,'filepath',pathtofiles{jj});
        tmpEEG = eeg_checkset(tmpEEG);
        
        disp(['working on channel selections for file: ' CURSUB]);
        
        if ~isempty(CURPOS)
            screenpos=[CURPOS(1) CURPOS(2)];
        end
        
        [res, jnk, okayhit]=inputgui( ...
            'geom', ...
            {...
            {6 22 [0 0] [1 1]} ... % this is just to control the size of the GUI {6 22 [0 0] [1 1]} ...
            {6 22 [0.05 0] [2 1]} ... % load config
            {6 22 [3.95 0] [2 1]} ... %2 saveas
            {6 22 [0.05 95] [1 1]} ... %3 BACK
            {6 22 [0.85 95] [1 1]} ... %3 NEXT
            
            }, ...
            'uilist', ...
            {...
            {'Style', 'text', 'tag','txt_ccfp','string',blanks(30)} ... %1 this is just to control the size of the GU
            {'Style', 'pushbutton', 'string', 'Load Parameters', ...
            'callback', evstring2} ... %2
            {'Style', 'pushbutton','string','Save channel selection file', ...
            'callback','global SAVEHIT CURPOS; CURPOS=get(gcf, ''Position''); SAVEHIT=1; close(gcf);'}, ...
            {'Style', 'pushbutton','string','BACK', ...
            'callback','global BACK CURPOS; CURPOS=get(gcf, ''Position''); BACK=1; close(gcf);'}, ...
            {'Style', 'pushbutton','string','NEXT', ...
            'callback','global NEXT CURPOS; CURPOS=get(gcf, ''Position''); NEXT=1; close(gcf);'}, ...
            }, ...
            'title', 'STATSLAB -- statslab()',...%, ...
            'eval',evstring1, ...
            'screenpos', screenpos ...
            );
        
        if SAVEHIT==1;
            
            if ~isempty(CHANCHOICES{ii,jj+1})
                CHANCHOICES{ii,jj+1}=[CHANCHOICES{ii,jj+1} CURCLICK];
            else
                CHANCHOICES{ii,jj+1}=CURCLICK;
            end
            [ParamName, ParamPath]=uiputfile('*.*','channel selection file');
            save(fullfile(ParamPath, ParamName),'CHANCHOICES');
            SAVEHIT=0;
            okayhit='cont';
        else
            
            if ~isempty(CHANCHOICES{ii,jj+1})
                CHANCHOICES{ii,jj+1}=[CHANCHOICES{ii,jj+1} CURCLICK];
            else
                CHANCHOICES{ii,jj+1}=CURCLICK;
            end
        end
        
        if isempty(BACK) && isempty(NEXT);
            
            % control while loop params
            if ii<length(condfiles_subs{jj})
                ii=ii+1;
                
            elseif ii==length(condfiles_subs{jj});
                ii=1;
                
                if jj<numconds;
                    jj=jj+1;
                    
                elseif jj==numconds
                    break
                end
            end
            
        elseif BACK==1;
            BACK=[];
            if ii>1;
                ii=ii-1;
            end
            if jj>1
                jj=jj-1;
            end
               
            okayhit='cont';
        elseif NEXT==1;
            NEXT=[];
            if ii<length(condfiles_subs{jj})
                ii=ii+1;
            elseif ii==length(condfiles_subs{jj})
                if jj<numconds
                    jj=jj+1;
                    ii=1;
                end
            end
            okayhit='cont';
        end
        
        % is cancel is hit or fig is closed
        if strcmp(okayhit,'retuninginputui') && ~isempty(INLOAD);
                INLOAD=[];
                OKAYPUSH=1;
                return
        elseif  isempty(okayhit) && ~isempty(INLOAD);
                INLOAD=[];
                return
        elseif strcmp(okayhit,'retuninginputui') && isempty(INLOAD);        
                break
        elseif isempty(okayhit) && ~isempty(OKAYPUSH);
                break
        elseif isempty(okayhit) && isempty(OKAYPUSH);
                break
        end    
    end
catch
    disp('error catch')
    clearvars -global CURCLICK tmpEEG CHANCHOICES SAVEHIT CURSUB LOADHIT CONDFILES_SUBS PATHTOFILES NUMCONDS NEXT BACK INLOAD OKAYPUSH CURPOS
end


% ending channel selection
if isempty(OKAYPUSH)
    
    chan_choices=[];
    clearvars -global CURCLICK tmpEEG CHANCHOICES SAVEHIT CURSUB LOADHIT CONDFILES_SUBS PATHTOFILES NUMCONDS NEXT BACK INLOAD OKAYPUSH CURPOS
    
elseif ~isempty(OKAYPUSH) % okay was hit
    chan_choices=CHANCHOICES;
    clearvars -global CURCLICK tmpEEG CHANCHOICES SAVEHIT CURSUB LOADHIT CONDFILES_SUBS PATHTOFILES NUMCONDS NEXT BACK INLOAD OKAYPUSH CURPOS
end

end