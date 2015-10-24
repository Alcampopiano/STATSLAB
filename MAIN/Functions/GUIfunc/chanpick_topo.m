function chanpick_topo(condfiles_subs,pathtofiles,numconds)

% Basic Graphical User Interface (GUI) without using GUIDE
% Available at https://dadorran.wordpress.com

% There are three basic areas to understand:
%   1. Layout    (how to position objects on the GUI)
%   2. Handles to Objects (used to modify object properties)
%   3. Callback functions (used by User Interface Objects)

%loadhit=[];
%if loadhit==1;
%    [ParamName ParamPath]=uigetfile('*.mat','choose channel selection file:','*.mat','multiselect','off');
%    tmp=load(fullfile(ParamPath, ParamName), '-mat');
%    field=fieldnames(tmp);
%    CHANCHOICES=tmp.(field{1});
%    INLOAD=1;
%    OKAYPUSH=[];
%else
% grid size
tabsize=max(cell2mat(cellfun(@length,condfiles_subs,'un',0)));

chanchoices=cell(tabsize,numconds*2);
j=1;
for i=1:numconds;
    [rowfname colfname]=size(condfiles_subs{i});
    for q=1:rowfname
        chanchoices{q,j}=condfiles_subs{i}{q};
        chanchoices{q,j+1}='';
    end
    j=j+2;
end
%end

%create a figure to house the GUI
f = figure('units','normalized','position',[.25 .25 .5 .65]);

% initialize linked data
data.chanarray=chanchoices;
data.condind=1;
data.subind=1;
guidata(gcf,data);

%create an editable textbox object
edit_box_h = uicontrol('style','edit',...
    'units', 'normalized',...
    'position', [0.45 0.01 0.1 0.05]);

backh = uicontrol('style', 'pushbutton',...
    'string', 'Back',...
    'units', 'normalized',...
    'position', [0.01 0.01 0.1 0.05],...
    'callback', {@backfunc});

nexth = uicontrol('style', 'pushbutton',...
    'string', 'Next',...
    'units', 'normalized',...
    'position', [0.11 0.01 0.1 0.05],...
    'callback', {@nextfunc});

cancelh = uicontrol('style', 'pushbutton',...
    'string', 'Cancel',...
    'units', 'normalized',...
    'position', [0.89 0.01 0.1 0.05],...
    'callback', {@cancelfunc});

okh = uicontrol('style', 'pushbutton',...
    'string', 'OK',...
    'units', 'normalized',...
    'position', [0.79 0.01 0.1 0.05],...
    'callback', {@okayfunc});

loadh = uicontrol('style', 'pushbutton',...
    'string', 'Load',...
    'units', 'normalized',...
    'position', [0.01 0.94 0.1 0.05],...
    'callback', {@loadfunc});


saveh = uicontrol('style', 'pushbutton',...
    'string', 'Save',...
    'units', 'normalized',...
    'position', [0.89 0.94 0.1 0.05],...
    'callback', {@savefunc});


sel=true;
while sel
    
    % invoke data link
    data=guidata(gcf);
    
    cursub=condfiles_subs{data.condind}{data.subind};
    tmpEEG = pop_loadset('filename',cursub,'filepath',pathtofiles{data.condind});
    tmpEEG = eeg_checkset(tmpEEG);
    statslab_topoplot([],tmpEEG.chanlocs, 'style', 'blank', 'drawaxis', 'on', 'electrodes', ...
        'labelpoint', 'plotrad', [], 'chaninfo', tmpEEG, 'nosedir' ,'+Y');
    disp(['working on channel selections for file: ' cursub]);
    
    % waitfor callback
    uiwait;
    
    % gather data linked to the figure
    data=guidata(gcf);
    
    % determine which button was pressed
    if strcmp(data.button, 'next')
        data.button=[];
        if data.subind<length(condfiles_subs{data.condind})
            data.subind=data.subind+1;
        elseif data.subind==length(condfiles_subs{data.condind})
            if data.condind<numconds
                data.condind=data.condind+1;
                data.subind=1;
            end
        end
        
    elseif strcmp(data.button, 'back')
        data.button=[];
        if data.subind>1;
            data.subind=data.subind-1;
        end
        if data.condind>1
            data.condind=cond-1;
        end
        
    elseif strcmp(data.button, 'load')
        
        
    elseif strcmp(data.button, 'save')
        
        
    elseif strcmp(data.button, 'okay')
        
        
    elseif strcmp(data.button, 'cancel')
        
        
        
        
    end
    
    % overwrite linked data
    guidata(gcf,data);

    % remove topoplot objects
    delete(gca);
    delete(gca);
end



%Slider object to control ellipse size
% uicontrol('style','Slider',...
%             'Min',0.5,'Max',2,'Value',1,...
%             'units','normalized',...
%             'position',[0.1    0.2    0.08    0.25],...
%             'callback',{@change_size,ellipse_h,ellipse_position });
%
% uicontrol('Style','text',...
%             'units','normalized',...
%             'position',[0    0.45    0.2    0.1],...
%             'String','Ellipse Size')

end



%%
function savechans(object_handle,event,chanchoices)
[ParamName, ParamPath]=uiputfile('*.*','channel selection file');
save(fullfile(ParamPath, ParamName),'chanchoices');
end


%%
function loadchans(object_handle,event,chanchoices)


end

%%
function okayfunc(object_handle,event,chanchoices)


end


%%
function cancelfunc(object_handle,event,chanchoices)


end


%%
function nextfunc(object_handle,event)

h=findobj(gcf,'Color','g');
labs=get(h,'String');

% invoke data link
data=guidata(gcf);

if isempty(data.chanarray{data.subind,data.condind*2}) % no channels have been selected
    data.chanarray{data.subind,data.condind*2}=labs;
    
else % some chans have been selected previously
    data.chanarray{data.subind,data.condind*2}=[data.chanarray{data.subind,data.condind*2}; labs];
    
end

data.button='next';
guidata(gcf,data);
uiresume(gcbf)
end

%%
function backfunc(object_handle,event)

h=findobj(gcf,'Color','g');
labs=get(h,'String');

% invoke data link
data=guidata(gcf);

if isempty(data.chanarray{data.subind,data.condind*2}) % no channels have been selected
    data.chanarray{data.subind,data.condind*2}=labs;
    
else % some chans have been selected previously
    data.chanarray{data.subind,data.condind*2}=[data.chanarray{data.subind,data.condind*2}; labs];
    
end

data.button='back';
guidata(gcf,data);
uiresume(gcbf)
end


%updated eg_fun used to demonstrate passing  iables
%copy paste this code into a file called eg_fun.m
% function eg_fun(object_handle, event, edit_handle, ellipse_handle)
%     str_entered = get(edit_handle, 'string');
%
%     if strcmp(str_entered, 'red')
%         col_val = [1 0 0];
%     elseif strcmp(str_entered, 'green')
%         col_val = [0 1 0];
%     elseif strcmp(str_entered, 'blue')
%          col_val = [0 0 1];
%     else
%         col_val = [0 0  0];
%     end
%     set(ellipse_handle, 'facecolor', col_val);
%
% %change_size code --------------------------------------------------
% %copy paste this code into a file called change_size.m
% function  change_size(objHandel, evt, annotation_handle, orig_pos)
%     slider_value = get(objHandel,'Value');
%     new_pos = orig_pos;
%     new_pos(3:4) = orig_pos(3:4)*slider_value;
%     set(annotation_handle, 'position', new_pos)
