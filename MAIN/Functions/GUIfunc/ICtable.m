function ICtable
f = gcf;
global ICCHOICES
t = uitable('Parent',f,...
    'units', 'normalized', ...
    'Position', [.05 .10 .9 .80], ...  
    'data',ICCHOICES,...
    'ColumnEditable', true, ...
    'CellEditCallback', 'global ICCHOICES; ICCHOICES = get(gco, ''Data'');');
end

