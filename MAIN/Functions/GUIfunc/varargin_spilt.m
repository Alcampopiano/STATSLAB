function [funcstr] = varargin_spilt(funcstr, infield )
% Varargin is in a cell array when taken from property grid, so turn it into one string for use in eval.
% Also add appropriate syntax to this string so that users do not have to add matlab syntax when using GUI.

% remove whitespace from ends
infield=strtrim(infield);

% parse varargin inputs from GUI
for i=1:length(infield)
    tmpvar{i}=['tmpvar', num2str(i)];
    
    % determine if its a string or numerical
    
    if any(isnan(str2double_statslab(infield{i}))); % will be NaN unless it is fully numeric
        
        if numel(strsplit(infield{i},' '))>1; %  a cell array of strings, becasue there are spaces
            
            % split it up
            split_cell=strsplit(infield{i},' ');
            
            % turn into a char string
            newstr='{';
            for j=1:length(split_cell);
                
                if j<length(split_cell)
                    
                    newstr=[newstr, '''', split_cell{j}, ''' '];
                else
                    newstr=[newstr, ' ''', split_cell{j}, '''', '}'];
                end
                
            end
            varstruct.(tmpvar{i})=newstr;
            
        else % no spaces so its a full string, but could be empty ('[]')
            
            if strcmp(infield{i}, '[]') || strcmp(infield{i}, '''') % account for empty as a string
                varstruct.(tmpvar{i})=infield{i};
            else
                newstr=['''', infield{i}, '''']; % add single qoutes
                varstruct.(tmpvar{i})=newstr;
                
            end 
        end
        
    else
        
        if i>1
            
            % check varargin for contrast keywords and add transpose operator
            if strcmpi(infield{i-1},'conA') || strcmpi(infield{i-1},'conB') || strcmpi(infield{i-1},'conAB')
                
                % entry should be numeric, so build appropriate string
                % with transpose
                newstr=['[', infield{i}, ']'''];
                varstruct.(tmpvar{i})=newstr;
                
            else
                
                % entry should be numeric, so build appropriate string
                newstr=['[', infield{i}, ']'];
                varstruct.(tmpvar{i})=newstr;
            end
            
        else
            
            % entry should be numeric, so build appropriate string
            newstr=['[', infield{i}, ']'];
            varstruct.(tmpvar{i})=newstr;
            
        end
  
    end
    
    if i<length(infield)
        
        funcstr=[funcstr varstruct.(tmpvar{i}) ','];
    else
        funcstr=[funcstr varstruct.(tmpvar{i}) ');'];
    end
    
    
end
end

