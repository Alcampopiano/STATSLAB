function [funcstr] = varargin_spilt(funcstr, infield )
% when varargin is in a cell array taken from property grid, and each varargin input needs to
% be glued to the end of a string, use this function to spilt varargin up
% and glue it to the end of another string.


% parse varargin inputs from GUI
for i=1:length(infield)
    tmpvar{i}=['tmpvar', num2str(i)];
    varstruct.(tmpvar{i})=infield{i};
    
    if i<length(infield)
        
        funcstr=[funcstr varstruct.(tmpvar{i}) ','];
    else
        funcstr=[funcstr varstruct.(tmpvar{i}) ');'];
    end
    
    
end

%
%     % concatenate varargin things
%     for i=1:length(statslab_propgrid.varargin_resample)
%         if i<length(statslab_propgrid.varargin_resample)
%
%             funcstr=[funcstr varstruct.(tmpvar{i}) ','];
%         else
%             funcstr=[funcstr varstruct.(tmpvar{i}) ');'];
%         end
%     end





end

