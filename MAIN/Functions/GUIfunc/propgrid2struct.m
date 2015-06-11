function contextconfig = propgrid2contextconfig(propgrid)

npg=length(propgrid.Properties);
%nlevels=npg/19;

contextconfig=[];
%for li=1:nlevels;
    for pi=1:npg;
        %brkInd=strfind(propgrid.Properties(pi).Name,'[');
        eval(['contextconfig.',propgrid.Properties(pi).Name,'=propgrid.Properties(pi).Value;']);
    end
%end 
