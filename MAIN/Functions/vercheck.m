function vercheck()
% check version

try
    % load file with version number
    if isunix && ~ismac
        verfile=load('/MAIN/Functions/statslabver/ver.mat');
        pathstr='Alcampopiano/STATSLAB/releases/tag/';
        pathfind=@(x) strfind(x,pathstr);
        evalstr='content=$(wget https://github.com/Alcampopiano/STATSLAB/releases -q -O -); echo $content';
        
    elseif ismac
        
        verfile=load('/MAIN/Functions/statslabver/ver.mat');
        pathstr='Alcampopiano/STATSLAB/releases/tag/';
        pathfind=@(x) strfind(x,pathstr);
        evalstr='content=$(curl -L https://github.com/Alcampopiano/STATSLAB/releases); echo $content';
        
    elseif ispc
       %disp('go here to download wget for windows, which I need to perform automatic checks for latest STATSLAB versions');
        pathstr='Alcampopiano/STATSLAB/releases/tag/';
        verfile=load('MAIN\Functions\statslabver\ver.mat');
        evalstr='powershell $webreq = Invoke-WebRequest https://github.com/Alcampopiano/STATSLAB/releases; echo $webreq.content';
        pathfind=@(x) strfind(x,pathstr);
    end
    
    curver=verfile.ver;
    
    % query git hub tag and compare to current ver number
    [jnk, gitinfo]=system(evalstr);
   
    % indentify unique tag zones
    pathinds=pathfind(gitinfo);
    pathlength=length(pathstr);
    
    
    for i=1:length(pathinds);
        
        temptag=gitinfo(pathinds(i)+pathlength:pathinds(i)+pathlength+3);
        
        % remove ending quote or character that wont be converted to number
        if isempty(str2num(temptag));
            
            % couldnt be converted, trim trailing character
            temptag=temptag(1:end-1);
            temptag=str2num(temptag);
            
        else
            temptag=str2num(temptag);
        end
        
        % gather tags
        tag(i)=temptag;
        
    end
    
    tag=max(tag);
    
    if tag==curver;
        disp(['You are using the newest version of STATSLAB -- v',num2str(curver)]);
        
    elseif curver<tag;
        
        disp(['A newer version of STATSLAB is available. To upgrade from v', num2str(curver),' to v',num2str(tag), ',']);
        
        if isunix && ~ismac
            disp('<a href="matlab: ! sudo gnome-open https://github.com/Alcampopiano/STATSLAB">https://github.com/Alcampopiano/STATSLAB</a>');
            
            % the below are alternatives to using gnome
            %sudo sensible-browser http://www.google.com
            %sudo xdg-open http://www.google.com')
            
            % the below might work on some linux machines
            %disp([' visit <a href="https://github.com/Alcampopiano/STATSLAB'',''-browser', '">https://github.com/Alcampopiano/STATSLAB</a>']);
            
        elseif ismac
            disp([' visit <a href="https://github.com/Alcampopiano/STATSLAB'',''-browser', '">https://github.com/Alcampopiano/STATSLAB</a>']);
        elseif ispc
            disp([' visit <a href="https://github.com/Alcampopiano/STATSLAB'',''-browser', '">https://github.com/Alcampopiano/STATSLAB</a>']);
            
        end
        
    end
    
catch
    disp('could not do automatic check for newer versions of STATSLAB.')
end





