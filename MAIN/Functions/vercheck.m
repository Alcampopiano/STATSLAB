function vercheck()
% check version


% load file with version number
if isunix
    verfile=load('/STATSLAB/ver.mat');
    slashfind=@(x) strfind(x,'/');
else ispc
    verfile=load('\STATSLAB\ver.mat');
    slashfind=@(x) strfind(x,'\');
end

curver=verfile.ver;

% query git hub tag and compare to current ver number
[jnk, gitinfo]=system('git ls-remote --tags git@github.com:Alcampopiano/STATSLAB.git');

% indentify unique tag zones
slashind=slashfind(gitinfo);

% find relavent slashes
slashind=slashind(2:4:end);


for i=1:length(slashind);
    tag(i)=str2num(gitinfo(slashind(i)+1:slashind(i)+4));
end

tag=max(tag);

if tag==curver;
    disp(['You are using the newest version of STATSLAB -- v',num2str(curver)]);
    
elseif curver<tag;
    disp('A newer version of STATSLAB is available');
    disp(['To upgrade from v',num2str(curver),' to v',num2str(tag),' visit <a href="https://github.com/Alcampopiano/STATSLAB'',''-browser', '">https://github.com/Alcampopiano/STATSLAB</a>']);  

end