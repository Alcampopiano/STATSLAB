function fid_m = mapwrite(data,fname,varargin)

g = struct(varargin{:});

try g.datsize; 		catch, g.datsize	 = [];                    end;
%try g.fid; 		catch, g.fid	 = [];                    end;

%Create file identifier...
if exist(fname,'file') %If the file already exist...
    %Open the file and give permission to append data...
    fid_m=fopen(fname,'a+');
else %If the file does not already exist...
    %Create new file...
    fid_m=fopen(fname,'w');
    
    %Write header...
    if isempty(g.datsize)
        headlen=ndims(data)+1;
        header=[headlen,size(data)];
    else
        headlen=length(g.datsize)+1;
        header=[headlen,g.datsize];
    end
    fwrite(fid_m,header,'int16');

end

% Write data...
fwrite(fid_m,data,'single');

% Close files...
fclose(fid_m);
