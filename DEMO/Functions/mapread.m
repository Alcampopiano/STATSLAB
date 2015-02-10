function map = mapread(fname,fieldname,varargin)

g = struct(varargin{:});

try g.datsize;    catch, g.datsize   = [];         end;
try g.offsetadj;  catch, g.offsetadj = 0;          end;

% Creat file identifier...
fid_m=fopen(fname,'r');

% Read header...
headlen=fread(fid_m,1,'int16');
datsize=fread(fid_m,headlen-1,'int16')';

fclose(fid_m);

if ~isempty(g.datsize);
    datsize=g.datsize;
end

% Map data
map = memmapfile(fname,'Offset',(headlen*2)+g.offsetadj*4,'Format',{'single',datsize,fieldname},'Repeat',1);