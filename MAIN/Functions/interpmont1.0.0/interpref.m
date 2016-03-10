%This function reads coordinates files and interpolates
%data to these locations from the currently loaded data in EEGLAB. The
%average of the interpolated channels is then calculated and subtracted
%from each of channels in the currently loaded data in EEGLAB.

%It is important that the sfp file containing the reference channel names
%and coordinates 1- contains the same fiducials as the loaded coordinates
%for the data field and 2- the channel names of the reference sites do not
%match any of the names from the data channels.

function EEG = interpref(EEG,coord_fname,varargin)

% INITIATE VARARGIN STRUCTURES...
try
    options = varargin;
    for index = 1:length(options)
        if iscell(options{index}) && ~iscell(options{index}{1}), options{index} = { options{index} }; end;
    end;
    if ~isempty( varargin ), g=struct(options{:});
    else g= []; end;
catch
    disp('ce_eegplot() error: calling convention {''key'', value, ... } error'); return;
end;

try g.chans; catch, g.chans=1:EEG.nbchan; end

tmpEEG=EEG;

EEG=pop_select(EEG,'channel',g.chans);
n_chans=EEG.nbchan;

reflocs=readlocs(coord_fname);
n_reflocs=length(reflocs);
for i=1:n_reflocs;
    reflocs(i).labels=['ref',reflocs(i).labels];
end

tmplocs=EEG.chanlocs;

EEG.chanlocs=[];
for i=1:n_chans;
    EEG.chanlocs(i).labels=tmplocs(i).labels;
    EEG.chanlocs(i).X=tmplocs(i).X;
    EEG.chanlocs(i).Y=tmplocs(i).Y;
    EEG.chanlocs(i).Z=tmplocs(i).Z;
end

for i=1:n_reflocs;
    EEG.data(n_chans+i,:,:)=zeros(size(EEG.data(1,:,:)));
    EEG.chanlocs(n_chans+i).labels=reflocs(i).labels;
    EEG.chanlocs(n_chans+i).X=reflocs(i).X;
    EEG.chanlocs(n_chans+i).Y=reflocs(i).Y;
    EEG.chanlocs(n_chans+i).Z=reflocs(i).Z;
end

%EEG.nbchan=length(EEG.chanlocs);
EEG.chanlocs = convertlocs(EEG.chanlocs,'cart2all');

%interoplate channels created from coord file.
EEG = eeg_interp(EEG,n_chans+1:n_chans+n_reflocs,'spherical');

%create the rereference channel from the interpolated coordfile channels.
RefChan=mean(EEG.data(n_chans+1:n_chans+n_reflocs,:,:),1);

%Restore original EEG structure...
EEG=tmpEEG;

%subtract the RefChan rom each of the channels in the data field.
for i=1:EEG.nbchan;
    EEG.data(i,:,:)=EEG.data(i,:,:)-RefChan;
end
