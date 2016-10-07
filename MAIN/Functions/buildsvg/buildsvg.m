function buildsvg(ah,axlim,STATS,savename,maptype,timeplot)

%get cdata
data=get(ah,'Cdata');

% define meshgrid
[X,Y] = meshgrid(1:size(data,2), 1:size(data,1));

%// Define a finer grid of points
[X2,Y2] = meshgrid(1:0.05:size(data,2), 1:0.05:size(data,1));


if strcmp(maptype,'jet')
    %// Interpolate the data and show the output
    outData = interp2(X, Y, data, X2, Y2, 'linear');
    
elseif strcmp(maptype,'bone') || strcmp(maptype,'bonesub')
    
    outData=data;
end

% make image of the data
ifh=figure;
imagesc(flipud(STATS.TF_times(timeplot)),STATS.TF_freqs, outData);
set(gca,'YDir','normal')
caxis(axlim);

if strcmp(maptype,'jet')
    cbfreeze(colorbar);
    %freezeColors;
elseif strcmp(maptype,'bone')
    colormap(ifh,bone(2));
    hb=cbfreeze(colorbar);
    set(hb,'YTick',[0 1]);
    
elseif strcmp(maptype,'bonesub')
    colormap(ifh,flipud(hot));
     hb=cbfreeze(colorbar); 
end

plot2svg(savename);
close(ifh);

end