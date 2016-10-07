function buildsvg_subject(STATS,options,test_stat,CI_diffup,sigvect,savename)

fh=figure; imagesc(flipud(STATS.xtimes(options.timeplot)),options.zaxis,test_stat)
ylim(options.zaxis);
set(gca,'YDir','normal');
caxis(options.yaxis);

hold on
jbfill(STATS.xtimes(options.timeplot),zeros(1,length(test_stat))+options.zaxis(2),CI_diffup,[1 1 1],[1 1 1],1,1);
hold on


for i=1:length(sigvect);
    if sigvect(i)==1;
        h=jbfill(STATS.xtimes(options.timeplot(i)),options.zaxis(2)/20,0,[.4 .4 .4],[.4 .4 .4],1,1);
        set(h,'LineWidth',3)
    end
end

plot2svg(savename);
close(fh);

% %get cdata
% data=get(ah,'Cdata');
% 
% % define meshgrid
% [X,Y] = meshgrid(1:size(data,2), 1:size(data,1));
% 
% %// Define a finer grid of points
% [X2,Y2] = meshgrid(1:0.05:size(data,2), 1:0.05:size(data,1));
% 
% 
% if strcmp(maptype,'jet')
%     %// Interpolate the data and show the output
%     outData = interp2(X, Y, data, X2, Y2, 'linear');
%     
% elseif strcmp(maptype,'bone') || strcmp(maptype,'bonesub')
%     
%     outData=data;
% end
% 
% % make image of the data
% ifh=figure;
% imagesc(flipud(STATS.xtimes(timeplot)), outData);
% set(gca,'YDir','normal')
% caxis(axlim);
% 
% if strcmp(maptype,'jet')
%     cbfreeze(colorbar);
%     %freezeColors;
% elseif strcmp(maptype,'bone')
%     colormap(ifh,bone(2));
%     hb=cbfreeze(colorbar);
%     set(hb,'YTick',[0 1]);
%     
% elseif strcmp(maptype,'bonesub')
%     colormap(ifh,flipud(hot));
%      hb=cbfreeze(colorbar); 
% end

end