function lowess_kd(XdataIn,YdataIn,nboot,span,nbins)


if isempty(XdataIn) && isempty(YdataIn)
    
    % select and load
    xvar=uigetfile('mat','Select X variable');
    yvar=uigetfile('mat','Select Y variable');
    x=load(xvar);
    y=load(yvar);
    
    % akward but effective way of getting variables out of structure
    xfld=fieldnames(x);
    x=x.(xfld{1});
    yfld=fieldnames(y);
    y=y.(yfld{1});
    
else
    x=XdataIn;
    y=YdataIn;  
end


% preallocate
xdata=zeros(1000,nboot);
ydata=zeros(1000,nboot);
kd_array=zeros(nbins+1,2^14);
kd_inds=zeros(nbins+1,2^14);
xspag=zeros(nboot,length(x));
yspag=zeros(nboot,length(x));

% preallocate new grid system for single surface plot
%ygrid=linspace(min(y),max(y),100);
%[colygrid]=size(ygrid,2);
%gridsystem=zeros(colygrid,nbins+1);

% there are many curve types
curve_type='lowess';

%{
% smoother on original data
hold on
smoother = smooth(x,y,span,curve_type);
[xsort,indx] = sort(x);
plot(xsort,smoother(indx),'Color','r','LineWidth',2);
hold on
hp=plot(x,y,'ko','MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize', 5);
%}

% begin bootstrap loop
for i=1:nboot;
    
    % draw a random sample
    bootvect=randi(length(x), 1,length(x));
    xboot=x(bootvect);
    yboot=y(bootvect);
    
    % calculate smoother, sort x and get indices
    smoother_boot = smoother(xboot,yboot,span,curve_type);
    [xboot_sort,xboot_inds] = sort(xboot);
    smoother_boot_sort=smoother_boot(xboot_inds);
    gather_xy_boot=[xboot_sort smoother_boot_sort];
    
    % gather up the original data, in case you want a spagetti plot
    xspag(i,:)=xboot_sort;
    yspag(i,:)=smoother_boot(xboot_inds);
    
    % interpolate
    distance = sqrt(sum(diff(gather_xy_boot,1,1).^2,2));  % Distance between subsequent points
    paracoord = [0; cumsum(distance)];               % Parametric coordinate
    coordnew = linspace(0,paracoord(end),1000).';   % 1000 evenly spaced points from 0 to s(end)
    x_new = interp1q(paracoord,gather_xy_boot(:,1),coordnew);     % Interpolate new x values
    y_new = interp1q(paracoord,gather_xy_boot(:,2),coordnew);     % Interpolate new y values
    %[coordall,sortindex] = sort([paracoord; coordnew]);  % Sort all the parametric coordinates
    %xcoords = [gather_xy_boot(:,1); x_new];   % Collect the x coordinates
    %xcoords = xcoords(sortindex);              % Sort the x coordinates
    %ycoords = [gather_xy_boot(:,2); y_new];   % Collect the y coordinate
    %ycoord = ycoords(sortindex);              % Sort the y coordinates
    xdata(:,i)=x_new;
    ydata(:,i)=y_new;
    
end

% bin data
botEdge=min(min(xdata));
topEdge=max(max(xdata));
binEdges = linspace(botEdge, topEdge, nbins+1);
[bincounts bininds] = histc(xdata,binEdges);

% create bin fields
[rowbin colbin]=size(bincounts);
[rowbinind colbinind]=size(bininds);
sep_bin_names=cell(rowbin,1);

% create bin fields
for i=1:rowbin;
    sep_bin_names{i}=['bin',num2str(i)];
end

% populate bins
for i=1:rowbin;
    tmp=bininds==i;
    %bin_namesX.(sep_bin_names{i})=xdata(tmp);
    bin_namesY.(sep_bin_names{i})=ydata(tmp);
    clear tmp
end

%{
  MATLAB native kernal density
% estimate kernal density based on binned data
for i=1:rowbin;
    [f,xi]=ksdensity(bin_namesY.(sep_bin_names{i}));
    kd_array(i,:)=f;
    kd_inds(i,:)=xi;
end
%}

% fill up weird emtpy bins that sometimes appear at the final bins
% TEST THIS
for i=1:rowbin;
    if isempty(bin_namesY.(sep_bin_names{i})) || size(bin_namesY.(sep_bin_names{i}),1)==1;
        bin_namesY.(sep_bin_names{i})=bin_namesY.(sep_bin_names{i-1});
    end
end

% BOTEV's, which is better than MATLAB's
% estimate kernal density based on binned data
for i=1:rowbin;
    [bandwidth,f,xi,cdf]=kde(bin_namesY.(sep_bin_names{i}),2^14);
    kd_array(i,:)=f;
    kd_inds(i,:)=xi;
end

% calculate confidence intervals



% plot options


% FIX THE INVERTED FIGURE THAT COMES FROM THIS
% I'm not even sure that you need to force a single surface/resoloution for
% the data to sit on. Anyway, it may be quicker, and smoother looking, but
% I'm not sure. It seemed to look very similar to the original way of doing
% things when I first tried it (aside from being inverted).
% Also, the backround has to be dealt with in a different way potentially
% becasue the vales in the background all take zero and are thus colored as
% the min color on the current color map. A mask could deal with this if
% you wanted the background to be another color. The non-downsampled data
% doesn't suffer from this as the data is independant of the figure's
% backround. Anyway, the below chunck of code can be ignored if you wish.
% Separate plot options are below for downsampled and non-downsampled data.

% for i=1:rowbin;
%     
%     % find nearest ygrid vals
%     lowval = min(kd_inds(i,:)); %value to find
%     highval= max(kd_inds(i,:));
%     tmplow = abs(ygrid-lowval);
%     tmphigh = abs(ygrid-highval);
%     [jnk idxlow] = min(tmplow); %index of closest value
%     [jnk idxhigh] = min(tmphigh); %index of closest value
%     new_length=idxhigh-idxlow+1;
%     disp(['newlen=', num2str(new_length)]);
%     
%     % build single surface by downsampling
%     dec_factor=(2^14)/new_length;
%     down_samp_array=decimate(kd_array(i,:),round(dec_factor))'; % DIFFERENT ON A MAC, IT IS A COLUMN VECTOR THERE.
%     [down_samp_col]=length(down_samp_array); % SO THIS NEEDS TO BE GETTING THE LENGTH
%     disp(['downsamp cols=', num2str(down_samp_col)]);
%     
%     from=colygrid-((idxlow-1)+(down_samp_col-1));
%     to=colygrid-(idxlow-1);
%     if from<1; disp(num2str(from)); trm=abs(from); from=1; down_samp_array(end-trm:end)=[]; disp(['from<1 i=', num2str(i)]);end
%     %if to>length(ygrid);  trm=to-(length(ygrid)); down_samp_array(end-trm:end)=[]; disp(['to>ygrid leng i=', num2str(i)]); end
%     gridsystem(from:to,i)=flipud(down_samp_array);
%     
%    %gridsystem2(idxlow:idxhigh,i)=down_samp_array;
%     
% end


% xc=linspace(min(x), max(x), size(gridsystem,2));
% yc=linspace(min(y), max(y), size(gridsystem,1));
% % CODE to plot on downsampled surface
% figure;surf(xc,yc,flipud(gridsystem), 'LineStyle','none'); shading interp;
% axis tight
% view(0,90)
% colormap(flipud(gray))
% 
% % plot data on top
% zz=max(max(kd_array));
% zz=zeros(1,size(x,1))+zz;
% hold on
% 
% % smoother on original data
% hold on
% smootherline = smoother(x,y,span,curve_type);
% [xsort,indx] = sort(x);
% plot3(xsort,smootherline(indx),zz,'Color',[1 1 1],'LineWidth',1);
% % %hold on
% % 
% % 
% % set data color and other properties
% hp=plot3(x,y,zz,'bo','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize', 5);
% hold on;
% %grid on
% axis tight



%{
% plot spagetti
% figure;
col=summer(nboot);
for i=1:nboot;
    plot(xspag(i,:),yspag(i,:),'Color',col(i,:))
    hold on
end
axis tight
grid on
%}

% if you want to change the backround of the figure (inside the axes)
% use...set(gca,'Color',[0 0 0]);


% gather indices which limit the kernal density measure
kdiCI=cell(1,rowbin);
kdaCI=cell(1,rowbin);

for i=1:rowbin;
    miny=min(bin_namesY.(sep_bin_names{i}));
    maxy=max(bin_namesY.(sep_bin_names{i}));
    [V minind]=min(abs(miny-kd_inds(i,:)));
    [V maxind]=min(abs(maxy-kd_inds(i,:)));
    
    %kd_inds(i,minind:maxind);
    kdi=kd_inds(i,:);
    kda=kd_array(i,:);
    mii=kdi<kdi(minind);
    mai=kdi>kdi(maxind);
    
    kdi(1,mii)=nan;
    kdi(1,mai)=nan;
    kda(1,mii)=nan;
    kda(1,mai)=nan;
    
    in=isnan(kdi);
    kdi(in)=[];
    kda(in)=[];
    
    
    % calculate CIs
    len=length(kdi);
    low=round(.05*len/2);
    up=len-low;
    kdi=kdi(low+1:up); % ycoord
    kda=kda(low+1:up); % zcoord
    kdiCI{i}=kdi;
    kdaCI{i}=kda;

end

% GENERAL PLOT CODE for data not downsampled
% will likely casue overhang in the figure
% so tidy up with ...xlim([9 50]) or whatever, query with -> xlim
%hold on
% figure;
% set(gca,'Color',[1 1 1]);
% for i=2:(size(kd_array,1)-2);
%     ycoord=kd_inds(i,:);
%     xcoord=zeros(1,size(ycoord,2))+binEdges(i);
%     zcoord=kd_array(i,:);
%     col = zcoord;  % This is the color, vary with x in this case.
%     surface([xcoord;xcoord],[ycoord;ycoord],[zcoord;zcoord],[col;col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2);
% end
% 
% axis tight
% colormap(flipud(gray))

figure;
set(gca,'Color',[1 1 1]);
for i=2:(size(kd_array,1)-2);
    
    ycoord=kdiCI{i};
    xcoord=zeros(1,size(ycoord,2))+binEdges(i);
    
    %zcoord=kd_array(i,:);
    zcoord=kdaCI{i};

    col = zcoord;  % This is the color, vary with x in this case.
    surface([xcoord;xcoord],[ycoord;ycoord],[zcoord;zcoord],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
end
axis tight
colormap(flipud(gray))

% plot data on top
zz=max(max(kd_array));
zz=zeros(1,size(x,1))+zz;
hold on

% smoother on original data
hold on
smootherline = smoother(x,y,span,curve_type);
[xsort,indx] = sort(x);
% smoother_boot_sort=smootherline(indx);
% gather_xy_boot=[xsort smoother_boot_sort];
% distance = sqrt(sum(diff(gather_xy_boot,1,1).^2,2));  % Distance between subsequent points
% paracoord = [0; cumsum(distance)];               % Parametric coordinate
% coordnew = linspace(0,paracoord(end),1000).';   % 1000 evenly spaced points from 0 to s(end)
% x_new = interp1q(paracoord,gather_xy_boot(:,1),coordnew);     % Interpolate new x values
% y_new = interp1q(paracoord,gather_xy_boot(:,2),coordnew);     % Interpolate new y values
% xdata_orig(1,:)=x_new;
% ydata_orig(1,:)=y_new;
% zzz=zeros(1,size(xdata_orig,2));
% zzz=zzz+max(max(kd_array));
%plot3(xdata_orig(1,:),ydata_orig(1,:),zzz,'Color',[.8 0 0],'LineWidth',1);
plot3(xsort,smootherline(indx),zz,'Color',[.8 0 0],'LineWidth',3);

% set data color and other properties
hp=plot3(x,y,zz,'o','MarkerEdgeColor',[0 0 .8],'MarkerFaceColor',[0 0 .8],'MarkerSize', 7);
hold on;
axis tight

%plot(x_new,y_new,'k*'); % Plot interpolated bound

% This changes color of grid lines, but also lables.
% set(gca, 'XColor', [.3 .3 .3])
% set(gca, 'ZColor', [.3 .3 .3])
% set(gca, 'YColor', [.3 .3 .3])

% use this which doesnt change everything, but does grid only, but you
% can't rotate it once this is applied
%gridcolor(gca,':', [1 1 1], {[.3 .3 .3]});









end

