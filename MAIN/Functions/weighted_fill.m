function weighted_fill(xval,yup, ylow, CI_color, colorlimit) 

%weighted_fill_smoother(xup, xlow, yup, ylow, CI_color, colorlimit, alpha_trans) 
% plot smoother data with weighted fill

ydiff=abs(yup-ylow); 
max_color = max(ydiff);
min_color = min(ydiff);
color_range = (max_color - min_color);
ydiff = (ydiff - min_color)/color_range; 
ydiff=ones(size(ydiff))-ydiff; 
C=ydiff;

%% Y1 and Y2 do not have to be sorted with lowess curves
%Y1=sort(yup);
%Y2=sort(ylow);

Y1=yup;
Y2=ylow;

%% Defining X, 
% this lines just average the interpolated upper and lower bounds for the xaxis values
%xx=horzcat(xup, xlow);
%avg_xx=mean(xx,2);
X=xval; 
%X=xvals_common;
%X=min(x):0.01:max(x); % WHAT AM I GOING TO DO WITH THIS? linespace?
%minx=min(xvals_common);
%maxx=max(xvals_common);
%X=minx:0.01:maxx;
%X=linspace(minx, maxx,length(xvals_common));
%figure;
%colorlimit = [0 0 0];
%CI_color = [1 1 1]; %white for low values
%whitebg([0 0 0]); % change backround of plot area
%FigHandle = figure('Color',[1 1 1]); % change border backround color here
%set(FigHandle, 'Position', [100, 100, 1049, 895]);


%% do visual weighting

for i = 1:length(X)-1
    h=fill([X(i), X(i+1), X(i+1), X(i)],[Y1(i), Y1(i+1), Y2(i+1), Y2(i)], ...
        C(i)*(CI_color-colorlimit)+colorlimit,...
        'EdgeColor', (C(i)*(CI_color-colorlimit)+colorlimit));
    set(h,'FaceAlpha',1) % 1 is a transparency thing, input argument?
    set(h,'EdgeAlpha',1)
    
    if i == 1
        hold on
    end
end
axis tight
grid on
box on
%axis([10,51.01,min(y),max(y)])


%FigHandle = figure('Color',[1 1 1]);
%set(FigHandle, 'Position', [100, 100, 1049, 895]);
%{
for i = 1:length(X)-1
    fill([X(i), X(i+1), X(i+1), X(i)],[Y1(i), Y1(i+1), Y2(i+1), Y2(i)], ...
        nanmean(C(i:i+1))*(CI_color-colorlimit)+colorlimit,...
        'EdgeColor', nanmean(C(i:i+1))*(CI_color-colorlimit)+colorlimit)
    
    if i == 1
        hold on
    end
end
grid on
axis tight
%}

%{
switch plot_data
    
    case 'scatter'
        hold on
        plot(x,y,'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize', 5);
        grid on
        axis tight
        
    case 'smoother'
        hold on
        [xsort,indx] = sort(x);
        plot(xsort,smoother(indx),'Color','r','LineWidth',2);
        
    case 'both'
        hold on
        plot(x,y,'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize', 5);
        hold on
        [xsort,indx] = sort(x);
        plot(xsort,smoother(indx),'Color','r','LineWidth',2);
        grid on
        axis tight       
end
%}

end




