function CI_PB_slope(x,y,CI_color,colorlimit,Xlabel,Ylabel)

%{

x and y = N X 1 column vectors of your data
CI_color = rgb triple e.g., [.4 .4 .4]
colorlimit = the color that wide CIs will fade to
             e.g., [1 1 1] -> fades to white
%}



[row_dat col_dat]=size(x);
n=row_dat;
meanx=mean(x);
%meany=mean(y);
stdx=std(x);
%stdy=std(y);

% get some stats that you need 
stats_reg=regstats(y,x,'linear',{'beta', 'mse'});

%standard error of residual
SE_resid=sqrt(stats_reg.mse);

%critical tvalue (alpha = .05)
tcrit=tinv(1-0.05/2,n-2);

xval = min(x):0.01:max(x);
yhat = stats_reg.beta(1)+stats_reg.beta(2)*xval;

%preallocating
[row col]=size(yhat);
yup=zeros(1,col);
ylow=zeros(1,col);
yup_p=zeros(1,col);
ylow_p=zeros(1,col);

for i=1:col;
        
    %Confidence bands
    
    div=((xval(i)-meanx)^2)/((n-1)*(stdx^2));
    sqpl=sqrt((1/n)+div);
    rng=abs(tcrit*SE_resid*sqpl);
    yup(i)=yhat(i)+rng;
    ylow(i)=yhat(i)-rng;
    
    %Prediction bands
  
    sqpl_p=sqrt((1+(1/n)+div));
    rng_p=abs(tcrit*SE_resid*sqpl_p);
    yup_p(i)=yhat(i)+rng_p;
    ylow_p(i)=yhat(i)-rng_p;
    
end

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1049, 895]);

%set(gca,'FontSize',20)

%hold on
weighted_fill(xval,yup, ylow, CI_color, colorlimit)

% plot scatter data
hold on
hp=plot(x,y,'ko','MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize', 7);     
hold on;
grid on
axis tight

set(gca,'FontSize',20)

% not shaded CI
%hold on
%h(2)=plot(xval,ylow,'k-.', 'LineWidth', 2);
%h(3)=plot(xval,yup,'k-.', 'LineWidth', 2);

% prediction bands
h(4)=plot(xval,ylow_p,'Color', CI_color);
h(5)=plot(xval,yup_p,'Color', CI_color);

% regression line
h(6)=plot(xval,yhat,'w','linewidth',2);
axis tight

% add other strings  in curly braces if you have more items in legend
%leg=legend(h([2 4]),{'Confidence band','Prediction band'});
%leg=legend(h,{'Confidence band'});
%set(leg,'Location','northwest')
%set(leg,'FontSize',18);
h(7)=ylabel(Ylabel, 'fontsize', 18);
h(8)=xlabel(Xlabel, 'fontsize', 18);

set(h([7 8]),'Interpreter', 'none');


end









