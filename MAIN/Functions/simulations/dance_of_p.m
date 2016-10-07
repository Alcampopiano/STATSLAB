function dance_of_p

%% Dance of pvalues normal curves

% created two population curves
u=0; %50
sig=1; %20
low=0;
hi=100;
step=.1;
es=0;  %10
g=0;
h=0;
n=20;
x1 = [low:step:hi];


%X=gh_simulator_shift(g,h,mu,sig,shift);
%pop1=gh_simulator_shift(0,0,0,sig,0);
pop1_big=gh_simulator_shift(g,h,u,sig,0);
[bandwidth,fgrp1,xigrp1,cdf]=kde(pop1_big,5000);
%pop1 = normpdf(x1,u,sig);
%pop1_big=normrnd(u,sig,1,100000);
%grp1=pop1_big;


pop2_big=gh_simulator_shift(g,h,u,sig,es);
[bandwidth,fgrp2,xigrp2,cdf]=kde(pop2_big,5000);
%pop2_big=normrnd(u+es,sig,1,100000);
x2 = [low+es:step:hi+es];
%pop1 = normpdf(x2,u+es,sig);
%pop2=gh_simulator_shift(0,0,50,20,es);
%grp2=pop2;

figure;
%plot(x1,grp1, 'r');
%[bandwidth,f,x,cdf]=kde(X,5000);
plot(xigrp1,fgrp1,'r', 'LineWidth' ,2);

%xlim([-5 5]);

hold on
%plot(x2,grp2, 'k');
plot(xigrp2,fgrp2,'k', 'LineWidth' ,2);
%xlim([0 100]);
axis tight

%% resampling and plotting
y=1:100;

%bootinds
lowind=round(.05*1000/2);
upind=1000-lowind;

for i=1:100;
    bv=randi(100000,1,n);
    samp1=pop1_big(bv);
    bv=randi(100000,1,n);
    samp2=pop2_big(bv);
    
    
    % compute Student's t
    [H(i),P(i),CI{i}]=ttest(samp1,samp2,.05, 'both');
    
    
    for j=1:1000;
    % compute percentile bootstrap
    bv=randi(n,1,n);
    samp1_bt=samp1(bv);
    bv=randi(n,1,n);
    samp2_bt=samp2(bv);
    
    % take random estimator based on boot samples
    tr1=trimmean(samp1_bt,40);
    tr2=trimmean(samp2_bt,40);
    dif(j)=tr1-tr2;
    end
    
    dif=sort(dif);
    CIbt{i}(1)=dif(lowind+1);
    CIbt{i}(2)=dif(upind);
    A=sum(dif>0);
    phat=A/1000;
    phat=min(phat, 1-phat);
    pgen(i)=2*phat;
    
end


% dealing with plots for traditional test
j=1;
for i=1:100;
    if CI{i}(1)<=-es && CI{i}(2)>=-es;
        a(j)=i;
        j=j+1;
    end
end

% plot pvalues
figure; barh(P);
pcrit=zeros(1,100)+.05;
hold on
plot(pcrit,y ,'r', 'LineWidth' ,1);
axis tight

% find CIs that don't include es
match=setdiff(y,a);

figure;
for i=1:100;
    
    if any(i==match)
        plot(CI{i}, [y(i) y(i)], 'Color', 'r', 'Marker' ,'.');
    else
        plot(CI{i}, [y(i) y(i)], 'Marker' ,'.');
    end
    hold on
end

% plot es
esx=u-(u+es);
esx=zeros(1,100)+esx;
plot(esx,y, 'k');



% dealing with plots for percentile boot test
j=1;
for i=1:100;
    if CIbt{i}(1)<=-es && CIbt{i}(2)>=-es;
        a(j)=i;
        j=j+1;
    end
end

% plot pvalues
figure; barh(pgen, 'g');
pcrit=zeros(1,100)+.05;
hold on
plot(pcrit,y ,'r', 'LineWidth' ,1);
axis tight

% find CIs that don't include es
match=setdiff(y,a);

figure;
for i=1:100;
    
    if any(i==match)
        plot(CIbt{i}, [y(i) y(i)], 'Color', 'r', 'Marker' ,'.');
    else
        plot(CIbt{i}, [y(i) y(i)], 'Color', 'g', 'Marker' ,'.');
    end
    hold on
end

% plot es
esx=u-(u+es);
esx=zeros(1,100)+esx;
plot(esx,y, 'k');


percboot=sum(pgen<=.05)
tradtest=sum(P<=.05)

end


function X=gh_simulator_shift(g,h,mu,sig,shift)

Z=normrnd(mu,sig,1,100000);


% desity of pops distribution
%[bandwidth,f,x,cdf]=kde(Z,5000);

%figure;
%plot(x,f,'k', 'LineWidth' ,2);
if g~=0;
    
    t1=exp(g*Z)-1;
    t2=g;
    t3=exp(h*Z.^2/2);
    X=(t1/t2).*t3;
    
elseif g==0;
    
    X=Z.*exp(h*Z.^2/2);
end

X=X+shift;

%hold on
% desity of pops distribution
%[bandwidth,f,x,cdf]=kde(X,5000);
%plot(x,f,'r', 'LineWidth' ,2);

%xlim([-5 5]);
end



