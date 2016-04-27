m = 1.65;
v = 2.16;
u = log((m^2)/sqrt(v+m^2));
sig = sqrt(log(v/(m^2)+1));
uo=1.65;

% log
pop=lognrnd(u,sig,1,10000);
n=20;

for i=1:50000;
    bv=randi(10000,1,n);
    samp=pop(bv);
    E=mean(samp)-uo;
    SE=std(samp)/(sqrt(n));
    T(i)=E/SE;   
end

%normal
R=normrnd(m,v,1,10000);

for i=1:50000;
    bv=randi(10000,1,n);
    samp=R(bv);
    E=mean(samp)-uo;
    SE=std(samp)/(sqrt(n));
    Tn(i)=E/SE;   
end


% desity of samp distribution
[bandwidth,f,xi,cdf]=kde(T,5000);
[bandwidth,fn,xin,cdf]=kde(Tn,5000);

a=f<.001;

an=fn<.001;

%remove zeros
f(a)=[];
xi(a)=[];
fn(an)=[];
xin(an)=[];


% plotting
figure;
plot(xin,fn,'r', 'LineWidth' ,2);
hold on
plot(xi,f,'k', 'LineWidth' ,2);
%xlim([-6 3.5]);

%% contaminated normal

% based on 10000 values
healthy=normrnd(0,1,1,90000);
schiz=normrnd(0,10.9,1,10000);

cont=[healthy schiz];
[bandwidth,fc,xic,cdf]=kde(cont,5000);

xr = [-3:.1:3];
nor = normpdf(xr,0,1);
[bandwidth,fn,xin,cdf]=kde(nor,5000);

% plotting
figure;
plot(xic,fc,'r', 'LineWidth' ,1);
hold on
plot(xr,nor,'k', 'LineWidth' ,1);
xlim([-3 3]);

%% resample twice from the contaminated normal (pretending it is the population)
% Compare the two resamples using the percentile boot test wrapped in monte carlo (there should be no difference between the two,
% becasue they are coming from the same population). See if mean, trimmed mean, or median is better in terms of type I error.

nboot=1000;
nmont=1000;
n=15;

%preallocate
delt_mean=zeros(1,nboot);
delt_trim20=zeros(1,nboot);
delt_trim10=zeros(1,nboot);
delt_med=zeros(1,nboot);

mean_pvalgen=zeros(1,nmont);
trim20_pvalgen=zeros(1,nmont);
trim10_pvalgen=zeros(1,nmont);
median_pvalgen=zeros(1,nmont);

% now we have our two samples. Conduct boot test
for m=1:nmont;
    
% take two samples from contaminated norm
bv=randi(100000,1,n);
grp1=cont(bv);
bv=randi(100000,1,n);
grp2=cont(bv);
    
    for j=1:nboot;
        
        bv=randi(n,1,n);
        bt_grp1=grp1(bv);
        bv=randi(n,1,n);
        bt_grp2=grp2(bv);
        
        % caluculate difference
        % mean
        delt_mean(j)=mean(bt_grp1)-mean(bt_grp2);
        delt_trim20(j)=trimmean(bt_grp1, 40)-trimmean(bt_grp2, 40);
        delt_trim10(j)=trimmean(bt_grp1, 10)-trimmean(bt_grp2, 10);
        delt_med(j)=median(bt_grp1)-median(bt_grp2);
    end   
 
        mean_pval=(sum(delt_mean>0)+.5*sum(delt_mean==0))/nboot;
        mean_pval=min(mean_pval,1-mean_pval);
        mean_pvalgen(m)=2*mean_pval;
    
 
        trim20_pval=(sum(delt_trim20>0)+.5*sum(delt_trim20==0))/nboot;
        trim20_pval=min(trim20_pval,1-trim20_pval);
        trim20_pvalgen(m)=2*trim20_pval;
        
        trim10_pval=(sum(delt_trim10>0)+.5*sum(delt_trim10==0))/nboot;
        trim10_pval=min(trim10_pval,1-trim10_pval);
        trim10_pvalgen(m)=2*trim10_pval;
    
 
        median_pval=(sum(delt_med>0)+.5*sum(delt_med==0))/nboot;
        median_pval=min(median_pval,1-median_pval);
        median_pvalgen(m)=2*median_pval;   
    
end

% tally up significant findings (type I error)
tal_mean=sum(mean_pvalgen<=0.05)/nmont
tal_trim20=sum(trim20_pvalgen<=0.05)/nmont
tal_trim10=sum(trim10_pvalgen<=0.05)/nmont
tal_median=sum(median_pvalgen<=0.05)/nmont

cont.mean=tal_mean;
cont.trim20=tal_trim20;
cont.trim10=tal_trim10;
cont.med=tal_median;

save('contaminated.mat','cont');

%% type I error for measures of location, based on log distributution.


% log
pop=lognrnd(u,sig,1,10000);
n=15;

%preallocate
delt_mean=zeros(1,nboot);
delt_trim20=zeros(1,nboot);
delt_trim10=zeros(1,nboot);
delt_med=zeros(1,nboot);

mean_pvalgen=zeros(1,nmont);
trim20_pvalgen=zeros(1,nmont);
trim10_pvalgen=zeros(1,nmont);
median_pvalgen=zeros(1,nmont);

% now we have our two samples. Conduct boot test
for m=1:nmont;
    
% take two samples from contaminated norm
bv=randi(10000,1,n);
grp1=pop(bv);
bv=randi(10000,1,n);
grp2=pop(bv);
    
    for j=1:nboot;
        
        bv=randi(n,1,n);
        bt_grp1=grp1(bv);
        bv=randi(n,1,n);
        bt_grp2=grp2(bv);
        
        % caluculate difference
        % mean
        delt_mean(j)=mean(bt_grp1)-mean(bt_grp2);
        delt_trim20(j)=trimmean(bt_grp1, 40)-trimmean(bt_grp2, 40);
        delt_trim10(j)=trimmean(bt_grp1, 10)-trimmean(bt_grp2, 10);
        delt_med(j)=median(bt_grp1)-median(bt_grp2);
    end   
 
        mean_pval=(sum(delt_mean>0)+.5*sum(delt_mean==0))/nboot;
        mean_pval=min(mean_pval,1-mean_pval);
        mean_pvalgen(m)=2*mean_pval;
    
 
        trim20_pval=(sum(delt_trim20>0)+.5*sum(delt_trim20==0))/nboot;
        trim20_pval=min(trim20_pval,1-trim20_pval);
        trim20_pvalgen(m)=2*trim20_pval;
        
        trim10_pval=(sum(delt_trim10>0)+.5*sum(delt_trim10==0))/nboot;
        trim10_pval=min(trim10_pval,1-trim10_pval);
        trim10_pvalgen(m)=2*trim10_pval;
    
 
        median_pval=(sum(delt_med>0)+.5*sum(delt_med==0))/nboot;
        median_pval=min(median_pval,1-median_pval);
        median_pvalgen(m)=2*median_pval;   
    
end

% tally up significant findings (type I error)
tal_mean=sum(mean_pvalgen<=0.05)/nmont
tal_trim20=sum(trim20_pvalgen<=0.05)/nmont
tal_trim10=sum(trim10_pvalgen<=0.05)/nmont
tal_median=sum(median_pvalgen<=0.05)/nmont

log.mean=tal_mean;
log.trim20=tal_trim20;
log.trim10=tal_trim10;
log.med=tal_median;

save('log.mat','log');

%% type I error rates for normal distribution

%normal
pop=normrnd(m,v,1,10000);
n=15;

%preallocate
delt_mean=zeros(1,nboot);
delt_trim20=zeros(1,nboot);
delt_trim10=zeros(1,nboot);
delt_med=zeros(1,nboot);

mean_pvalgen=zeros(1,nmont);
trim20_pvalgen=zeros(1,nmont);
trim10_pvalgen=zeros(1,nmont);
median_pvalgen=zeros(1,nmont);

% now we have our two samples. Conduct boot test
for m=1:nmont;
    
% take two samples from contaminated norm
bv=randi(10000,1,n);
grp1=pop(bv);
bv=randi(10000,1,n);
grp2=pop(bv);
    
    for j=1:nboot;
        
        bv=randi(n,1,n);
        bt_grp1=grp1(bv);
        bv=randi(n,1,n);
        bt_grp2=grp2(bv);
        
        % caluculate difference
        % mean
        delt_mean(j)=mean(bt_grp1)-mean(bt_grp2);
        delt_trim20(j)=trimmean(bt_grp1, 40)-trimmean(bt_grp2, 40);
        delt_trim10(j)=trimmean(bt_grp1, 10)-trimmean(bt_grp2, 10);
        delt_med(j)=median(bt_grp1)-median(bt_grp2);
    end   
 
        mean_pval=(sum(delt_mean>0)+.5*sum(delt_mean==0))/nboot;
        mean_pval=min(mean_pval,1-mean_pval);
        mean_pvalgen(m)=2*mean_pval;
    
 
        trim20_pval=(sum(delt_trim20>0)+.5*sum(delt_trim20==0))/nboot;
        trim20_pval=min(trim20_pval,1-trim20_pval);
        trim20_pvalgen(m)=2*trim20_pval;
        
        trim10_pval=(sum(delt_trim10>0)+.5*sum(delt_trim10==0))/nboot;
        trim10_pval=min(trim10_pval,1-trim10_pval);
        trim10_pvalgen(m)=2*trim10_pval;
    
 
        median_pval=(sum(delt_med>0)+.5*sum(delt_med==0))/nboot;
        median_pval=min(median_pval,1-median_pval);
        median_pvalgen(m)=2*median_pval;   
    
end

% tally up significant findings (type I error)
tal_mean=sum(mean_pvalgen<=0.05)/nmont
tal_trim20=sum(trim20_pvalgen<=0.05)/nmont
tal_trim10=sum(trim10_pvalgen<=0.05)/nmont
tal_median=sum(median_pvalgen<=0.05)/nmont

% get a sample of the standard error from the last bootstrap test
norm.SE=std(delt_mean,1);
norm.SE=std(delt_trim20,1);
norm.SE=std(delt_trim10,1);
norm.SE=std(delt_med,1);

norm.mean=tal_mean;
norm.trim20=tal_trim20;
norm.trim10=tal_trim10;
norm.med=tal_median;

save('norm.mat','norm');


%% comparisons testing type II error and different dists
% contaminated pops true differences
nmont=1000;
nboot=1000;

n=15;

% effect size here is roughly .5, (mean(cont2)-mean(cont1))/std(cont1,1)
% based on 10000 values
healthy=normrnd(0,1,1,90000);
schiz=normrnd(0,10.9,1,10000);
cont1=[healthy schiz];

% based on 10000 values
healthy=normrnd(1.8,1,1,90000);
schiz=normrnd(1.8,10.9,1,10000);
cont2=[healthy schiz];

% 
% [bandwidth,fc,xic,cdf]=kde(cont1,5000);
% [bandwidth,fn,xin,cdf]=kde(cont2,5000);
% 
% % plotting
% xr = [-3:.1:3];
% figure;
% plot(xic,fc,'r', 'LineWidth' ,1);
% hold on
% plot(xin,fn,'k', 'LineWidth' ,1);
% xlim([-3 3]);


%preallocate
delt_mean=zeros(1,nboot);
delt_trim20=zeros(1,nboot);
delt_trim10=zeros(1,nboot);
delt_med=zeros(1,nboot);

mean_pvalgen=zeros(1,nmont);
trim20_pvalgen=zeros(1,nmont);
trim10_pvalgen=zeros(1,nmont);
median_pvalgen=zeros(1,nmont);


% now we have our two samples. Conduct boot test
for m=1:nmont;
    
% take two samples from contaminated norm
bv=randi(100000,1,n);
grp1=cont1(bv);
bv=randi(100000,1,n);
grp2=cont2(bv);
    
    for j=1:nboot;
        
        bv=randi(n,1,n);
        bt_grp1=grp1(bv);
        bv=randi(n,1,n);
        bt_grp2=grp2(bv);
        
        % caluculate difference
        % mean
        delt_mean(j)=mean(bt_grp1)-mean(bt_grp2);
        delt_trim20(j)=trimmean(bt_grp1, 40)-trimmean(bt_grp2, 40);
        delt_trim10(j)=trimmean(bt_grp1, 20)-trimmean(bt_grp2, 20);
        delt_med(j)=median(bt_grp1)-median(bt_grp2);
        

    end   
 
        mean_pval=(sum(delt_mean>0)+.5*sum(delt_mean==0))/nboot;
        mean_pval=min(mean_pval,1-mean_pval);
        mean_pvalgen(m)=2*mean_pval;
    
 
        trim20_pval=(sum(delt_trim20>0)+.5*sum(delt_trim20==0))/nboot;
        trim20_pval=min(trim20_pval,1-trim20_pval);
        trim20_pvalgen(m)=2*trim20_pval;
        
        trim10_pval=(sum(delt_trim10>0)+.5*sum(delt_trim10==0))/nboot;
        trim10_pval=min(trim10_pval,1-trim10_pval);
        trim10_pvalgen(m)=2*trim10_pval;
    
 
        median_pval=(sum(delt_med>0)+.5*sum(delt_med==0))/nboot;
        median_pval=min(median_pval,1-median_pval);
        median_pvalgen(m)=2*median_pval;  
        
        [H,P,CI,stats]=ttest2(bt_grp1,bt_grp2,.05, 'both', 'equal');
        %delta_t(m)=stats(1);
        P(m)=P;
end

% tally up significant findings (type I error)
tal_mean=sum(mean_pvalgen<=0.05)
tal_trim20=sum(trim20_pvalgen<=0.05)
tal_trim10=sum(trim10_pvalgen<=0.05)
tal_median=sum(median_pvalgen<=0.05)

tal_t=sum(P<=0.05)

%% lognormal populations, true differences



%% normal pops true differences





%% Dance of pvalues normal curves

% created two population curves
u=50;
sig=20;
low=0;
hi=100;
step=.1;
es=10;
x1 = [low:step:hi];
pop1 = normpdf(x1,u,sig);
pop1_big=normrnd(u,sig,1,100000);

grp1=pop1;
[bandwidth,fgrp1,xigrp1,cdf]=kde(grp1,5000);


pop2_big=normrnd(u+es,sig,1,100000);

x2 = [low+es:step:hi+es];
pop1 = normpdf(x2,u+es,sig);
grp2=pop1;
[bandwidth,fgrp2,xigrp2,cdf]=kde(grp2,5000);

figure;
plot(x1,grp1, 'r');
hold on
plot(x2,grp2, 'k');
%xlim([0 100]);
axis tight

%% resampling and plotting
y=1:100;
n=30;
for i=1:100;
    bv=randi(100000,1,n);
    samp1=pop1_big(bv);
    bv=randi(100000,1,n);
    samp2=pop2_big(bv);
    
    [H(i),P(i),CI{i}]=ttest(samp1,samp2,.05, 'both');
end

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



%% Dance of pvalues contaminated normal curves


% % created two population curves
% 
% pop1=normrnd(0,1,1,90000);
% popalt=normrnd(0,10.9,1,10000);
% 
% grp1=[pop1 popalt];
% [bandwidth,fgrp1,xigrp1,cdf]=kde(grp1,5000);
% 
% pop2=normrnd(1,1,1,90000);
% popalt=normrnd(1,10.9,1,10000);
% 
% grp2=[pop2 popalt];
% [bandwidth,fgrp2,xigrp2,cdf]=kde(grp2,5000);
% 
% figure;
% plot(xigrp1,fgrp1,'r', 'LineWidth' ,1);
% hold on
% plot(xigrp2,fgrp2,'k', 'LineWidth' ,1);
% xlim([-3 4]);
% 
% %%
% n=25;
% for i=1:100;
%     bv=randi(100000,1,n);
%     samp1=grp1(bv);
%     bv=randi(100000,1,n);
%     samp2=grp2(bv);
%     
%     [H(i),P(i),CI{i}]=ttest(samp1,samp2,.05, 'both');
% end
% 
% figure;
% y=1:100;
% for i=1:100;
% 
% plot(CI{i}, [y(i) y(i)])
% hold on
% end
% 
% clear a
% dif=u-(u+es);
% j=1;
% for i=1:100;
%     if CI{i}(1)<=dif && CI{i}(2)>=dif;
%         a(j)=i;
%         j=j+1;
%     end
% end


%% sounds
psound=cell(length(P),2);

P1=find(P<.001);
P2=find(P>=.001 & P<.01);
P3=find(P>=.01 & P<.05);

P4=find(P>=.05 & P<.1);
P5=find(P>=.1 & P<.5);
P6=find(P>=.5);


psound(P1,1)={'G4'};
psound(P2,1)={'A4'};
psound(P3,1)={'B4'};

psound(P4,1)={'G3'};
psound(P5,1)={'A3'};
psound(P6,1)={'B3'};

psound(P1,2)={.5};

psound(P2,2)={.5};
psound(P3,2)={.5};

psound(P4,2)={.5};
psound(P5,2)={.5};
psound(P6,2)={.5};
[ys]=nplay(psound);
wavwrite(ys,8192,'mypvalues.wav');
     
% f=261.63
% Fs=8192; % Sampling frequency.
% t=0:1/Fs:.25;
% y=sin(2*pi*f*t);
% 
% for i=1:100
% % --- Sound off ---
% %sound(repmat(y,1,100),Fs)
% beep2(261.63,.25)
% end













