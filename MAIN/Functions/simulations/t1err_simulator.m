function [sim]=t1err_simulator(sampsize,varargin)
% measure the type I error of one-sample t-test for a number of underlying population distibutions
%
% inputs:
% sampsize -> [numeric] sample size
%
% options:
%
% 'popdist' -> 'normal', 'log', 'contaminated', 'skewkurt (depreciated)', 'gandh', 'onewild', 'slash', 'exp', or 'chi2contam'
% 'popstd' -> [numeric] the standard deviation of the population
% 'g' -> [numeric] range 0-1 amount of skew
% 'h' -> [numeric] range 0-1, increase tail thickness
%
% depreciated
% 'contam' -> [numeric] the standard deviation of the contaminated portion of the distribution. Only used with contaminated option
% 'skew' -> [numeric] amount of skew, negative values mean negative skew
% 'kurt' -> [numeric] amount of kurtosis (outliers). kurt=3 is kurtosis of normal curve

%% defaults and general params
nboot=10000;
options.popdist='normal';
options.popstd=10;
%options.contam=100; % aim for 10x the SD of the normal population
%options.skew=0;
%options.kurt=3;
options.g=0;
options.h=0;
options.df=4; % for chi
m = 1;
v=options.popstd;

popsize=100000;
n=sampsize;

optionNames = fieldnames(options);

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
    inpName = pair{1};
    
    if any(strcmp(inpName,optionNames))
        
        % overwrite default options
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name',inpName)
    end
end

diststr=options.popdist;

%% cases
switch options.popdist
    
    case 'normal'
        disp('normal dist chosen');
        pop=normrnd(m,v,1,popsize);
        
        
    case 'log'
        
        disp('log dist chosen');
        %m=1.649; % wont work if 0;
        
        % set this so that it matches with wlicox (1998) where he does not specify the VAR of the lognormal dist
        % he mentions that the actual T dist has mean of -5, and this VAR value leads to that mean of T
        m=1.649;
        v=4.5;
        u = log((m^2)/sqrt(v+m^2));
        sig = sqrt(log(v/(m^2)+1));
        pop=lognrnd(u,sig,1,popsize);
        
        
    case 'contaminated'
        disp('contaminated dist chosen');
        
        % set to match wilcox 2003;
        m=0;
        v=1;
        options.contam=10;
        
        perc=round(.1*popsize);
        healthy=normrnd(m,v,1,popsize-perc);
        schiz=normrnd(m,options.contam,1,perc);
        pop=[healthy schiz];
        
        % add case to control skew and kurtosis
        %     case 'skewkurt'
        %         disp('skewkurt option chosen');
        %        % kurt (norm dist has kurt 3) >3 means more outliers than normal
        %        % skew is zero (centered) for a normal or symetric curve)
        %        % negative skew value = negatively skewed range roughly [-1 1] for u=o, sd=1?
        %          moments = {0,options.popstd,options.skew,options.kurt};
        %         [pop,type] = pearsrnd(moments{:},100000,1);
        
    case 'gandh'
        m=0;
        disp('g and h class distribution');
        pop=gh_simulator(options.g,options.h);
        
    case 'onewild'
        disp('one wild night');
        pop=normrnd(m,v,1,popsize);
        
        % Contributions to Survey Sampling and Applied Statistics (David HA, 1978)
        % Robustness of Location Estimators in the Presence of an Outlier, Pages 235-250
        
    case 'slash'
        disp('welcome to the jungle');
        pop=normrnd(m,v,1,popsize);
        bv=rand(1,popsize);
        pop=pop./bv;
        
        
    case 'exp'
        m=1; % 0 wont work
        disp('exponential dist chosen');
        pop=exprnd(m,1,popsize);
        
    case 'chi2contam'
        disp('wilcox and keselman, 2003')
        %y = chi2pdf(x,4);
        pop=chi2rnd(options.df,1,popsize);
        uni=rand(1,popsize);
        ind=uni<.1;
        pop(ind)=pop(ind)*10;
        %m=mean(pop);
        m=7.6; % the mean of the dist according to Wilcox
end

%% resampling
ts=zeros(1,nboot);
if strcmp(options.popdist, 'onewild')
    
    for j=1:nboot;
        
        % make it wild
        bv=randi(popsize,1,n);
        bt_grp1=pop(bv);
        bs=randi(n,1,1);
        bt_grp1(bs)=bt_grp1(bs)*10;
        
        % T based on chosen distibution
        Ep=mean(bt_grp1)-m;
        SE=std(bt_grp1)/(sqrt(n));
        ts(j)=Ep/SE;
    end
    
    
else
    
    for j=1:nboot;
        
        % based on chosen distibution
        bv=randi(popsize,1,n);
        bt_grp1=pop(bv);
        Ep=mean(bt_grp1)-m;
        SE=std(bt_grp1)/(sqrt(n));
        ts(j)=Ep/SE;
        
        % based on normal distribution
        %         bv=randi(popsize,1,n);
        %         samp=R(bv);
        %         En=mean(samp)-m;
        %         SE=std(samp)/(sqrt(n));
        %         tn(j)=En/SE;
    end
    
end

%% type I error
% find critical values for each tail of T dist with n-1 df
% one-tailed test quantiles!
tcrit_l=tinv(.05,n-1);
tcrit_r=tinv(.95,n-1);

% looking at the actual dist of T (bootstrapped),
% find probability of values exceeding the conventional critical values
pr=invprctile(ts,tcrit_r)/100;
pr=1-pr;
pl=invprctile(ts,tcrit_l)/100;
totp=pr+pl;

% CI of the bootstrapped estimated ts
CIlow=round(.05*nboot/2)+1;
CIup=nboot-CIlow-1;
temp=sort(ts(1,:));
conflow=temp(CIlow);
confup=temp(CIup);

%% keep track of parameters
sim.(diststr).CIassumed_95=[tinv(.025,n-1) tinv(.975,n-1)]; 
sim.(diststr).CI_95=[conflow confup]; 
sim.(diststr).t1err_05.rtail=pr;
sim.(diststr).t1err_05.rtail=pl;
sim.(diststr).t1err.total=totp;
sim.params.nboot=nboot;
sim.params.sampsize=sampsize;
sim.params.popdist=options.popdist;
sim.params.popstd=options.popstd;
sim.params.skew=options.g;
sim.params.kurt=options.h;
%sim.params.contam=options.contam;
save([diststr, 't1err_simulator.mat'],'sim');

%% figures

% errors
bar(1,pr,'r');
hold on
bar(2,pl,'k');
title('actual type I error','Interpreter', 'none');
legend('right tail','left tail');

% desity of actual distribution
[bandwidth,f,x,cdf]=kde(ts,500);

%plotting
figure;
subplot(2,1,1)
plot(x,f,'r', 'LineWidth' ,2);
hold on

xax = -5:0.1:5;
Tdist = tpdf(xax,n-1);
plot(xax,Tdist,'-k','linewidth',2);
title(['population ', diststr],'Interpreter', 'none');
legend('actual T','assumed T');
xlim([-5 5]);

subplot(2,1,2)
plot([conflow confup],[.25 .25],'.-r')
hold on
plot([tinv(.025,n-1) tinv(.975,n-1)],[.75 .75],'.-k');
xlim([-5 5]);
ylim([0 1]);

%
% figure;
% bar(1,sim.(diststr).FDR.mean,'r');
% hold on
% bar(2,sim.(diststr).FDR.trim20,'k');
%
% title(['type I error ',diststr],'Interpreter', 'none');
% legend('mean','trim20');

end



















