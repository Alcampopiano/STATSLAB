function [sim]=SE_simulator(sampsize,varargin)

% measure the standard error using bootstrap estimate for mean, median, and
% 20% trimmed mean as a function of various population distributions

% inputs:
% sampsize -> [numeric] sample size

% options:

% 'popdist' -> 'normal', 'log', 'contaminated', 'skewkurt (depreciated)', 'gandh', 'onewild', 'slash', 'exp', or 'chi2contam'
% 'popstd' -> [numeric] the standard deviation of the population
% 'g' -> [numeric] range 0-1 amount of skew
% 'h' -> [numeric] range 0-1, increase tail thickness
% 
% depreciated
% 'contam' -> [numeric] the standard deviation of the contaminated portion of the distribution. Only used with contaminated option
% 'skew' -> [numeric] amount of skew, negative values mean negative skew
% 'kurt' -> [numeric] amount of kurtosis (outliers). kurt=3 is kurtosis of normal curve

% general params
nboot=50000;
options.popdist='normal';
options.popstd=10;
options.contam=10; % aim for 10x the SD of the normal population
%options.skew=0;
%options.kurt=3;
options.g=0;
options.h=0;
options.df=4; % for chi
m = 0;
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

switch options.popdist
    
    case 'normal'
        disp('normal dist chosen');
        pop=normrnd(m,v,1,popsize);
        
        
    case 'log'
        disp('log dist chosen');
        m=1; % wont work if 0;
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
        %schiz=normrnd(m,v/contam,1,perc);% if pop=1, pop2 std should = 10.9
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
        disp('g-and-h class distribution');
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
        % for i=1:nboot;
        pop=chi2rnd(options.df,1,popsize);
        uni=rand(1,popsize);
        ind=uni<.1;
        pop(ind)=pop(ind)*10;
        
end

bt_mean=zeros(1,nboot);
bt_trim20=zeros(1,nboot);
bt_med=zeros(1,nboot);

if strcmp(options.popdist, 'onewild')
    
    for j=1:nboot;
        bv=randi(popsize,1,n);
        bt_grp1=pop(bv);
        bs=randi(n,1,1);
        bt_grp1(bs)=bt_grp1(bs)*10;
        bt_mean(j)=mean(bt_grp1);
        bt_trim20(j)=trimmean(bt_grp1, 40);
        bt_med(j)=median(bt_grp1);
    end
    
    
else
    
    for j=1:nboot;
        bv=randi(popsize,1,n);
        bt_grp1=pop(bv);
        bt_mean(j)=mean(bt_grp1);
        bt_trim20(j)=trimmean(bt_grp1, 40);
        bt_med(j)=median(bt_grp1);
    end
    
end

%calculate SE
SE_m=std(bt_mean,1);
SE_t20=std(bt_trim20,1);
SE_me=std(bt_med,1);


sim.(diststr).SE.mean=SE_m;
sim.(diststr).SE.trim20=SE_t20;
sim.(diststr).SE.med=SE_me;

% keep track of parameters
sim.params.nboot=nboot;
sim.params.sampsize=sampsize;
sim.params.popdist=options.popdist;
sim.params.popstd=options.popstd;
%sim.params.skew=options.skew;
%sim.params.kurt=options.kurt;
sim.params.contam=options.contam;

save([diststr, '_simulations.mat'],'sim');

%% figures

%SE
figure;
bar(1,sim.(diststr).SE.mean,'r');
hold on
bar(2,sim.(diststr).SE.trim20,'k');
bar(3,sim.(diststr).SE.med,'m');

title(['SE ',diststr],'Interpreter', 'none');
legend('mean','trim20','median');

% desity of pops distribution
[bandwidth,f,x,cdf]=kde(pop,5000);

%plotting
figure;
plot(x,f,'r', 'LineWidth' ,2);
title(['population ', diststr],'Interpreter', 'none');

end



















