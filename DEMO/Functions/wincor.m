function [results] = wincor(nboot,tr,Xdata,Ydata);
tic

 
% Winsorized and bootstrapped correlation (rw).
% 
%     Input arguments:
%     nboot = number of bootstrap resamples (e.g., 1000)
%     tr = percentage of winsorizing (e.g., .2 for 20%)
% 
% You will be prompted to load two .mat files one at a time (X and Y data). 
% These are N X M matrices where rows = subjects, and columns = variables. 
% Each Y column is compared to each X column. The results structure is
% organized by Y comparisons (eg., Y1 vs X1 X2...Xn, Y2 vs X1 X2 Xn etc)


if isempty(Xdata) && isempty(Ydata)

    x_fname=uigetfile('*.mat','select X variable (EEG condition data)', 'MultiSelect','off');
    y_fname=uigetfile('*.mat','select Y variable (potential EEG correlate(s))', 'MultiSelect','off');

%% load data for x and y variables (they can be nxm matrices). Convert to cell arrays if input is a structure.

    x=load(x_fname);
    if isstruct(x);
        xy_data(1,1)=struct2cell(x);
    else
        xy_data{1,1}=x;
    end


    y=load(y_fname);
    if isstruct(y);
        xy_data(1,2)=struct2cell(y);
    else
        xy_data{1,2}=y;
    end
    
    
else
    
    % else if wincor is called from another function that gave X and Y data as inputs
    xy_data{1,1}=Xdata;
    xy_data{1,2}=Ydata;
    
end






%% preallocate sizes
[row_y col_y]=size(xy_data{1,2});
[row_x col_x]=size(xy_data{1,1});
rw_surrogates=zeros(nboot,1);
up_ind=floor(.975*nboot+.5);
low_ind=floor(.025*nboot+.5);
bootvect=randi(row_y,nboot,row_y); % predefine inds for resampling, so every comparison uses same resamples

%% build results structure

for i=1:col_y;
    field_name{i,1}=['Y', num2str(i)];
end

for i=1:col_y;
    results.(field_name{i})=struct('rw',{[]},'CI',{[]},'p',{[]},'distribution',{[]});
    results.(field_name{i}).CI=zeros(2,col_x);
    results.(field_name{i}).distribution=zeros(nboot,1);
end

%% arrange data into an nx2 array

h1 = waitbar(0,'1','Name','comparing X against each Y variable','Position',[1100 486 550 40]);
childh1 = get(h1, 'Children');
set(childh1, 'Position',[5 10 538 15]);

h2 = waitbar(0,'1','Name','winsorizing & bootstrapping','Position',[1100 545 550 40]);
childh2 = get(h2, 'Children');
set(childh2, 'Position',[5 10 538 15]);


for ycomp=1:col_y;
    
    for xcomp=1:col_x;
        
        % reset data before each new arrangment
        data=zeros(row_y,2); % allocating for original data
        
        % arranging into a two column matrix
        data(:,1)=xy_data{1,1}(:,xcomp);
        data(:,2)=xy_data{1,2}(:,ycomp);
        
        %% calculate rw -> winsorize data and do pearson's r
        %% THIS IS JUST FOR PLOTTING PURPOSES
        [windata] = winsorize(data,tr);
        [rw]=corr(windata);
        
        
        %% bootstrapping
        for bootcurrent=1:nboot;
            
            % resampling using random indices
            bootdata=data(bootvect(bootcurrent,:),:);
            
            %% calculate rw (based on bootstrapped data)
            [windata_boot] = winsorize(bootdata,tr);
            [rw_boot]=corr(windata_boot);
            rw_surrogates(bootcurrent,1)=rw_boot(2,1);
            
            clear windata_boot rw_boot
            
            %waitbar(bootcurrent/nboot,h2,sprintf('%12s',[num2str(bootcurrent),'/',num2str(nboot)]))
        end
        
        
        %% Calculate CIs
        rw_bootsort=sort(rw_surrogates);
        conf_up=rw_bootsort(up_ind);
        conf_low=rw_bootsort(low_ind);
        
        %% Calculate generalized p values
        pval_init=sum(rw_surrogates<0)/nboot;
        pval_gen=2*min(pval_init,1-pval_init);
        
        %% Populate results structure
        
        % store original rw
        results.(field_name{ycomp}).rw(1,xcomp)=rw(2,1);
        
        % CIs
        results.(field_name{ycomp}).CI(1,xcomp)=conf_low;
        results.(field_name{ycomp}).CI(2,xcomp)=conf_up;
        
        % p values
        results.(field_name{ycomp}).p(1,xcomp)=pval_gen;
        
        % populate distributions
        results.(field_name{ycomp}).distribution(:,xcomp)=rw_surrogates(:,1);
       
    waitbar(xcomp/col_x,h2,sprintf('%12s',[num2str(xcomp),'/',num2str(col_x)]))    
    end
    
    waitbar(ycomp/col_y,h1,sprintf('%12s',[num2str(ycomp),'/',num2str(col_y)]))
end

close(h1,h2);
toc
end

