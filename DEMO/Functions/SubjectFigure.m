function [STATS]=SubjectFigure(STATS,infodisplay,plottype,varargin)
% build figure for group results
% varargin can be 'all' which plots all contrasts for all factors, or
% specify what you want with name/val pairs:
% ...'FactorA', [1 2 3], 'FactorB', [1], 'FactorAB', [1])
% Defaults for contrasts to use are empty, so ommiting one of the factors will not produce any plot for that one.
% infodisplay is a flag 1 or 0 that will spit out the contrast matrices and condition lables for you to look at in the command window
% infodisplay default is set to 0.

%plottype ='wave', or, 'diff', or 'CI_MOE'
% see documentation for other optional inputs


if isempty(STATS)
    [fnamestats]=uigetfile('*.mat','Select the STATS structure for this analysis');
    load(fnamestats);
else
    load(STATS);
end

% special case where using 'all' plots all possible contrasts
if any(strcmp(varargin,'all'));
    
    % this finds out if the design was factorial or not and sets default
    % options accordingly
    if strcmp(STATS.design,'ww');
        
        options = struct('FactorA', 1:size(STATS.subject_results.subject_1.factor_A.contrasts,2), ...
            'FactorB', 1:size(STATS.subject_results.subject_1.factor_B.contrasts,2), ...
            'FactorAB', 1:size(STATS.subject_results.subject_1.factor_AxB.contrasts,2));
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.subject_results.subject_1.factor_A.contrasts)
            disp('FactorB'); disp(STATS.subject_results.subject_1.factor_B.contrasts)
            disp('FactorAB'); disp(STATS.subject_results.subject_1.factor_AxB.contrasts)
        end
        
    elseif strcmp(STATS.design,'w');
        
        options = struct('FactorA', 1:size(STATS.subject_results.subject_1.factor_A.contrasts,2));
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.subject_results.subject_1.factor_A.contrasts)
        end
        
    elseif strcmp(STATS.design,'bw'); % considered not factorial in single-subject cases
        
        options = struct('FactorA', 1:size(STATS.subject_results.Factor_A1.subject_1.factor_A.contrasts,2));
        
        if infodisplay
            disp('j level labels'); disp(STATS.jlabels)
            disp('k level labels'); disp(STATS.klabels)
            disp('contrasts'); disp(STATS.subject_results.Factor_A1.subject_1.factor_A.contrasts)
        end
        
    end
    
    % set default plot options
    options.diffcol=[0 0 0];
    options.CIcol=[.5 .5 .5];
    options.yaxis='auto';
    options.caxis='auto';
    options.zaxis=[0 inf];
    
    % overwrite options with plot options
    % read the acceptable names
    optionNames = fieldnames(options);
    
    % get rid of 'all' option in varargin
    ind=find(strcmp('all', varargin));
    varargin(ind)=[];
    
    
    for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
        inpName = pair{1};
        
        if any(strcmp(inpName,optionNames))
            
            % overwrite default options
            options.(inpName) = pair{2};
        else
            error('%s is not a recognized parameter name',inpName)
        end
    end
    
    
    % else check for name/val agreement and set defaults to empty
else
    nargs = length(varargin);
    if round(nargs/2)~=nargs/2
        error('need propertyName/propertyValue pairs for optional inputs')
    end
    
    % this finds out if the design was factorial or not and sets default
    % options accordingly
    if strcmp(STATS.design,'ww');
        
        options = struct('FactorA', [],'FactorB', [], 'FactorAB', []);
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.subject_results.subject_1.factor_A.contrasts)
            disp('FactorB'); disp(STATS.subject_results.subject_1.factor_B.contrasts)
            disp('FactorAB'); disp(STATS.subject_results.subject_1.factor_AxB.contrasts)
        end
        
    elseif strcmp(STATS.design,'w');
        
        options = struct('FactorA', []);
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.subject_results.subject_1.factor_A.contrasts)
        end
        
    elseif strcmp(STATS.design,'bw'); % considered not factorial in single-subject cases
        
        options = struct('FactorA', []);
        
        if infodisplay
            disp('j level labels'); disp(STATS.jlabels)
            disp('k level labels'); disp(STATS.klabels)
            disp('contrasts'); disp(STATS.subject_results.Factor_A1.subject_1.factor_A.contrasts)
        end
        
    end
    
    % set default plot options
    options.diffcol=[0 0 0];
    options.CIcol=[.5 .5 .5];
    options.yaxis='auto';
    options.caxis='auto';
    options.zaxis=[0 inf];
    
    % read the acceptable names
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
    
end

% update STATS structure
STATS.subplotoptions=options;

switch STATS.design
    case 'w'
        
        if ~isempty(options.FactorA);
            
            % get subject fields
            subnames = fieldnames(STATS.subject_results);
            
            numfigs=size(options.FactorA,2);
            
            for i=1:numfigs
                k=1;
                m=1;
                
                for j=1:length(STATS.condnames);
                    
                    % get the condition waveforms
                    if STATS.subject_results.subject_1.factor_A.contrasts(j,options.FactorA(i))==1
                        %plot1st(k,:)=STATS.condwaves_trim(j,:);
                        lab1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.subject_results.subject_1.factor_A.contrasts(j,options.FactorA(i))==-1
                        %plot2nd(m,:)=STATS.condwaves_trim(j,:);
                        lab2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                end
                
                % concatenate, if needed, the legend lables
                lab1st=strjoin(lab1stlist,'+');
                lab2nd=strjoin(lab2ndlist,'+');
                
                switch plottype
                    case 'wave'
                        for q=1:length(subnames);
                            
                            % extract CI
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,:);
                            
                            % extract difference wave
                            plotdiff=STATS.subject_results.(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                            
                            % begin plotting
                            figure(i);
                            subplot(length(subnames),1,q);
                            
                            % plot zeroline
                            plot(STATS.xtimes,zeros(1,length(STATS.xtimes)),'r','LineWidth',1);
                            hold on
                            
                            % plot CI
                            jbfill(STATS.xtimes,CI_up,CI_low,options.CIcol,options.CIcol,1,1);
                            hold on
                            
                            % plot diff wave
                            plot(STATS.xtimes,plotdiff,'Color',options.diffcol);
                            ylim(gca,options.yaxis);
                            grid on
                            
                        end
                        
                    case 'CI_MOE'
                        
                        
                        for q=1:length(subnames);
                            
                            % extract MOE
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,:);
                            CI_diffup=CI_up-CI_low;
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                            
                            % begin plotting
                            figure(i);
                            subplot(length(subnames),1,q);
                            surf([STATS.xtimes;STATS.xtimes], ones(2,length(STATS.xtimes)), [zeros(1,length(CI_diffup));CI_diffup],[test_stat;test_stat], 'LineStyle', 'none', 'FaceColor', 'interp');
                            zlim(gca,options.zaxis);
                            hold on
                            
                            % find where CI includes zero
                            sigvect=ones(1,STATS.numpnts);
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,:);
                            
                            for m=1:STATS.numpnts;
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs));
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(1));
                            Z=zeros(1,length(timevect));
                            scatter3(timevect,zeroplace,Z,10,'s','filled','MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
                            
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                            
                            view(0,0)
                            hold on
                            grid on
                        end
                        
                    case 'diff'
                        
                        X=linspace(min(STATS.xtimes),max(STATS.xtimes),length(STATS.xtimes));
                        for q=1:length(subnames);
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                            
                            % begin plotting
                            figure(i);
                            subplot(length(subnames),1,q);
                            imagesc(X,[],test_stat(1,:))
                            hold on
                            set(gca,'YTickLabel',[]);
                            
                            % find where CI includes zero
                            sigvect=ones(1,STATS.numpnts);
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,:);
                            
                            for m=1:STATS.numpnts;
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs))+1.5;
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(1));
                            scatter(timevect,zeroplace,10, 'k', 'MarkerFaceColor', 'k')
                            
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                            
                        end
                end
                
                if ~strcmp(plottype,'wave') && ~strcmp(options.caxis,'auto')
                    ha=axes('visible', 'off');
                    set(ha, 'Units', 'Normalized', 'Position', [.9, 0.05, .015, .9]);
                    colorbar('FontSize',15);
                    caxis(options.caxis)
                end
                
                % add a title to the parent figure, and control its position
                h = axes('visible','off');
                title([lab1st,'-',lab2nd],'parent',h,'visible','on','Position',[.5 1.05 0]);
                clear lab1st lab2nd lab1stlist lab2ndlist
            end
        end
        
    case 'ww'
        % factor A
        numfigsA=size(options.FactorA,2);
        numfigsB=size(options.FactorB,2);
        
        
        if ~isempty(options.FactorA);
            
            % get subject fields
            subnames = fieldnames(STATS.subject_results);
            
            numfigsA=size(options.FactorA,2);
            
            for i=1:numfigsA
                k=1;
                m=1;
                
                for j=1:length(STATS.condnames);
                    
                    % get the condition waveforms
                    if STATS.subject_results.subject_1.factor_A.contrasts(j,options.FactorA(i))==1
                        %plot1st(k,:)=STATS.condwaves_trim(j,:);
                        lab1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.subject_results.subject_1.factor_A.contrasts(j,options.FactorA(i))==-1
                        %plot2nd(m,:)=STATS.condwaves_trim(j,:);
                        lab2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                end
                
                % concatenate, if needed, the legend lables
                lab1st=strjoin(lab1stlist,'+');
                lab2nd=strjoin(lab2ndlist,'+');
                
                switch plottype
                    case 'wave'
                        for q=1:length(subnames);
                            
                            % extract CI
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,:);
                            
                            % extract difference wave
                            plotdiff=STATS.subject_results.(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                            
                            % begin plotting
                            figure(i);
                            subplot(length(subnames),1,q);
                            
                            % plot zeroline
                            plot(STATS.xtimes,zeros(1,length(STATS.xtimes)),'r','LineWidth',1);
                            hold on
                            
                            % plot CI
                            jbfill(STATS.xtimes,CI_up,CI_low,options.CIcol,options.CIcol,1,1);
                            hold on
                            
                            % plot diff wave
                            plot(STATS.xtimes,plotdiff,'Color',options.diffcol);
                            
                            ylim(gca,options.yaxis);
                            grid on
                            
                        end
                        
                    case 'CI_MOE'
                        
                        for q=1:length(subnames);
                            
                            % extract MOE
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,:);
                            CI_diffup=CI_up-CI_low;
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                            
                            % begin plotting
                            figure(i);
                            subplot(length(subnames),1,q);
                            surf([STATS.xtimes;STATS.xtimes], ones(2,length(STATS.xtimes)), [zeros(1,length(CI_diffup));CI_diffup],[test_stat;test_stat], 'LineStyle', 'none', 'FaceColor', 'interp');
                            zlim(gca,options.zaxis);
                            hold on
                            
                            % find where CI includes zero
                            sigvect=ones(1,STATS.numpnts);
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,:);
                            
                            for m=1:STATS.numpnts;
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs));
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(1));
                            Z=zeros(1,length(timevect));
                            scatter3(timevect,zeroplace,Z,10,'s','filled','MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                            
                            view(0,0);
                            hold on
                            grid on
                        end
                        
                    case 'diff'
                        X=linspace(min(STATS.xtimes),max(STATS.xtimes),length(STATS.xtimes));
                        for q=1:length(subnames);
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                            
                            % begin plotting
                            figure(i);
                            subplot(length(subnames),1,q);
                            imagesc(X,[],test_stat(1,:))
                            hold on
                            set(gca,'YTickLabel',[]);
                            
                            % find where CI includes zero
                            sigvect=ones(1,STATS.numpnts);
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,:);
                            
                            for m=1:STATS.numpnts;
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs))+1.5;
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(1));
                            scatter(timevect,zeroplace,10, 'k', 'MarkerFaceColor', 'k')
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                        end
                end
                
                if ~strcmp(plottype,'wave') && ~strcmp(options.caxis,'auto')
                    ha=axes('visible', 'off');
                    set(ha, 'Units', 'Normalized', 'Position', [.9, 0.05, .015, .9]);
                    colorbar('FontSize',15);
                    caxis(options.caxis)
                end
                
                % add a title to the parent figure, and control its position
                h = axes('visible','off');
                title([lab1st,'-',lab2nd],'parent',h,'visible','on','Position',[.5 1.05 0]);
                clear lab1st lab2nd lab1stlist lab2ndlist
            end
        end
        
        % factor B
        if ~isempty(options.FactorB);
            
            % get subject fields
            subnames = fieldnames(STATS.subject_results);
            
            numfigsB=size(options.FactorB,2);
            
            for i=1:numfigsB
                k=1;
                m=1;
                
                for j=1:length(STATS.condnames);
                    
                    % get the condition waveforms
                    if STATS.subject_results.subject_1.factor_B.contrasts(j,options.FactorB(i))==1
                        %plot1st(k,:)=STATS.condwaves_trim(j,:);
                        lab1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.subject_results.subject_1.factor_B.contrasts(j,options.FactorB(i))==-1
                        %plot2nd(m,:)=STATS.condwaves_trim(j,:);
                        lab2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                end
                
                % concatenate, if needed, the legend lables
                lab1st=strjoin(lab1stlist,'+');
                lab2nd=strjoin(lab2ndlist,'+');
                
                switch plottype
                    case 'wave'
                        for q=1:length(subnames);
                            
                            % extract CI
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(2,:);
                            
                            % extract difference wave
                            plotdiff=STATS.subject_results.(subnames{q}).factor_B.test_stat(options.FactorB(i),:);
                            
                            % begin plotting
                            figure(numfigsA+i);
                            subplot(length(subnames),1,q);
                            
                            % plot zeroline
                            plot(STATS.xtimes,zeros(1,length(STATS.xtimes)),'r','LineWidth',1);
                            hold on
                            
                            % plot CI
                            jbfill(STATS.xtimes,CI_up,CI_low,options.CIcol,options.CIcol,1,1);
                            hold on
                            
                            % plot diff wave
                            plot(STATS.xtimes,plotdiff,'Color',options.diffcol);
                            
                            ylim(gca, options.yaxis);
                            grid on
                            
                        end
                        
                    case 'CI_MOE'
                        
                        
                        for q=1:length(subnames);
                            
                            % extract MOE
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(2,:);
                            CI_diffup=CI_up-CI_low;
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_B.test_stat(options.FactorB(i),:);
                            
                            % begin plotting
                            figure(numfigsA+i);
                            subplot(length(subnames),1,q);
                            surf([STATS.xtimes;STATS.xtimes], ones(2,length(STATS.xtimes)), [zeros(1,length(CI_diffup));CI_diffup],[test_stat;test_stat], 'LineStyle', 'none', 'FaceColor', 'interp');
                            zlim(gca,options.zaxis);
                            hold on
                            
                            % find where CI includes zero
                            sigvect=ones(1,STATS.numpnts);
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(2,:);
                            
                            for m=1:STATS.numpnts;
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs));
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(1));
                            Z=zeros(1,length(timevect));
                            scatter3(timevect,zeroplace,Z,10,'s','filled','MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                            
                            view(0,0)
                            hold on
                            grid on
                        end
                        
                    case 'diff'
                        X=linspace(min(STATS.xtimes),max(STATS.xtimes),length(STATS.xtimes));
                        for q=1:length(subnames);
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_B.test_stat(options.FactorB(i),:);
                            
                            % begin plotting
                            figure(numfigsA+i);
                            subplot(length(subnames),1,q);
                            imagesc(X,[],test_stat(1,:))
                            hold on
                            set(gca,'YTickLabel',[]);
                            
                            % find where CI includes zero
                            sigvect=ones(1,STATS.numpnts);
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(2,:);
                            
                            for m=1:STATS.numpnts;
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs))+1.5;
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(1));
                            scatter(timevect,zeroplace,10, 'k', 'MarkerFaceColor', 'k')
                            
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                            
                        end
                end
                
                if ~strcmp(plottype,'wave') && ~strcmp(options.caxis,'auto')
                    ha=axes('visible', 'off');
                    set(ha, 'Units', 'Normalized', 'Position', [.9, 0.05, .015, .9]);
                    colorbar('FontSize',15);
                    caxis(options.caxis)
                end
                
                % add a title to the parent figure, and control its position
                h = axes('visible','off');
                title([lab1st,'-',lab2nd],'parent',h,'visible','on','Position',[.5 1.05 0]);
                clear lab1st lab2nd lab1stlist lab2ndlist
            end
        end
        % factor AB
        if ~isempty(options.FactorAB);
            
            % get subject fields
            subnames = fieldnames(STATS.subject_results);
            
            numfigsAB=size(options.FactorAB,2);
            
            for i=1:numfigsAB
                k=1;
                m=1;
                
                for j=1:length(STATS.condnames);
                    
                    % get condition waveforms and labels
                    if STATS.inferential_results.factor_AxB.contrasts(j,options.FactorAB(i))==1
                        %plot1st(k,:)=STATS.condwaves_trim(j,:);
                        lab1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.inferential_results.factor_AxB.contrasts(j,options.FactorAB(i))==-1
                        %plot2nd(m,:)=STATS.condwaves_trim(j,:);
                        lab2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                end
                
                % concatenate, if needed, the legend lables
                lab1st=strjoin({lab1stlist{1} lab2ndlist{1}},'-');
                lab2nd=strjoin({lab2ndlist{2} lab1stlist{2}} ,'-');
                
                switch plottype
                    case 'wave'
                        for q=1:length(subnames);
                            
                            % extract CI
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_AxB.CI{options.FactorAB(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_AxB.CI{options.FactorAB(i),1}(2,:);
                            
                            % extract difference wave
                            plotdiff=STATS.subject_results.(subnames{q}).factor_AxB.test_stat(options.FactorAB(i),:);
                            %test_stat=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                            
                            % begin plotting
                            figure(numfigsA+numfigsB+i);
                            subplot(length(subnames),1,q);
                            
                            % plot zeroline
                            plot(STATS.xtimes,zeros(1,length(STATS.xtimes)),'r','LineWidth',1);
                            hold on
                            
                            % plot CI
                            jbfill(STATS.xtimes,CI_up,CI_low,options.CIcol,options.CIcol,1,1);
                            hold on
                            
                            % plot diff wave
                            plot(STATS.xtimes,plotdiff,'Color',options.diffcol);
                            
                            ylim(gca, options.yaxis);
                            grid on
                            
                        end
                        
                    case 'CI_MOE'
                        
                        
                        for q=1:length(subnames);
                            
                            % extract MOE
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_AxB.CI{options.FactorAB(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_AxB.CI{options.FactorAB(i),1}(2,:);
                            CI_diffup=CI_up-CI_low;
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_AxB.test_stat(options.FactorAB(i),:);
                            
                            % begin plotting
                            figure(numfigsA+numfigsB+i);
                            subplot(length(subnames),1,q);
                            surf([STATS.xtimes;STATS.xtimes], ones(2,length(STATS.xtimes)), [zeros(1,length(CI_diffup));CI_diffup],[test_stat;test_stat], 'LineStyle', 'none', 'FaceColor', 'interp');
                            zlim(gca,options.zaxis);
                            hold on
                            
                            % find where CI includes zero
                            sigvect=ones(1,STATS.numpnts);
                            
                            for m=1:STATS.numpnts;
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs));
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(1));
                            Z=zeros(1,length(timevect));
                            scatter3(timevect,zeroplace,Z,10,'s','filled','MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                            
                            view(0,0)
                            hold on
                            grid on
                            
                        end
                        
                    case 'diff'
                        
                        X=linspace(min(STATS.xtimes),max(STATS.xtimes),length(STATS.xtimes));
                        for q=1:length(subnames);
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_AxB.test_stat(options.FactorAB(i),:);
                            
                            % begin plotting
                            figure(numfigsA+numfigsB+i);
                            subplot(length(subnames),1,q);
                            imagesc(X,[],test_stat(1,:))
                            hold on
                            set(gca,'YTickLabel',[]);
                            
                            % find where CI includes zero
                            sigvect=ones(1,STATS.numpnts);
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_AxB.CI{options.FactorAB(i),1}(1,:);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_AxB.CI{options.FactorAB(i),1}(2,:);
                            
                            for m=1:STATS.numpnts;
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs))+1.5;
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(1));
                            scatter(timevect,zeroplace,10, 'k', 'MarkerFaceColor', 'k')
                            
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                        end
                end
                
                if ~strcmp(plottype,'wave') && ~strcmp(options.caxis,'auto')
                    ha=axes('visible', 'off');
                    set(ha, 'Units', 'Normalized', 'Position', [.9, 0.05, .015, .9]);
                    colorbar('FontSize',15);
                    caxis(options.caxis)
                end
                
                % add a title to the parent figure, and control its position
                h = axes('visible','off');
                title(['(',lab1st,')','-','(',lab2nd,')'],'parent',h,'visible','on','Position',[.5 1.05 0]);
                clear lab1st lab2nd lab1stlist lab2ndlist
            end
        end
        
    case 'bw' % do a series of 1-way figures
        
        if ~isempty(options.FactorA);
            
            % get factor fields
            factnames = fieldnames(STATS.subject_results);
            
            % loop through levels of factor A
            y=1;
            for v=1:length(factnames);
                
                % get subject fields
                subnames = fieldnames(STATS.subject_results.(factnames{v}));
                
                numfigs=size(options.FactorA,2);
                
                for i=1:numfigs
                    k=1;
                    m=1;
                    
                    for j=1:length(STATS.klabels);
                        
                        % get the condition waveforms
                        if STATS.subject_results.(factnames{v}).subject_1.factor_A.contrasts(j,options.FactorA(i))==1
                            %plot1st(k,:)=STATS.condwaves_trim(j,:);
                            lab1stlist{k}=STATS.klabels{j};
                            k=k+1;
                        elseif STATS.subject_results.(factnames{v}).subject_1.factor_A.contrasts(j,options.FactorA(i))==-1
                            %plot2nd(m,:)=STATS.condwaves_trim(j,:);
                            lab2ndlist{m}=STATS.klabels{j};
                            m=m+1;
                        end
                        
                    end
                    
                    
                    % concatenate, if needed, the legend lables
                    lab1st=strjoin(lab1stlist,'+');
                    lab2nd=strjoin(lab2ndlist,'+');
                    titlestr=[STATS.jlabels{v},' @ ',lab1st,' - ',lab2nd];
                    
                    switch plottype
                        
                        case 'wave'
                            
                            for q=1:length(subnames);
                                
                                % extract CI
                                CI_low(1,:)=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,:);
                                CI_up(1,:)=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,:);
                                %CI_diffup=CI_up-CI_low;
                                %CI_difflow=CI_diffup*-1;
                                
                                % extract difference wave
                                plotdiff=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                                %test_stat=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                                
                                % begin plotting
                                figure(y);
                                subplot(length(subnames),1,q);
                                
                                % plot zeroline
                                plot(STATS.xtimes,zeros(1,length(STATS.xtimes)),'r','LineWidth',1);
                                hold on
                                
                                % plot CI
                                jbfill(STATS.xtimes,CI_up,CI_low,options.CIcol,options.CIcol,1,1);
                                hold on
                                
                                % plot diff wave
                                plot(STATS.xtimes,plotdiff,'Color',options.diffcol);
                                ylim(gca,options.yaxis);
                                grid on
                                
                            end
                            
                        case 'CI_MOE'
                            
                            
                            for q=1:length(subnames);
                                
                                % extract MOE
                                CI_low(1,:)=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,:);
                                CI_up(1,:)=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,:);
                                CI_diffup=(CI_up-CI_low)/2;
                                
                                % extract difference wave
                                test_stat=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                                
                                % begin plotting
                                figure(y);
                                subplot(length(subnames),1,q);
                                surf([STATS.xtimes;STATS.xtimes], ones(2,length(STATS.xtimes)), [zeros(1,length(CI_diffup));CI_diffup],[test_stat;test_stat], 'LineStyle', 'none', 'FaceColor', 'interp');
                                zlim(gca,options.zaxis);
                                hold on
                                
                                % find where CI includes zero
                                sigvect=ones(1,STATS.numpnts);
                                
                                for m=1:STATS.numpnts;
                                    
                                    if CI_low(m)<=0 && CI_up(m)>=0
                                        sigvect(m)=0;
                                    end
                                    
                                end
                                
                                timelocs=find(sigvect>0);
                                zeroplace=zeros(1,length(timelocs));
                                timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(1));
                                Z=zeros(1,length(timevect));
                                scatter3(timevect,zeroplace,Z,10,'s','filled','MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
                                
                                % give everyone a colorbar
                                if strcmp(options.caxis,'auto');
                                    colorbar;
                                    
                                else
                                    caxis(options.caxis);
                                end
                                
                                view(0,0)
                                hold on
                                grid on
                            end
                            
                        case 'diff'
                            
                            X=linspace(min(STATS.xtimes),max(STATS.xtimes),length(STATS.xtimes));
                            for q=1:length(subnames);
                                
                                % extract difference wave
                                test_stat=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                                
                                % begin plotting
                                figure(y);
                                subplot(length(subnames),1,q);
                                imagesc(X,[],test_stat(1,:))
                                hold on
                                set(gca,'YTickLabel',[]);
                                
                                % find where CI includes zero
                                sigvect=ones(1,STATS.numpnts);
                                CI_low(1,:)=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,:);
                                CI_up(1,:)=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,:);
                                
                                for m=1:STATS.numpnts;
                                    
                                    if CI_low(m)<=0 && CI_up(m)>=0
                                        sigvect(m)=0;
                                    end
                                    
                                end
                                
                                timelocs=find(sigvect>0);
                                zeroplace=zeros(1,length(timelocs))+1.5;
                                timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(1));
                                scatter(timevect,zeroplace,10, 'k', 'MarkerFaceColor', 'k')
                                
                                
                                
                                % give everyone a colorbar
                                if strcmp(options.caxis,'auto');
                                    colorbar;
                                    
                                else
                                    caxis(options.caxis);
                                end
                                
                            end
                    end
                    
                    % add a title and colorbar to the parent figure, and control its position
                    
                    if ~strcmp(plottype,'wave') && ~strcmp(options.caxis,'auto')
                        ha=axes('visible', 'off');
                        set(ha, 'Units', 'Normalized', 'Position', [.9, 0.05, .015, .9]);
                        colorbar('FontSize',15);
                        caxis(options.caxis)
                    end
                    
                    h = axes('visible','off');
                    ht=title(titlestr,'parent',h,'visible','on','Position',[.5 1.05 0]);
                    set(ht,'Interpreter', 'none');
                    y=y+1;
                    clear lab1st lab2nd titlestr lab1stlist lab2ndlist
                end
                
                
            end
            
        end
end

disp('******* Saving STATS structure *******')
save(['STATS_',STATS.savestring,'.mat'],'STATS');


end









