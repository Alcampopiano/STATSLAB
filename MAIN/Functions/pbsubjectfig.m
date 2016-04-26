function [STATS]=pbsubjectfig(STATS,infodisplay,varargin)



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
    options.plottype='CI_MOE';
    options.diffcol=[0 0 0];
    options.CIcol=[.5 .5 .5];
    options.yaxis='auto';
    options.caxis='auto';
    options.zaxis=[0 inf];
    options.timeplot=1:length(STATS.xtimes);
    options.topos='no';
    
    if any(strcmp(varargin,'timeplot'));
        timems=find(strcmp(varargin,'timeplot'));
        ms_input=[varargin{timems+1}(1) varargin{timems+1}(2)]; % get the ms values for setting xaxis later on
        MStoTF_min=round((varargin{timems+1}(1)/1000-STATS.xmin)/(STATS.xmax-STATS.xmin) * (STATS.numpnts-1))+1;
        MStoTF_max=round((varargin{timems+1}(2)/1000-STATS.xmin)/(STATS.xmax-STATS.xmin) * (STATS.numpnts-1))+1;
        varargin{timems+1}=MStoTF_min:MStoTF_max;  
        
    else
        ms_input(1)=min(STATS.xtimes);
        ms_input(2)=max(STATS.xtimes);    
    end
    
    
    
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
    options.plottype='CI_MOE';
    options.diffcol=[0 0 0];
    options.CIcol=[.5 .5 .5];
    options.yaxis='auto';
    options.caxis='auto';
    options.zaxis=[0 inf];
    options.timeplot=1:length(STATS.xtimes);
    options.topos='no';
    
    if any(strcmp(varargin,'timeplot'));
        timems=find(strcmp(varargin,'timeplot'));
        ms_input=[varargin{timems+1}(1) varargin{timems+1}(2)]; % get the ms values for setting xaxis later on
        MStoTF_min=round((varargin{timems+1}(1)/1000-STATS.xmin)/(STATS.xmax-STATS.xmin) * (STATS.numpnts-1))+1;
        MStoTF_max=round((varargin{timems+1}(2)/1000-STATS.xmin)/(STATS.xmax-STATS.xmin) * (STATS.numpnts-1))+1;
        varargin{timems+1}=MStoTF_min:MStoTF_max;
        
    else
        ms_input(1)=min(STATS.xtimes);
        ms_input(2)=max(STATS.xtimes);
    end
    
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

% interp for later topographies or not
if ~strcmp(options.topos, 'no');
    [STATS]=topobuild(STATS,options.topos, 'subject');
    disp('*** finished interpolating and computing topography files. This only needs to be done once ***');
    save(['STATS_',STATS.savestring,'.mat'],'STATS');
end

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
                lab1st=strjoin_statslab(lab1stlist,'+');
                lab2nd=strjoin_statslab(lab2ndlist,'+');
                
                switch options.plottype
                    case 'wave'
                        for q=1:length(subnames);
                            
                            % extract CI
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,options.timeplot);
                            
                            % extract difference wave
                            plotdiff=STATS.subject_results.(subnames{q}).factor_A.test_stat(options.FactorA(i),options.timeplot);
                            
                            % begin plotting
                            figure(i);
                            subplot(length(subnames),1,q);
                            
                            % plot zeroline
                            plot(STATS.xtimes(options.timeplot),zeros(1,length(STATS.xtimes(options.timeplot))),'r','LineWidth',1);
                            hold on
                            
                            % plot CI
                            jbfill(STATS.xtimes(options.timeplot),CI_up,CI_low,options.CIcol,options.CIcol,1,1);
                            hold on
                            
                            % plot diff wave
                            plot(STATS.xtimes(options.timeplot),plotdiff,'Color',options.diffcol);
                            axis tight
                            ylim(gca,options.yaxis);
                            grid on
                            
                            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q});
                            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q});
                            
                        end
                        
                    case 'CI_MOE'
                        
                        
                        for q=1:length(subnames);
                            
                            % extract MOE
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,options.timeplot);
                            CI_diffup=CI_up-CI_low;
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_A.test_stat(options.FactorA(i),options.timeplot);
                            
                            % begin plotting
                            figure(i);
                            subplot(length(subnames),1,q);
                            surf([STATS.xtimes(options.timeplot);STATS.xtimes(options.timeplot)], ones(2,length(STATS.xtimes(options.timeplot))), [zeros(1,length(CI_diffup));CI_diffup],[test_stat;test_stat], 'LineStyle', 'none', 'FaceColor', 'interp');
                            zlim(gca,options.zaxis);
                            hold on
                            
                            % find where CI includes zero
                            sigvect=ones(1,length(STATS.xtimes(options.timeplot)));
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,options.timeplot);
                            
                            for m=1:length(sigvect);
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs));
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(options.timeplot(1)));
                            Z=zeros(1,length(timevect));
                            scatter3(timevect,zeroplace,Z,10,'s','filled','MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
                            
                            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q});
                            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q});
                            
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                            
                            view(0,0)
                            xlim([ms_input(1) ms_input(2)]); % took these values at the start so that I could force the proper Xlim 
                            hold on
                            grid on
                        end
                        
                    case 'diff'
                        
                        X=linspace(min(STATS.xtimes(options.timeplot)),max(STATS.xtimes(options.timeplot)),length(STATS.xtimes(options.timeplot)));
                        for q=1:length(subnames);
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_A.test_stat(options.FactorA(i),options.timeplot);
                            
                            % begin plotting
                            figure(i);
                            subplot(length(subnames),1,q);
                            imagesc(X,[],test_stat(1,:))
                            hold on
                            set(gca,'YTickLabel',[]);
                            
                            % find where CI includes zero
                            sigvect=ones(1,length(STATS.xtimes(options.timeplot)));
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,options.timeplot);
                            
                            for m=1:length(sigvect);
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs))+1.5;
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(options.timeplot(1)));
                            scatter(timevect,zeroplace,10, 'k', 'MarkerFaceColor', 'k')
                            
                            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q});
                            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q});
                            
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                            
                        end
                end
                
                if ~strcmp(options.plottype,'wave') && ~strcmp(options.caxis,'auto')
                    ha=axes('visible', 'off');
                    set(ha, 'Units', 'Normalized', 'Position', [.9, 0.05, .015, .9]);
                    colorbar('FontSize',10);
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
                lab1st=strjoin_statslab(lab1stlist,'+');
                lab2nd=strjoin_statslab(lab2ndlist,'+');
                
                switch options.plottype
                    case 'wave'
                        for q=1:length(subnames);
                            
                            % extract CI
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,options.timeplot);
                            
                            % extract difference wave
                            plotdiff=STATS.subject_results.(subnames{q}).factor_A.test_stat(options.FactorA(i),options.timeplot);
                            
                            % begin plotting
                            figure(i);
                            subplot(length(subnames),1,q);
                           
                            
                            % plot zeroline
                            plot(STATS.xtimes(options.timeplot),zeros(1,length(STATS.xtimes(options.timeplot))),'r','LineWidth',1);
                            hold on
                            
                            % plot CI
                            jbfill(STATS.xtimes(options.timeplot),CI_up,CI_low,options.CIcol,options.CIcol,1,1);
                            hold on
                            
                            % plot diff wave
                            plot(STATS.xtimes(options.timeplot),plotdiff,'Color',options.diffcol);
                            
                            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q});
                            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q});
                            
                            axis tight
                            ylim(gca,options.yaxis);
                            grid on
                            
                        end
                        
                    case 'CI_MOE'
                        
                        for q=1:length(subnames);
                            
                            % extract MOE
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,options.timeplot);
                            CI_diffup=CI_up-CI_low;
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_A.test_stat(options.FactorA(i),options.timeplot);
                            
                            % begin plotting
                            figure(i);
                            subplot(length(subnames),1,q);
                            surf([STATS.xtimes(options.timeplot);STATS.xtimes(options.timeplot)], ones(2,length(STATS.xtimes(options.timeplot))), [zeros(1,length(CI_diffup));CI_diffup],[test_stat;test_stat], 'LineStyle', 'none', 'FaceColor', 'interp');
                            axis tight
                            zlim(gca,options.zaxis);
                            hold on
                            
                            % find where CI includes zero
                            sigvect=ones(1,length(STATS.xtimes(options.timeplot)));
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,options.timeplot);
                            
                            for m=1:length(sigvect);
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs));
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(options.timeplot(1)));
                            Z=zeros(1,length(timevect));
                            scatter3(timevect,zeroplace,Z,10,'s','filled','MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
                            
                            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q});
                            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q});
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                            
                            view(0,0);
                            xlim([ms_input(1) ms_input(2)]); % took these values at the start so that I could force the proper Xlim 
                            hold on
                            grid on
                        end
                        
                    case 'diff'
                        X=linspace(min(STATS.xtimes(options.timeplot)),max(STATS.xtimes(options.timeplot)),length(STATS.xtimes(options.timeplot)));
                        for q=1:length(subnames);
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_A.test_stat(options.FactorA(i),options.timeplot);
                            
                            % begin plotting
                            figure(i);
                            subplot(length(subnames),1,q);
                            imagesc(X,[],test_stat(1,:))
                            hold on
                            set(gca,'YTickLabel',[]);
                            
                            % find where CI includes zero
                            sigvect=ones(1,length(STATS.xtimes(options.timeplot)));
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,options.timeplot);
                            
                            for m=1:length(sigvect);
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs))+1.5;
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(options.timeplot(1)));
                            scatter(timevect,zeroplace,10, 'k', 'MarkerFaceColor', 'k')
                            
                            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q});
                            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q});
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                        end
                end
                
                if ~strcmp(options.plottype,'wave') && ~strcmp(options.caxis,'auto')
                    ha=axes('visible', 'off');
                    set(ha, 'Units', 'Normalized', 'Position', [.9, 0.05, .015, .9]);
                    colorbar('FontSize',10);
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
                lab1st=strjoin_statslab(lab1stlist,'+');
                lab2nd=strjoin_statslab(lab2ndlist,'+');
                
                switch options.plottype
                    case 'wave'
                        for q=1:length(subnames);
                            
                            % extract CI
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(2,options.timeplot);
                            
                            % extract difference wave
                            plotdiff=STATS.subject_results.(subnames{q}).factor_B.test_stat(options.FactorB(i),options.timeplot);
                            
                            % begin plotting
                            figure(numfigsA+i);
                            subplot(length(subnames),1,q);
                            
                            % plot zeroline
                            plot(STATS.xtimes(options.timeplot),zeros(1,length(STATS.xtimes(options.timeplot))),'r','LineWidth',1);
                            axis tight
                            hold on
                            
                            % plot CI
                            jbfill(STATS.xtimes(options.timeplot),CI_up,CI_low,options.CIcol,options.CIcol,1,1);
                            hold on
                            
                            % plot diff wave
                            plot(STATS.xtimes(options.timeplot),plotdiff,'Color',options.diffcol);
                            
                            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q});
                            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q});
                            
                            ylim(gca, options.yaxis);
                            grid on
                            
                        end
                        
                    case 'CI_MOE'
                        
                        
                        for q=1:length(subnames);
                            
                            % extract MOE
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(2,options.timeplot);
                            CI_diffup=CI_up-CI_low;
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_B.test_stat(options.FactorB(i),options.timeplot);
                            
                            % begin plotting
                            figure(numfigsA+i);
                            subplot(length(subnames),1,q);
                            surf([STATS.xtimes(options.timeplot);STATS.xtimes(options.timeplot)], ones(2,length(STATS.xtimes(options.timeplot))), [zeros(1,length(CI_diffup));CI_diffup],[test_stat;test_stat], 'LineStyle', 'none', 'FaceColor', 'interp');
                            axis tight
                            zlim(gca,options.zaxis);
                            hold on
                            
                            % find where CI includes zero
                            sigvect=ones(1,length(STATS.xtimes(options.timeplot)));
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(2,options.timeplot);
                            
                            for m=1:length(sigvect);
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs));
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(options.timeplot(1)));
                            Z=zeros(1,length(timevect));
                            scatter3(timevect,zeroplace,Z,10,'s','filled','MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
                            
                            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q});
                            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q});
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                            
                            view(0,0)
                            xlim([ms_input(1) ms_input(2)]); % took these values at the start so that I could force the proper Xlim 
                            hold on
                            grid on
                        end
                        
                    case 'diff'
                        X=linspace(min(STATS.xtimes(options.timeplot)),max(STATS.xtimes(options.timeplot)),length(STATS.xtimes(options.timeplot)));
                        for q=1:length(subnames);
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_B.test_stat(options.FactorB(i),options.timeplot);
                            
                            % begin plotting
                            figure(numfigsA+i);
                            subplot(length(subnames),1,q);
                            imagesc(X,[],test_stat(1,:))
                            hold on
                            set(gca,'YTickLabel',[]);
                            
                            % find where CI includes zero
                            sigvect=ones(1,length(STATS.xtimes(options.timeplot)));
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_B.CI{options.FactorB(i),1}(2,options.timeplot);
                            
                            for m=1:length(sigvect);
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs))+1.5;
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(options.timeplot(1)));
                            scatter(timevect,zeroplace,10, 'k', 'MarkerFaceColor', 'k')
                            
                            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q});
                            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q});
                            
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                            
                        end
                end
                
                if ~strcmp(options.plottype,'wave') && ~strcmp(options.caxis,'auto')
                    ha=axes('visible', 'off');
                    set(ha, 'Units', 'Normalized', 'Position', [.9, 0.05, .015, .9]);
                    colorbar('FontSize',10);
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
                    if STATS.sample_results.factor_AxB.contrasts(j,options.FactorAB(i))==1
                        %plot1st(k,:)=STATS.condwaves_trim(j,:);
                        lab1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.sample_results.factor_AxB.contrasts(j,options.FactorAB(i))==-1
                        %plot2nd(m,:)=STATS.condwaves_trim(j,:);
                        lab2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                end
                
                % concatenate, if needed, the legend lables
                lab1st=strjoin_statslab({lab1stlist{1} lab2ndlist{1}},'-');
                lab2nd=strjoin_statslab({lab2ndlist{2} lab1stlist{2}} ,'-');
                
                switch options.plottype
                    case 'wave'
                        for q=1:length(subnames);
                            
                            % extract CI
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_AxB.CI{options.FactorAB(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_AxB.CI{options.FactorAB(i),1}(2,options.timeplot);
                            
                            % extract difference wave
                            plotdiff=STATS.subject_results.(subnames{q}).factor_AxB.test_stat(options.FactorAB(i),options.timeplot);
                            %test_stat=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                            
                            % begin plotting
                            figure(numfigsA+numfigsB+i);
                            subplot(length(subnames),1,q);
                            
                            % plot zeroline
                            plot(STATS.xtimes(options.timeplot),zeros(1,length(STATS.xtimes(options.timeplot))),'r','LineWidth',1);
                            axis tight
                            hold on
                            
                            % plot CI
                            jbfill(STATS.xtimes(options.timeplot),CI_up,CI_low,options.CIcol,options.CIcol,1,1);
                            hold on
                            
                            % plot diff wave
                            plot(STATS.xtimes(options.timeplot),plotdiff,'Color',options.diffcol);
                            
                            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q});
                            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q});
                            
                            ylim(gca, options.yaxis);
                            grid on
                            
                        end
                        
                    case 'CI_MOE'
                        
                        
                        for q=1:length(subnames);
                            
                            % extract MOE
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_AxB.CI{options.FactorAB(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_AxB.CI{options.FactorAB(i),1}(2,options.timeplot);
                            CI_diffup=CI_up-CI_low;
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_AxB.test_stat(options.FactorAB(i),options.timeplot);
                            
                            % begin plotting
                            figure(numfigsA+numfigsB+i);
                            subplot(length(subnames),1,q);
                            surf([STATS.xtimes(options.timeplot);STATS.xtimes(options.timeplot)], ones(2,length(STATS.xtimes(options.timeplot))), [zeros(1,length(CI_diffup));CI_diffup],[test_stat;test_stat], 'LineStyle', 'none', 'FaceColor', 'interp');
                            axis tight
                            zlim(gca,options.zaxis);
                            hold on
                            
                            % find where CI includes zero
                            sigvect=ones(1,length(STATS.xtimes(options.timeplot)));
                            
                            for m=1:length(sigvect);
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs));
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(options.timeplot(1)));
                            Z=zeros(1,length(timevect));
                            scatter3(timevect,zeroplace,Z,10,'s','filled','MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
                            
                            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q});
                            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q});
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                            
                            view(0,0)
                            xlim([ms_input(1) ms_input(2)]); % took these values at the start so that I could force the proper Xlim 
                            hold on
                            grid on
                            
                        end
                        
                    case 'diff'
                        
                        X=linspace(min(STATS.xtimes(options.timeplot)),max(STATS.xtimes(options.timeplot)),length(STATS.xtimes(options.timeplot)));
                        for q=1:length(subnames);
                            
                            % extract difference wave
                            test_stat=STATS.subject_results.(subnames{q}).factor_AxB.test_stat(options.FactorAB(i),options.timeplot);
                            
                            % begin plotting
                            figure(numfigsA+numfigsB+i);
                            subplot(length(subnames),1,q);
                            imagesc(X,[],test_stat(1,:))
                            hold on
                            set(gca,'YTickLabel',[]);
                            
                            % find where CI includes zero
                            sigvect=ones(1,length(STATS.xtimes(options.timeplot)));
                            CI_low(1,:)=STATS.subject_results.(subnames{q}).factor_AxB.CI{options.FactorAB(i),1}(1,options.timeplot);
                            CI_up(1,:)=STATS.subject_results.(subnames{q}).factor_AxB.CI{options.FactorAB(i),1}(2,options.timeplot);
                            
                            for m=1:length(sigvect);
                                
                                if CI_low(m)<=0 && CI_up(m)>=0
                                    sigvect(m)=0;
                                end
                                
                            end
                            
                            timelocs=find(sigvect>0);
                            zeroplace=zeros(1,length(timelocs))+1.5;
                            timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(options.timeplot(1)));
                            scatter(timevect,zeroplace,10, 'k', 'MarkerFaceColor', 'k')
                            
                            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q});
                            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q});
                            
                            
                            % give everyone a colorbar
                            if strcmp(options.caxis,'auto');
                                colorbar;
                                
                            else
                                caxis(options.caxis);
                            end
                        end
                end
                
                if ~strcmp(options.plottype,'wave') && ~strcmp(options.caxis,'auto')
                    ha=axes('visible', 'off');
                    set(ha, 'Units', 'Normalized', 'Position', [.9, 0.05, .015, .9]);
                    colorbar('FontSize',10);
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
                    lab1st=strjoin_statslab(lab1stlist,'+');
                    lab2nd=strjoin_statslab(lab2ndlist,'+');
                    titlestr=[STATS.jlabels{v},' @ ',lab1st,' - ',lab2nd];
                    
                    switch options.plottype
                        
                        case 'wave'
                            
                            for q=1:length(subnames);
                                
                                % extract CI
                                CI_low(1,:)=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,options.timeplot);
                                CI_up(1,:)=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,options.timeplot);
                                %CI_diffup=CI_up-CI_low;
                                %CI_difflow=CI_diffup*-1;
                                
                                % extract difference wave
                                plotdiff=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.test_stat(options.FactorA(i),options.timeplot);
                                %test_stat=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.test_stat(options.FactorA(i),:);
                                
                                % begin plotting
                                figure(y);
                                subplot(length(subnames),1,q);
                                
                                % plot zeroline
                                plot(STATS.xtimes(options.timeplot),zeros(1,length(STATS.xtimes(options.timeplot))),'r','LineWidth',1);
                                axis tight
                                hold on
                                
                                % plot CI
                                jbfill(STATS.xtimes(options.timeplot),CI_up,CI_low,options.CIcol,options.CIcol,1,1);
                                hold on
                                
                                % plot diff wave
                                plot(STATS.xtimes(options.timeplot),plotdiff,'Color',options.diffcol);
                                ylim(gca,options.yaxis);
                                grid on
                                
                                set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q,v});
                                set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q,v});
                                
                            end
                            
                        case 'CI_MOE'
                            
                            
                            for q=1:length(subnames);
                                
                                % extract MOE
                                CI_low(1,:)=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,options.timeplot);
                                CI_up(1,:)=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,options.timeplot);
                                CI_diffup=(CI_up-CI_low)/2;
                                
                                % extract difference wave
                                test_stat=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.test_stat(options.FactorA(i),options.timeplot);
                                
                                % begin plotting
                                figure(y);
                                subplot(length(subnames),1,q);
                                surf([STATS.xtimes(options.timeplot);STATS.xtimes(options.timeplot)], ones(2,length(STATS.xtimes(options.timeplot))), [zeros(1,length(CI_diffup));CI_diffup],[test_stat;test_stat], 'LineStyle', 'none', 'FaceColor', 'interp');
                                axis tight
                                zlim(gca,options.zaxis);
                                hold on
                                
                                % find where CI includes zero
                                sigvect=ones(1,length(STATS.xtimes(options.timeplot)));
                                
                                for m=1:length(sigvect);
                                    
                                    if CI_low(m)<=0 && CI_up(m)>=0
                                        sigvect(m)=0;
                                    end
                                    
                                end
                                
                                timelocs=find(sigvect>0);
                                zeroplace=zeros(1,length(timelocs));
                                timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(options.timeplot(1)));
                                Z=zeros(1,length(timevect));
                                scatter3(timevect,zeroplace,Z,10,'s','filled','MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
                                
                                set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q,v});
                                set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q,v});
                                
                                % give everyone a colorbar
                                if strcmp(options.caxis,'auto');
                                    colorbar;
                                    
                                else
                                    caxis(options.caxis);
                                end
                                
                                view(0,0)
                                xlim([ms_input(1) ms_input(2)]); % took these values at the start so that I could force the proper Xlim 
                                hold on
                                grid on
                            end
                            
                        case 'diff'
                            
                            X=linspace(min(STATS.xtimes(options.timeplot)),max(STATS.xtimes(options.timeplot)),length(STATS.xtimes(options.timeplot)));
                            for q=1:length(subnames);
                                
                                % extract difference wave
                                test_stat=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.test_stat(options.FactorA(i),options.timeplot);
                                
                                % begin plotting
                                figure(y);
                                subplot(length(subnames),1,q);
                                imagesc(X,[],test_stat(1,:))
                                hold on
                                set(gca,'YTickLabel',[]);
                                
                                % find where CI includes zero
                                sigvect=ones(1,length(STATS.xtimes(options.timeplot)));
                                CI_low(1,:)=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.CI{options.FactorA(i),1}(1,options.timeplot);
                                CI_up(1,:)=STATS.subject_results.(factnames{v}).(subnames{q}).factor_A.CI{options.FactorA(i),1}(2,options.timeplot);
                                
                                for m=1:length(sigvect);
                                    
                                    if CI_low(m)<=0 && CI_up(m)>=0
                                        sigvect(m)=0;
                                    end
                                    
                                end
                                
                                timelocs=find(sigvect>0);
                                zeroplace=zeros(1,length(timelocs))+1.5;
                                timevect=(timelocs*(1000/STATS.srate))-abs(STATS.xtimes(options.timeplot(1)));
                                scatter(timevect,zeroplace,10, 'k', 'MarkerFaceColor', 'k')
                                
                                set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS,q,v});
                                set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS,q,v});
                                
                                
                                
                                % give everyone a colorbar
                                if strcmp(options.caxis,'auto');
                                    colorbar;
                                    
                                else
                                    caxis(options.caxis);
                                end
                                
                            end
                    end
                    
                    % add a title and colorbar to the parent figure, and control its position
                    
                    if ~strcmp(options.plottype,'wave') && ~strcmp(options.caxis,'auto')
                        ha=axes('visible', 'off');
                        set(ha, 'Units', 'Normalized', 'Position', [.9, 0.05, .015, .9]);
                        colorbar('FontSize',10);
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

function mouseclick_callback(gcbo,eventdata,STATS,sub,atjlvl)
        
        try

            SIZEBOX=250; % some arbitrary size of some box
            
            if strcmp(STATS.design, 'bw')
                
                atcond=(atjlvl-1)*STATS.levels(2)+1; % to determine which cell to start from
                col=atcond+STATS.levels(2)-1; % the last cell to plot
                col_size=STATS.levels(2);
                
                % get dimensions.
                rowcols(2) = ceil(sqrt(col_size)); % EEGpage is number of subjects/topos
                rowcols(1) = ceil(col_size/rowcols(2));
                
            else
                
                [row col]=size(STATS.grouptopofiles);
                atcond=1;
                
                % get dimensions.
                rowcols(2) = ceil(sqrt(col)); % EEGpage is number of subjects/topos
                rowcols(1) = ceil(col/rowcols(2));
                
            end

            
            %get the point that was clicked on
            cP = get(gca,'Currentpoint');
            ms_plot = cP(1,1);
            %y = cP(1,2);
            
            s=1;
            for r=atcond:col % loop for each subject?
                
                if r==atcond;
                    % build eventual destination figure
                    curfig = figure('paperpositionmode', 'auto', 'visible', 'off');
                    pos = get(curfig,'Position');
                    posx = max(0, pos(1)+(pos(3)-SIZEBOX*rowcols(2))/2);
                    posy = pos(2)+pos(4)-SIZEBOX*rowcols(1);
                    set(curfig,'Position', [posx posy  SIZEBOX*rowcols(2)  SIZEBOX*rowcols(1)]);
                    
                end
                curax = subplot( rowcols(1), rowcols(2), mod(r-1, rowcols(1)*rowcols(2))+1);
                set(curax, 'visible', 'off')
                
                % loading subject
                data=load(STATS.subtopofiles{r}{sub});
                disp(STATS.subtopofiles{r}{sub});
                
                MStoTF=round((ms_plot/1000-data.EEG.xmin)/(data.EEG.xmax-data.EEG.xmin) * (data.EEG.pnts-1))+1;
                maplim=max(max(abs(data.EEG.data(:,MStoTF))));
                
                if isempty(maplim);
                    pop_topoplot(data.EEG, 1, ms_plot, [], 0,'shading','interp','colorbar','off');
                    %set(gcf, 'visible','off');
                    htopo(s)=gca;
                    ctopo(s,:)=caxis;
                else
                    pop_topoplot(data.EEG, 1, ms_plot, [], 0,'shading','interp','colorbar','off','maplimits', [-maplim maplim]);
                    %set(gcf, 'visible','off');
                    htopo(s)=gca;
                    ctopo(s,:)=caxis;
                end
                
                oh=findobj(curax); % find and get rid of EEGLABs subplot titles
                alltext=findall(oh,'Type','text');
                delete(alltext);
                text(.5,-.1,num2str(STATS.condnames{r}),'Units','normalized','Interpreter', 'none'); % add subject numbers to bottom centre of subplots
                
                if isempty(maplim)
                    colorbar;
                end

                if ~isempty(maplim)
                    if r==col % last subject
                        
                        maplim=max(max(abs(ctopo)));
                        
                        % set limits
                        for qq=1:length(htopo);
                            caxis(htopo(qq), [-maplim maplim])
                        end
                        
                        
                        hax=axes('visible', 'off');
                        set(hax, 'Units', 'Normalized', 'Position', [.88, 0.25, .025, .5]);
                        colorbar('FontSize',15);
                        caxis([-maplim maplim]);
       
%                         curax_pos=get(curax,'position');
%                         colorbar('location','eastoutside');
%                         set(curax,'position',curax_pos);
                    end
                end
                s=s+1;
                
            end % end of r loop
            htit = axes('visible','off');
            title(['Subject #', num2str(sub), ' at ', num2str(ms_plot),'ms'],'parent',htit,'visible','on');
            
            %title([condlabs{i}, ' from ', num2str(ms_plot),'ms'],'parent',h,'visible','on');
            
            % get and set title handle
            %thandle = get(gca,'Title');
            %set(thandle,'String',s);
            % finally change the position of our red plus, and make it
            % visible.
            %set(cursor_handle,'Xdata',x,'Ydata',y,'visible','on')
        catch
            disp('no topographies available at this time');
        end
        
    end


end









