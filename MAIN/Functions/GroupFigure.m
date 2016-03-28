function [STATS]=GroupFigure(STATS,infodisplay,varargin)

% This function plots group-level data. 
% For ERSP and ITC measures, basic plots are created where the group-level time X frequency for the specified contrast is plotted. 
% All figures have clickable features, so click waveforms/frequencies to display more information (e.g., topographies, CIs around frequency bands).
% 
% Inputs:
% 
% ***infodisplay***
% Check box to display your condition labels and contrast matrices when plotting figures ***end***
% 
% 
% ***varargin***
% Options are specified in pairs (key -> val)
% 
% timeplot ->
% 
% 	[numeric] - the range of xaxis values to plot. Default scales to length of epoch.
% 
% topos ->
% 	 no (default) - do not calculate topographies
% 
% 	channel location file - calculate and allow plotting of topographies (only need to do once).
% 
% FactorA, FactorB, FactorAB, all (use the "all" option without a key/value pair) ->
% 
% 	[numeric] - specify which contrasts you wish to plot for each factor and the interaction
% 
% 	all (default) - plot all contrasts that were specified when statistics were calculated. 
%                   This option does not need to be paired with another option, just use it on its own 
%                   (i.e., do not pair with FactorA, FactorB, or FactorAB).
% 
% 
% For example,
% 
% FactorA
% 2
% timeplot
% -200 600
% topos
% mychanlocsfile.sfp
% 
% This set of options would plot the second Factor A comparison (otherwise just use the “all” option to plot all contrasts). 
% The type of plot would be a standard difference wave in blue for each subject bounded by a CI. 
% The channel locations file would be use to precompute topographies. 
% This allows you to click on the waveforms at any latency for any subject and see the topographies for each condition.
% 
% Using GroupFigure at the commandline:
% 
% [STATS]=GroupFigure_sample('STATS_mythesis.mat',1,'all','timeplot',[-200 600], 'topos' ,'no');
% ***end***


if isempty(STATS)
    [fnamestats]=uigetfile('*.mat','Select the STATS structure for this analysis');
    load(fnamestats);
else 
    load(STATS);
end

% set history
[hist_str]=statslab_history(['STATS_', STATS.savestring, '.mat'],infodisplay,varargin);
STATS.history.GroupFigure=hist_str;

% special case where using 'all' plots all possible contrasts
if length(varargin)==1 && strcmp(varargin{1},'all');
    
    % this finds out if the design was factorial or not and sets default
    % options accordingly
    if size(fieldnames(STATS.inferential_results),1)==3;
        isfactorial=1;
        
        options = struct('FactorA', 1:size(STATS.inferential_results.factor_A.contrasts,2), ...
            'FactorB', 1:size(STATS.inferential_results.factor_B.contrasts,2), ...
            'FactorAB', 1:size(STATS.inferential_results.factor_AxB.contrasts,2));
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.inferential_results.factor_A.contrasts)
            disp('FactorB'); disp(STATS.inferential_results.factor_B.contrasts)
            disp('FactorAB'); disp(STATS.inferential_results.factor_AxB.contrasts)
        end
        
    elseif size(fieldnames(STATS.inferential_results),1)==1;
        isfactorial=0;
        options = struct('FactorA', 1:size(STATS.inferential_results.factor_A.contrasts,2));
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.inferential_results.factor_A.contrasts)
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
    if size(fieldnames(STATS.inferential_results),1)==3;
        isfactorial=1;
        options = struct('FactorA', [],'FactorB', [], 'FactorAB', []);
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.inferential_results.factor_A.contrasts)
            disp('FactorB'); disp(STATS.inferential_results.factor_B.contrasts)
            disp('FactorAB'); disp(STATS.inferential_results.factor_AxB.contrasts)
        end
        
    elseif size(fieldnames(STATS.inferential_results),1)==1;
        isfactorial=0;
        options = struct('FactorA', []);
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.inferential_results.factor_A.contrasts)
        end
        
        
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
STATS.plotoptions=options;

% switch cases for factorial vs. single-factor
switch isfactorial
    case 0
        
        numfigs=size(options.FactorA,2);
        
        for i=1:numfigs
            k=1;
            m=1;
            
            for j=1:length(STATS.condnames);
                
                % get the condition waveforms
                if STATS.inferential_results.factor_A.contrasts(j,options.FactorA(i))==1
                    plot1st(k,:)=STATS.condwaves_trim(j,:);
                    leg1stlist{k}=STATS.condnames{j};
                    k=k+1;
                elseif STATS.inferential_results.factor_A.contrasts(j,options.FactorA(i))==-1
                    plot2nd(m,:)=STATS.condwaves_trim(j,:);
                    leg2ndlist{m}=STATS.condnames{j};
                    m=m+1;
                end
                
            end
            
            
            
            % reduce, if needed, the condition waveforms
            % should this be SUM or MEAN if pooling across levels?
            plot1st=mean(plot1st,1);
            plot2nd=mean(plot2nd,1);
            plotdiff=STATS.inferential_results.factor_A.test_stat(options.FactorA(i),:);
            
            % concatenate, if needed, the legend lables
            leg1st=strjoin_statslab(leg1stlist,'+');
            leg2nd=strjoin_statslab(leg2ndlist,'+');
            
            % begin plotting
            figure(i);
            subplot(2,1,1)
            h(1)=plot(STATS.xtimes,plot1st,'r', 'LineWidth',3);
            hold on
            h(2)=plot(STATS.xtimes,plot2nd,'b','LineWidth',3);
            set(gca,'FontSize',20)
            %title('CS vs CT')
            lh=legend(leg1st,leg2nd);
            % get rid of subscripts that occur when there are underscores
            set(lh,'Interpreter', 'none');
            axis tight
            grid on
            
            subplot(2,1,2)
            h(3)=plot(STATS.xtimes,plotdiff,'k');
            CIup=STATS.inferential_results.factor_A.CI{options.FactorA(i)}(2,:);
            CIlow=STATS.inferential_results.factor_A.CI{options.FactorA(i)}(1,:);
            h(4)=jbfill(STATS.xtimes,CIup, CIlow, [.5 .5 .5], [.5 .5 .5], 1, .6);
            axis tight
            grid on
            
            % clear variables before next iteration
            clear plot1st plot2nd plotdiff leg1st leg2nd leg2ndlist leg1stlist lh
            
            
            
        end
        
        
    case 1
        
        if ~isempty(options.FactorA);
            numfigs=size(options.FactorA,2);
            
            for i=1:numfigs
                k=1;
                m=1;
                
                for j=1:length(STATS.condnames);
                    
                    % get the condition waveforms
                    if STATS.inferential_results.factor_A.contrasts(j,options.FactorA(i))==1
                        plot1st(k,:)=STATS.condwaves_trim(j,:);
                        leg1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.inferential_results.factor_A.contrasts(j,options.FactorA(i))==-1
                        plot2nd(m,:)=STATS.condwaves_trim(j,:);
                        leg2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                end
                
                
                
                % reduce, if needed, the condition waveforms
                % should this be SUM or MEAN if pooling across levels?
                plot1st=mean(plot1st,1);
                plot2nd=mean(plot2nd,1);
                plotdiff=STATS.inferential_results.factor_A.test_stat(options.FactorA(i),:);
                
                % concatenate, if needed, the legend lables
                leg1st=strjoin_statslab(leg1stlist,'+');
                leg2nd=strjoin_statslab(leg2ndlist,'+');
                
                % begin plotting
                figure;
                subplot(2,1,1)
                h(1)=plot(STATS.xtimes,plot1st,'r', 'LineWidth',3);
                hold on
                h(2)=plot(STATS.xtimes,plot2nd,'b','LineWidth',3);
                set(gca,'FontSize',20)
                %title('CS vs CT')
                lh=legend(leg1st,leg2nd);
                % get rid of subscripts that occur when there are underscores
                set(lh,'Interpreter', 'none');
                axis tight
                grid on
                
                subplot(2,1,2)
                h(3)=plot(STATS.xtimes,plotdiff,'k');
                CIup=STATS.inferential_results.factor_A.CI{options.FactorA(i)}(2,:);
                CIlow=STATS.inferential_results.factor_A.CI{options.FactorA(i)}(1,:);
                h(4)=jbfill(STATS.xtimes,CIup, CIlow, [.5 .5 .5], [.5 .5 .5], 1, .6);
                axis tight
                grid on
                
                % clear variables before next iteration
                clear plot1st plot2nd plotdiff leg1st leg2nd leg2ndlist leg1stlist
                
                
                
            end
        end
        
        if ~isempty(options.FactorB);
            numfigs=size(options.FactorB,2);
            
            for i=1:numfigs
                k=1;
                m=1;
                
                for j=1:length(STATS.condnames);
                    
                    % get the condition waveforms
                    if STATS.inferential_results.factor_B.contrasts(j,options.FactorB(i))==1
                        plot1st(k,:)=STATS.condwaves_trim(j,:);
                        leg1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.inferential_results.factor_B.contrasts(j,options.FactorB(i))==-1
                        plot2nd(m,:)=STATS.condwaves_trim(j,:);
                        leg2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                end
                
                
                
                % reduce, if needed, the condition waveforms
                % should this be SUM or MEAN if pooling across levels?
                plot1st=mean(plot1st,1);
                plot2nd=mean(plot2nd,1);
                plotdiff=STATS.inferential_results.factor_B.test_stat(options.FactorB(i),:);
                
                % concatenate, if needed, the legend lables
                leg1st=strjoin_statslab(leg1stlist,'+');
                leg2nd=strjoin_statslab(leg2ndlist,'+');
                
                % begin plotting
                figure;
                subplot(2,1,1)
                h(1)=plot(STATS.xtimes,plot1st,'r', 'LineWidth',3);
                hold on
                h(2)=plot(STATS.xtimes,plot2nd,'b','LineWidth',3);
                set(gca,'FontSize',20)
                %title('CS vs CT')
                lh=legend(leg1st,leg2nd);
                % get rid of subscripts that occur when there are underscores
                set(lh,'Interpreter', 'none');
                axis tight
                grid on
                
                subplot(2,1,2)
                h(3)=plot(STATS.xtimes,plotdiff,'k');
                CIup=STATS.inferential_results.factor_B.CI{options.FactorB(i)}(2,:);
                CIlow=STATS.inferential_results.factor_B.CI{options.FactorB(i)}(1,:);
                h(4)=jbfill(STATS.xtimes,CIup, CIlow, [.5 .5 .5], [.5 .5 .5], 1, .6);
                axis tight
                grid on
                
                % clear variables before next iteration
                clear plot1st plot2nd plotdiff leg1st leg2nd leg2ndlist leg1stlist
                
                
                
            end
        end
        
        if ~isempty(options.FactorAB);
            numfigs=size(options.FactorAB,2);
            
            for i=1:numfigs
                k=1;
                m=1;
                
                for j=1:length(STATS.condnames);
                    
                    % get condition waveforms and labels
                    if STATS.inferential_results.factor_AxB.contrasts(j,options.FactorAB(i))==1
                        plot1st(k,:)=STATS.condwaves_trim(j,:);
                        leg1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.inferential_results.factor_AxB.contrasts(j,options.FactorAB(i))==-1
                        plot2nd(m,:)=STATS.condwaves_trim(j,:);
                        leg2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                
                end

                % create interaction difference waveforms 
                plot1stdiff=plot1st(1,:)-plot2nd(1,:);
                plot2nddiff=plot2nd(2,:)-plot1st(2,:);
                plotdiff=STATS.inferential_results.factor_AxB.test_stat(options.FactorAB(i),:);
                
                % concatenate, if needed, the legend lables
                leg1st=strjoin_statslab({leg1stlist{1} leg2ndlist{1}},'-');
                leg2nd=strjoin_statslab({leg2ndlist{2} leg1stlist{2}} ,'-');
                
                % begin plotting
                figure;
                subplot(2,1,1)
                h(1)=plot(STATS.xtimes,plot1stdiff,'r', 'LineWidth',3);
                hold on
                h(2)=plot(STATS.xtimes,plot2nddiff,'b','LineWidth',3);
                set(gca,'FontSize',20)
                %title('CS vs CT')
                lh=legend(leg1st,leg2nd);
                % get rid of subscripts that occur when there are underscores
                set(lh,'Interpreter', 'none');
                axis tight
                grid on
                
                subplot(2,1,2)
                h(3)=plot(STATS.xtimes,plotdiff,'k');
                CIup=STATS.inferential_results.factor_AxB.CI{options.FactorAB(i)}(2,:);
                CIlow=STATS.inferential_results.factor_AxB.CI{options.FactorAB(i)}(1,:);
                h(4)=jbfill(STATS.xtimes,CIup, CIlow, [.5 .5 .5], [.5 .5 .5], 1, .6);
                axis tight
                grid on
                
                % clear variables before next iteration
                clear plot1st plot2nd plotdiff leg1st leg2nd leg2ndlist leg1stlist
                
                
                
            end
        end
   
end

disp('******* Saving STATS structure *******')
save(['STATS_',STATS.savestring,'.mat'],'STATS');









end









