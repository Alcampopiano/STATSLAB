function [STATS]=pbgroupfig_sample(STATS,infodisplay,varargin)
% build figure for group results
% varargin can be 'all' which plots all contrasts for all factors, or
% specify what you want with name/val pairs:
% ...'FactorA', [1 2 3], 'FactorB', [1], 'FactorAB', [1])
% Defaults for contrasts to use are empty, so ommiting one of the factors will not produce any plot for that one.
% infodisplay is a flag 1 or 0 that will spit out the contrast matrices and condition lables for you to look at in the command window
% infodisplay default is set to 0.

% 
% if isempty(STATS)
%     [fnamestats]=uigetfile('*.mat','Select the STATS structure for this analysis');
%     load(fnamestats);
% else 
%     load(STATS);
% end

% set history
% [hist_str]=statslab_history(['STATS_', STATS.savestring, '.mat'],infodisplay,varargin);
% STATS.history.GroupFigure=hist_str;

% special case where using 'all' plots all possible contrasts
if any(strcmp(varargin,'all')); 
    
    % something silly to make varargin pairable
    % varar=cell(1,2);
    % varar{1}=varargin{1};
    % varargin=varar;
    
    % this finds out if the design was factorial or not and sets default
    % options accordingly
    if size(fieldnames(STATS.sample_results),1)==3;
        isfactorial=1;
        
        options = struct('FactorA', 1:size(STATS.sample_results.factor_A.contrasts,2), ...
            'FactorB', 1:size(STATS.sample_results.factor_B.contrasts,2), ...
            'FactorAB', 1:size(STATS.sample_results.factor_AxB.contrasts,2));
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.sample_results.factor_A.contrasts)
            disp('FactorB'); disp(STATS.sample_results.factor_B.contrasts)
            disp('FactorAB'); disp(STATS.sample_results.factor_AxB.contrasts)
        end
        
    elseif size(fieldnames(STATS.sample_results),1)==1;
        isfactorial=0;
        options = struct('FactorA', 1:size(STATS.sample_results.factor_A.contrasts,2));
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.sample_results.factor_A.contrasts)
        end
        
    end
    
    % get rid of 'all' in varargin option and convert ms entries into TFs
    remall=strcmp(varargin,'all');
    varargin(remall)=[];
    if any(strcmp(varargin,'timeplot'));
        timems=find(strcmp(varargin,'timeplot'));
        MStoTF_min=round((varargin{timems+1}(1)/1000-STATS.xmin)/(STATS.xmax-STATS.xmin) * (STATS.numpnts-1))+1;
        MStoTF_max=round((varargin{timems+1}(2)/1000-STATS.xmin)/(STATS.xmax-STATS.xmin) * (STATS.numpnts-1))+1;
        varargin{timems+1}=MStoTF_min:MStoTF_max;
        
    end
    
    
    nargs = length(varargin);
    if round(nargs/2)~=nargs/2
        error('need propertyName/propertyValue pairs for optional inputs')
    end
    
    % add other options
    options.timeplot=1:length(STATS.xtimes);
    
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

    % else check for name/val agreement and set defaults to empty
else
    nargs = length(varargin);
    if round(nargs/2)~=nargs/2
        error('need propertyName/propertyValue pairs for optional inputs')
    end
    
    % this finds out if the design was factorial or not and sets default
    % options accordingly
    if size(fieldnames(STATS.sample_results),1)==3;
        isfactorial=1;
        options = struct('FactorA', [],'FactorB', [], 'FactorAB', []);
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.sample_results.factor_A.contrasts)
            disp('FactorB'); disp(STATS.sample_results.factor_B.contrasts)
            disp('FactorAB'); disp(STATS.sample_results.factor_AxB.contrasts)
        end
        
    elseif size(fieldnames(STATS.sample_results),1)==1;
        isfactorial=0;
        options = struct('FactorA', []);
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.sample_results.factor_A.contrasts)
        end
        
        
    end
    
    
   % add other options
   options.timeplot=1:length(STATS.xtimes);
   if any(strcmp(varargin,'timeplot'));
       timems=find(strcmp(varargin,'timeplot'));
       MStoTF_min=round((varargin{timems+1}(1)/1000-STATS.xmin)/(STATS.xmax-STATS.xmin) * (STATS.numpnts-1))+1;
       MStoTF_max=round((varargin{timems+1}(2)/1000-STATS.xmin)/(STATS.xmax-STATS.xmin) * (STATS.numpnts-1))+1;
       varargin{timems+1}=MStoTF_min:MStoTF_max;
       
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
                if STATS.sample_results.factor_A.contrasts(j,options.FactorA(i))==1
                    plot1st(k,:)=STATS.condwaves(j,options.timeplot);
                    leg1stlist{k}=STATS.condnames{j};
                    k=k+1;
                elseif STATS.sample_results.factor_A.contrasts(j,options.FactorA(i))==-1
                    plot2nd(m,:)=STATS.condwaves(j,options.timeplot);
                    leg2ndlist{m}=STATS.condnames{j};
                    m=m+1;
                end
                
            end
            
            
            
            % reduce, if needed, the condition waveforms
            % should this be SUM or MEAN if pooling across levels?
            plot1st=mean(plot1st,1);
            plot2nd=mean(plot2nd,1);
            plotdiff=STATS.sample_results.factor_A.test_stat(options.FactorA(i),options.timeplot);
            
            % concatenate, if needed, the legend lables
            leg1st=strjoin_statslab(leg1stlist,'+');
            leg2nd=strjoin_statslab(leg2ndlist,'+');
            
            % begin plotting
            figure(i);
            subplot(2,1,1)
            h(1)=plot(STATS.xtimes(options.timeplot),plot1st,'r', 'LineWidth',3);
            hold on
            h(2)=plot(STATS.xtimes(options.timeplot),plot2nd,'b','LineWidth',3);
            set(gca,'FontSize',20)
            %title('CS vs CT')
            lh=legend(leg1st,leg2nd);
            % get rid of subscripts that occur when there are underscores
            set(lh,'Interpreter', 'none');
            axis tight
            grid on
            
            subplot(2,1,2)
            h(3)=plot(STATS.xtimes(options.timeplot),plotdiff,'k');
            CIup=STATS.sample_results.factor_A.CI{options.FactorA(i)}(2,options.timeplot);
            CIlow=STATS.sample_results.factor_A.CI{options.FactorA(i)}(1,options.timeplot);
            h(4)=jbfill(STATS.xtimes(options.timeplot),CIup, CIlow, [.5 .5 .5], [.5 .5 .5], 1, .6);
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
                    if STATS.sample_results.factor_A.contrasts(j,options.FactorA(i))==1
                        plot1st(k,:)=STATS.condwaves(j,options.timeplot);
                        leg1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.sample_results.factor_A.contrasts(j,options.FactorA(i))==-1
                        plot2nd(m,:)=STATS.condwaves(j,options.timeplot);
                        leg2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                end
                
                
                
                % reduce, if needed, the condition waveforms
                % should this be SUM or MEAN if pooling across levels?
                plot1st=mean(plot1st,1);
                plot2nd=mean(plot2nd,1);
                plotdiff=STATS.sample_results.factor_A.test_stat(options.FactorA(i),options.timeplot);
                
                % concatenate, if needed, the legend lables
                leg1st=strjoin_statslab(leg1stlist,'+');
                leg2nd=strjoin_statslab(leg2ndlist,'+');
                
                % begin plotting
                figure;
                subplot(2,1,1)
                h(1)=plot(STATS.xtimes(options.timeplot),plot1st,'r', 'LineWidth',3);
                set(gca,'ButtonDownFcn', @mouseclick_callback)
                set(get(gca,'Children'),'ButtonDownFcn', @mouseclick_callback)
                
                hold on
                h(2)=plot(STATS.xtimes(options.timeplot),plot2nd,'b','LineWidth',3);
                set(gca,'FontSize',20)
                %title('CS vs CT')
                lh=legend(leg1st,leg2nd);
                % get rid of subscripts that occur when there are underscores
                set(lh,'Interpreter', 'none');
                axis tight
                grid on
                
                subplot(2,1,2)
                h(3)=plot(STATS.xtimes(options.timeplot),plotdiff,'k');
                CIup=STATS.sample_results.factor_A.CI{options.FactorA(i)}(2,options.timeplot);%USE_THIS
                %CIup=STATS.sample_results.factor_A.CI{options.FactorA(i)}(2,513:922);
                
               
                
                CIlow=STATS.sample_results.factor_A.CI{options.FactorA(i)}(1,options.timeplot);
                %CIlow=STATS.sample_results.factor_A.CI{options.FactorA(i)}(1,513:922);
                
                h(4)=jbfill(STATS.xtimes(options.timeplot),CIup, CIlow, [.5 .5 .5], [.5 .5 .5], 1, .6);
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
                    if STATS.sample_results.factor_B.contrasts(j,options.FactorB(i))==1
                        plot1st(k,:)=STATS.condwaves(j,options.timeplot);
                        leg1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.sample_results.factor_B.contrasts(j,options.FactorB(i))==-1
                        plot2nd(m,:)=STATS.condwaves(j,options.timeplot);
                        leg2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                end
                
                
                
                % reduce, if needed, the condition waveforms
                % should this be SUM or MEAN if pooling across levels?
                plot1st=mean(plot1st,1);
                plot2nd=mean(plot2nd,1);
                plotdiff=STATS.sample_results.factor_B.test_stat(options.FactorB(i),options.timeplot);
                
                % concatenate, if needed, the legend lables
                leg1st=strjoin_statslab(leg1stlist,'+');
                leg2nd=strjoin_statslab(leg2ndlist,'+');
                
                % begin plotting
                figure;
                subplot(2,1,1)
                h(1)=plot(STATS.xtimes(options.timeplot),plot1st,'r', 'LineWidth',3);
                hold on
                h(2)=plot(STATS.xtimes(options.timeplot),plot2nd,'b','LineWidth',3);
                set(gca,'FontSize',20)
                %title('CS vs CT')
                lh=legend(leg1st,leg2nd);
                % get rid of subscripts that occur when there are underscores
                set(lh,'Interpreter', 'none');
                axis tight
                grid on
                
                subplot(2,1,2)
                h(3)=plot(STATS.xtimes(options.timeplot),plotdiff,'k');
                CIup=STATS.sample_results.factor_B.CI{options.FactorB(i)}(2,options.timeplot);
                CIlow=STATS.sample_results.factor_B.CI{options.FactorB(i)}(1,options.timeplot);
                h(4)=jbfill(STATS.xtimes(options.timeplot),CIup, CIlow, [.5 .5 .5], [.5 .5 .5], 1, .6);
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
                    if STATS.sample_results.factor_AxB.contrasts(j,options.FactorAB(i))==1
                        plot1st(k,:)=STATS.condwaves(j,options.timeplot);
                        leg1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.sample_results.factor_AxB.contrasts(j,options.FactorAB(i))==-1
                        plot2nd(m,:)=STATS.condwaves(j,options.timeplot);
                        leg2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                
                end

                % create interaction difference waveforms 
                plot1stdiff=plot1st(1,:)-plot2nd(1,:);
                plot2nddiff=plot2nd(2,:)-plot1st(2,:);
                plotdiff=STATS.sample_results.factor_AxB.test_stat(options.FactorAB(i),options.timeplot);
                
                % concatenate, if needed, the legend lables
                leg1st=strjoin_statslab({leg1stlist{1} leg2ndlist{1}},'-');
                leg2nd=strjoin_statslab({leg2ndlist{2} leg1stlist{2}} ,'-');
                
                % begin plotting
                figure;
                subplot(2,1,1)
                h(1)=plot(STATS.xtimes(options.timeplot),plot1stdiff,'r', 'LineWidth',3);
                hold on
                h(2)=plot(STATS.xtimes(options.timeplot),plot2nddiff,'b','LineWidth',3);
                set(gca,'FontSize',20)
                %title('CS vs CT')
                lh=legend(leg1st,leg2nd);
                % get rid of subscripts that occur when there are underscores
                set(lh,'Interpreter', 'none');
                axis tight
                grid on
                
                subplot(2,1,2)
                h(3)=plot(STATS.xtimes(options.timeplot),plotdiff,'k');
                CIup=STATS.sample_results.factor_AxB.CI{options.FactorAB(i)}(2,options.timeplot);
                CIlow=STATS.sample_results.factor_AxB.CI{options.FactorAB(i)}(1,options.timeplot);
                h(4)=jbfill(STATS.xtimes(options.timeplot),CIup, CIlow, [.5 .5 .5], [.5 .5 .5], 1, .6);
                axis tight
                grid on
                
                % clear variables before next iteration
                clear plot1st plot2nd plotdiff leg1st leg2nd leg2ndlist leg1stlist
                
                
                
            end
        end
   
end

disp('******* Saving STATS structure *******')
save(['STATS_',STATS.savestring,'.mat'],'STATS');

function mouseclick_callback(gcbo,eventdata)
      % the arguments are not important here, they are simply required for
      % a callback function. we don't even use them in the function,
      % but Matlab will provide them to our function, we we have to
      % include them.
      %
      % first we get the point that was clicked on
      cP = get(gca,'Currentpoint');
      x = cP(1,1);
      y = cP(1,2);
      % Now we find out which mouse button was clicked, and whether a
      % keyboard modifier was used, e.g. shift or ctrl
      %switch get(gcf,'SelectionType')
          %case 'normal' % Click left mouse button.
              s = sprintf('left: (%1.4g, %1.4g) level = %1.4g',x,y, x.*exp(-x.^2-y.^2));
          %case 'alt'    % Control - click left mouse button or click right mouse button.
          %   s = sprintf('right: (%1.4g, %1.4g level = %1.4g)',x,y, x.*exp(-x.^2-y.^2));
          %case 'extend' % Shift - click left mouse button or click both left and right mouse buttons.
          %    s = sprintf('2-click: (%1.4g, %1.4g level = %1.4g)',x,y, x.*exp(-x.^2-y.^2));
          %case 'open'   % Double-click any mouse button.
              %s = sprintf('double click: (%1.4g, %1.4g) level = %1.4g',x,y, x.*exp(-x.^2-y.^2));
      %end
      % get and set title handle
      thandle = get(gca,'Title');
      set(thandle,'String',s);
      % finally change the position of our red plus, and make it
      % visible.
      %set(cursor_handle,'Xdata',x,'Ydata',y,'visible','on')
end


end









