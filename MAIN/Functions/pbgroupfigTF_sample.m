function [STATS]=pbgroupfigTF_sample(STATS,infodisplay,varargin)

% freq band fields
fields=fieldnames(STATS.sample_results);
ersp=0;
itc=0;

% special case where using 'all' plots all possible contrasts
if any(strcmp(varargin,'all'));
    
    % this finds out if the design was factorial or not and sets default
    % options accordingly
    if size(fieldnames(STATS.sample_results.(fields{1})),1)==3;
        isfactorial=1;
        
        options = struct('FactorA', 1:size(STATS.sample_results.(fields{1}).factor_A.contrasts,2), ...
            'FactorB', 1:size(STATS.sample_results.(fields{1}).factor_B.contrasts,2), ...
            'FactorAB', 1:size(STATS.sample_results.(fields{1}).factor_AxB.contrasts,2));
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.sample_results.(fields{1}).factor_A.contrasts)
            disp('FactorB'); disp(STATS.sample_results.(fields{1}).factor_B.contrasts)
            disp('FactorAB'); disp(STATS.sample_results.(fields{1}).factor_AxB.contrasts)
        end
        
    elseif size(fieldnames(STATS.sample_results.(fields{1})),1)==1;
        isfactorial=0;
        options = struct('FactorA', 1:size(STATS.sample_results.(fields{1}).factor_A.contrasts,2));
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.sample_results.(fields{1}).factor_A.contrasts)
        end
        
    end
    
    
    % get rid of 'all' in varargin option and convert ms entries into TFs
    remall=strcmp(varargin,'all');
    varargin(remall)=[];
    if any(strcmp(varargin,'timeplot'));
        timems=find(strcmp(varargin,'timeplot'));
        tri = delaunayn(STATS.TF_times');
        ind=dsearchn(STATS.TF_times',tri,[varargin{timems+1}(1) varargin{timems+1}(2)]');
        %         MStoTF_min=STATS.TF_times(ind(1));
        %         MStoTF_max=STATS.TF_times(ind(2));
        varargin{timems+1}=ind(1):(ind(2));
    end
    
    nargs = length(varargin);
    if round(nargs/2)~=nargs/2
        error('need propertyName/propertyValue pairs for optional inputs')
    end
    
    % add other options
    options.timeplot=1:length(STATS.TF_times);
    
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
    if size(fieldnames(STATS.sample_results.(fields{1})),1)==3;
        isfactorial=1;
        options = struct('FactorA', [],'FactorB', [], 'FactorAB', []);
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.sample_results.(fields{1}).factor_A.contrasts)
            disp('FactorB'); disp(STATS.sample_results.(fields{1}).factor_B.contrasts)
            disp('FactorAB'); disp(STATS.sample_results.(fields{1}).factor_AxB.contrasts)
        end
        
    elseif size(fieldnames(STATS.sample_results.(fields{1})),1)==1;
        isfactorial=0;
        options = struct('FactorA', []);
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.sample_results.(fields{1}).factor_A.contrasts)
        end
        
        
    end
    
    
    % add other options
    options.timeplot=1:length(STATS.TF_times);
    
    if any(strcmp(varargin,'timeplot'));
        timems=find(strcmp(varargin,'timeplot'));
        tri = delaunayn(STATS.TF_times');
        ind=dsearchn(STATS.TF_times',tri,[varargin{timems+1}(1) varargin{timems+1}(2)]');
        %         MStoTF_min=STATS.TF_times(ind(1));
        %         MStoTF_max=STATS.TF_times(ind(2));
        varargin{timems+1}=ind(1):(ind(2));
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

            % stats
            for q=1:STATS.freqbins;
                pvals(q,:)=(STATS.sample_results.(fields{q}).factor_A.pval(options.FactorA(i),:)>.05);
            end
            
            for j=1:length(STATS.condnames);
                
                % get the condition waveforms
                if STATS.sample_results.(fields{1}).factor_A.contrasts(j,options.FactorA(i))==1
                    
                    c=load(['group_TFwaves_',STATS.savestring,'_',STATS.condnames{j}]);
                    plot1st(:,:,k)=c.condition;
                    leg1stlist{k}=STATS.condnames{j};
                    k=k+1;
                elseif STATS.sample_results.(fields{1}).factor_A.contrasts(j,options.FactorA(i))==-1
                    
                    c=load(['group_TFwaves_',STATS.savestring,'_',STATS.condnames{j}]);
                    plot2nd(:,:,m)=c.condition;
                    leg2ndlist{m}=STATS.condnames{j};
                    m=m+1;
                end
                
            end
            
            
            % reduce, if needed, the condition waveforms
            % should this be SUM or MEAN if pooling across levels?
            plot1st=mean(plot1st,3);
            plot2nd=mean(plot2nd,3);
            plotdiff=plot1st-plot2nd;
            %plotdiff=STATS.sample_results.factor_A.test_stat(options.FactorA(i),options.timeplot);
            
            % concatenate, if needed, the legend lables
            leg1st=strjoin_statslab(leg1stlist,'+');
            leg2nd=strjoin_statslab(leg2ndlist,'+');
            
            % set caxis for handles based on data type, and color maps
            if strcmp(STATS.measure,'ersp');
                ersp=1;
                conds_ax=max([max(max(abs(plot1st))) max(max(abs(plot2nd)))]);
                conds_ax=[-conds_ax conds_ax];
                
            elseif strcmp(STATS.measure,'itc');
                itc=1;
                conds_ax=max([max(max(abs(plot1st))) max(max(abs(plot2nd)))]);
                conds_axmin=min([min(min(abs(plot1st))) min(min(abs(plot2nd)))]);
                conds_ax=[conds_axmin conds_ax];
            end
            
            % difference limits
            conds_axdiff=max(max(abs(plotdiff)));
            conds_axdiff=[-conds_axdiff conds_axdiff];
            
            % begin plotting
            figure;
            hsub(1)=subplot(4,1,1);
            h(1)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plot1st(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
            if ersp==1; colormap(jet); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
            if itc==1; colormap(flipud(hot)); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
            
            hsub(2)=subplot(4,1,2);
            h(2)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plot2nd(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
            if ersp==1; colormap(jet); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
            if itc==1; colormap(flipud(hot)); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
            
            hsub(3)=subplot(4,1,3);
            h(3)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plotdiff(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
            colormap(jet); caxis(conds_axdiff); cbfreeze(colorbar); freezeColors;
            
            set(gca,'ButtonDownFcn', {@mouseclick_callback,STATS});
            set(allchild(gca),'buttondownfcn',{@mouseclick_callback,STATS, [leg1st,'-',leg2nd], options.timeplot});
            
            hsub(4)=subplot(4,1,4);
            h(4)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pvals(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
            colormap(bone(2)); caxis([0 1]); hb=cbfreeze(colorbar); set(hb,'YTick',[0 1]); freezeColors; 
            
            % add legend and font
            set(gca,'FontSize',20)
            h=axes('visible','off');
            lh=title([leg1st,'-',leg2nd],'parent',h,'visible','on','Position',[.5 1.05 0]);
            
            % get rid of subscripts that occur when there are underscores
            set(lh,'Interpreter', 'none');
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
                
                % stats
                for q=1:STATS.freqbins;
                    % diffsA1(i,:)=STATS.sample_results.(bands{i}).factor_A.test_stat(1,:);
                    pvals(q,:)=(STATS.sample_results.(fields{q}).factor_A.pval(options.FactorA(i),:)>.05);
                    % pvalAB(i,:)=(STATS.sample_results.(bands{i}).factor_AxB.pval(1,:)>.05);
                end
                
                for j=1:length(STATS.condnames);
                    
                    % get the condition waveforms
                    if STATS.sample_results.(fields{1}).factor_A.contrasts(j,options.FactorA(i))==1
                        
                        c=load(['group_TFwaves_',STATS.savestring,'_',STATS.condnames{j}]);
                        plot1st(:,:,k)=c.condition;
                        leg1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.sample_results.(fields{1}).factor_A.contrasts(j,options.FactorA(i))==-1
                        
                        c=load(['group_TFwaves_',STATS.savestring,'_',STATS.condnames{j}]);
                        plot2nd(:,:,m)=c.condition;
                        leg2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                end
                
                
                % reduce, if needed, the condition waveforms
                % should this be SUM or MEAN if pooling across levels?
                plot1st=mean(plot1st,3);
                plot2nd=mean(plot2nd,3);
                plotdiff=plot1st-plot2nd;
                %plotdiff=STATS.sample_results.factor_A.test_stat(options.FactorA(i),options.timeplot);
                
                % concatenate, if needed, the legend lables
                leg1st=strjoin_statslab(leg1stlist,'+');
                leg2nd=strjoin_statslab(leg2ndlist,'+');
                
                % set caxis for handles based on data type, and color maps
                if strcmp(STATS.measure,'ersp');
                    ersp=1;
                    conds_ax=max([max(max(abs(plot1st))) max(max(abs(plot2nd)))]);
                    conds_ax=[-conds_ax conds_ax];
                    
                elseif strcmp(STATS.measure,'itc');
                    itc=1;
                    conds_ax=max([max(max(abs(plot1st))) max(max(abs(plot2nd)))]);
                    conds_axmin=min([min(min(abs(plot1st))) min(min(abs(plot2nd)))]);
                    conds_ax=[conds_axmin conds_ax];
                end
                
                % difference limits
                conds_axdiff=max(max(abs(plotdiff)));
                conds_axdiff=[-conds_axdiff conds_axdiff];
                
                % begin plotting
                figure;
                hsub(1)=subplot(4,1,1);
                h(1)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plot1st(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                if ersp==1; colormap(jet); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
                if itc==1; colormap(flipud(hot)); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
                
                hsub(2)=subplot(4,1,2);
                h(2)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plot2nd(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                if ersp==1; colormap(jet); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
                if itc==1; colormap(flipud(hot)); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
                
                hsub(3)=subplot(4,1,3);
                h(3)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plotdiff(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                colormap(jet); caxis(conds_axdiff); cbfreeze(colorbar); freezeColors;
                
                hsub(4)=subplot(4,1,4);
                h(4)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pvals(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                colormap(bone(2)); caxis([0 1]); hb=cbfreeze(colorbar); set(hb,'YTick',[0 1]); freezeColors;
                
                % add legend and font
                set(gca,'FontSize',20);
                h=axes('visible','off');
                lh=title([leg1st,'-',leg2nd],'parent',h,'visible','on','Position',[.5 1.05 0]);
                
                % get rid of subscripts that occur when there are underscores
                set(lh,'Interpreter', 'none');
                axis tight
                grid on
                
                % clear variables before next iteration
                clear plot1st plot2nd plotdiff leg1st leg2nd leg2ndlist leg1stlist lh
                
              
            end
        end
        
        if ~isempty(options.FactorB);
            numfigs=size(options.FactorB,2);
 
            for i=1:numfigs
                k=1;
                m=1;
                
                
                % stats
                for q=1:STATS.freqbins;
                    pvals(q,:)=(STATS.sample_results.(fields{q}).factor_B.pval(options.FactorB(i),:)>.05);
                end
                
                for j=1:length(STATS.condnames);
                    
                    % get the condition waveforms
                    if STATS.sample_results.(fields{1}).factor_B.contrasts(j,options.FactorB(i))==1
                        
                        c=load(['group_TFwaves_',STATS.savestring,'_',STATS.condnames{j}]);
                        plot1st(:,:,k)=c.condition;
                        leg1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.sample_results.(fields{1}).factor_B.contrasts(j,options.FactorB(i))==-1
                        
                        c=load(['group_TFwaves_',STATS.savestring,'_',STATS.condnames{j}]);
                        plot2nd(:,:,m)=c.condition;
                        leg2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end
                    
                end
                
                
                % reduce, if needed, the condition waveforms
                % should this be SUM or MEAN if pooling across levels?
                plot1st=mean(plot1st,3);
                plot2nd=mean(plot2nd,3);
                plotdiff=plot1st-plot2nd;
                %plotdiff=STATS.sample_results.factor_A.test_stat(options.FactorA(i),options.timeplot);
                
                % concatenate, if needed, the legend lables
                leg1st=strjoin_statslab(leg1stlist,'+');
                leg2nd=strjoin_statslab(leg2ndlist,'+');
                
                % set caxis for handles based on data type, and color maps
                if strcmp(STATS.measure,'ersp');
                    ersp=1;
                    conds_ax=max([max(max(abs(plot1st))) max(max(abs(plot2nd)))]);
                    conds_ax=[-conds_ax conds_ax];
                    
                elseif strcmp(STATS.measure,'itc');
                    itc=1;
                    conds_ax=max([max(max(abs(plot1st))) max(max(abs(plot2nd)))]);
                    conds_axmin=min([min(min(abs(plot1st))) min(min(abs(plot2nd)))]);
                    conds_ax=[conds_axmin conds_ax];
                end
                
                % difference limits
                conds_axdiff=max(max(abs(plotdiff)));
                conds_axdiff=[-conds_axdiff conds_axdiff];
                
                % begin plotting
                figure;
                hsub(1)=subplot(4,1,1);
                h(1)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plot1st(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                if ersp==1; colormap(jet); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
                if itc==1; colormap(flipud(hot)); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
                
                hsub(2)=subplot(4,1,2);
                h(2)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plot2nd(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                if ersp==1; colormap(jet); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
                if itc==1; colormap(flipud(hot)); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
                
                hsub(3)=subplot(4,1,3);
                h(3)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plotdiff(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                colormap(jet); caxis(conds_axdiff); cbfreeze(colorbar); freezeColors;
                
                hsub(4)=subplot(4,1,4);
                h(4)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pvals(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                colormap(bone(2)); caxis([0 1]); hb=cbfreeze(colorbar); set(hb,'YTick',[0 1]); freezeColors;
                
                % add legend and font
                set(gca,'FontSize',20)
                h=axes('visible','off');
                lh=title([leg1st,'-',leg2nd],'parent',h,'visible','on','Position',[.5 1.05 0]);
                
                % get rid of subscripts that occur when there are underscores
                set(lh,'Interpreter', 'none');
                axis tight
                grid on
                
                % clear variables before next iteration
                clear plot1st plot2nd plotdiff leg1st leg2nd leg2ndlist leg1stlist lh
                
            end
        end
        
        if ~isempty(options.FactorAB);
            numfigs=size(options.FactorAB,2);
            
            for i=1:numfigs
                k=1;
                m=1;
                
                % stats
                for q=1:STATS.freqbins;
                    plotdiff(q,:)=STATS.sample_results.(fields{q}).factor_AxB.test_stat(options.FactorAB(i),:);
                    pvals(q,:)=(STATS.sample_results.(fields{q}).factor_AxB.pval(options.FactorAB(i),:)>.05);
                end
                
                for j=1:length(STATS.condnames);
                    
                    
                    % get the condition waveforms
                    if STATS.sample_results.(fields{1}).factor_AxB.contrasts(j,options.FactorAB(i))==1
                        
                        c=load(['group_TFwaves_',STATS.savestring,'_',STATS.condnames{j}]);
                        plot1st(:,:,k)=c.condition;
                        leg1stlist{k}=STATS.condnames{j};
                        k=k+1;
                    elseif STATS.sample_results.(fields{1}).factor_AxB.contrasts(j,options.FactorAB(i))==-1
                        
                        c=load(['group_TFwaves_',STATS.savestring,'_',STATS.condnames{j}]);
                        plot2nd(:,:,m)=c.condition;
                        leg2ndlist{m}=STATS.condnames{j};
                        m=m+1;
                    end 
                    
                end
                
                plot1stdiff=plot1st(:,:,1)-plot2nd(:,:,1);
                plot2nddiff=plot2nd(:,:,2)-plot1st(:,:,2);
                %plotdiff=STATS.sample_results.factor_AxB.test_stat(options.FactorAB(i),options.timeplot);
                
                % create interaction difference waveforms
                %plot1stdiff=plot1st(1,:)-plot2nd(1,:);
                %plot2nddiff=plot2nd(2,:)-plot1st(2,:);
                %plotdiff=STATS.sample_results.factor_AxB.test_stat(options.FactorAB(i),options.timeplot);
                
                % concatenate, if needed, the legend lables
                leg1st=strjoin_statslab({leg1stlist{1} leg2ndlist{1}},'-');
                leg2nd=strjoin_statslab({leg2ndlist{2} leg1stlist{2}} ,'-');
                
                % set caxis for handles based on data type, and color maps
                if strcmp(STATS.measure,'ersp');
                    %ersp=1;
                    conds_ax=max([max(max(abs(plot1st))) max(max(abs(plot2nd)))]);
                    conds_ax=[-conds_ax conds_ax];
                    
                elseif strcmp(STATS.measure,'itc');
                    %itc=1;
                    conds_ax=max([max(max(abs(plot1st))) max(max(abs(plot2nd)))]);
                    conds_axmin=min([min(min(abs(plot1st))) min(min(abs(plot2nd)))]);
                    conds_ax=[conds_axmin conds_ax];
                end
                
                % difference limits
                conds_axdiff=max(max(abs(plotdiff)));
                conds_axdiff=[-conds_axdiff conds_axdiff];
                
                figure;
                hsub(1)=subplot(2,1,1);
                h(1)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plotdiff(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                colormap(jet); caxis(conds_axdiff); cbfreeze(colorbar); freezeColors;
                
                hsub(2)=subplot(2,1,2);
                h(2)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pvals(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                colormap(bone(2)); caxis([0 1]); hb=cbfreeze(colorbar); set(hb,'YTick',[0 1]); freezeColors;
                
                % add legend and font
                set(gca,'FontSize',20)
                h=axes('visible','off');
                lh=title(['(',leg1st,')','-','(',leg2nd,')'],'parent',h,'visible','on','Position',[.5 1.05 0]);
                % get rid of subscripts that occur when there are underscores
                set(lh,'Interpreter', 'none');
                axis tight
                grid on
                
                % clear variables before next iteration
                clear plot1st plot2nd plotdiff leg1st leg2nd leg2ndlist leg1stlist lh
               
            end
        end       
end

disp('******* Saving STATS structure *******')
save(['STATS_',STATS.savestring,'.mat'],'STATS');

function mouseclick_callback(gcbo,eventdata,STATS, titl, timeplot)
        
        try

            %SIZEBOX=250; % some arbitrary size of some box
            
%             if strcmp(STATS.design, 'bw')
%                 
%                 atcond=(atjlvl-1)*STATS.levels(2)+1; % to determine which cell to start from
%                 col=atcond+STATS.levels(2)-1; % the last cell to plot
%                 col_size=STATS.levels(2);
%                 
%                 % get dimensions.
%                 rowcols(2) = ceil(sqrt(col_size)); % EEGpage is number of subjects/topos
%                 rowcols(1) = ceil(col_size/rowcols(2));
%                 
%             else
%                 
%                 [row col]=size(STATS.grouptopofiles);
%                 atcond=1;
%                 
%                 % get dimensions.
%                 rowcols(2) = ceil(sqrt(col)); % EEGpage is number of subjects/topos
%                 rowcols(1) = ceil(col/rowcols(2));
%                 
%             end

            %flds=fieldnames(STATS.sample_results);
            %get the point that was clicked on
            cP = get(gca,'Currentpoint');
            % ms_plot = cP(1,1);
            freqclick = cP(1,2);
            [V I]=min(abs(STATS.TF_freqs-freqclick));
            
          
            
            
%             s=1;
%             for r=atcond:col % loop for each subject?
%                 
%                 if r==atcond;
%                     % build eventual destination figure
%                     curfig = figure('paperpositionmode', 'auto', 'visible', 'off');
%                     pos = get(curfig,'Position');
%                     posx = max(0, pos(1)+(pos(3)-SIZEBOX*rowcols(2))/2);
%                     posy = pos(2)+pos(4)-SIZEBOX*rowcols(1);
%                     set(curfig,'Position', [posx posy  SIZEBOX*rowcols(2)  SIZEBOX*rowcols(1)]);
%                     
%                 end
%                 curax = subplot( rowcols(1), rowcols(2), mod(r-1, rowcols(1)*rowcols(2))+1);
%                 set(curax, 'visible', 'off')
%                 
%                 % loading subject
%                 data=load(STATS.subtopofiles{r}{sub});
%                 disp(STATS.subtopofiles{r}{sub});
%                 
%                 MStoTF=round((ms_plot/1000-data.EEG.xmin)/(data.EEG.xmax-data.EEG.xmin) * (data.EEG.pnts-1))+1;
%                 maplim=max(max(abs(data.EEG.data(:,MStoTF))));
%                 
%                 if isempty(maplim);
%                     pop_topoplot(data.EEG, 1, ms_plot, [], 0,'shading','interp','colorbar','off');
%                     %set(gcf, 'visible','off');
%                     htopo(s)=gca;
%                     ctopo(s,:)=caxis;
%                 else
%                     pop_topoplot(data.EEG, 1, ms_plot, [], 0,'shading','interp','colorbar','off','maplimits', [-maplim maplim]);
%                     %set(gcf, 'visible','off');
%                     htopo(s)=gca;
%                     ctopo(s,:)=caxis;
%                 end
%                 
%                 oh=findobj(curax); % find and get rid of EEGLABs subplot titles
%                 alltext=findall(oh,'Type','text');
%                 delete(alltext);
%                 text(.5,-.1,num2str(STATS.condnames{r}),'Units','normalized','Interpreter', 'none'); % add subject numbers to bottom centre of subplots
%                 
%                 if isempty(maplim)
%                     colorbar;
%                 end
% 
%                 if ~isempty(maplim)
%                     if r==col % last subject
%                         
%                         maplim=max(max(abs(ctopo)));
%                         
%                         % set limits
%                         for qq=1:length(htopo);
%                             caxis(htopo(qq), [-maplim maplim])
%                         end
%                         
%                         
%                         hax=axes('visible', 'off');
%                         set(hax, 'Units', 'Normalized', 'Position', [.88, 0.25, .025, .5]);
%                         colorbar('FontSize',15);
%                         caxis([-maplim maplim]);
%        
% %                         curax_pos=get(curax,'position');
% %                         colorbar('location','eastoutside');
% %                         set(curax,'position',curax_pos);
%                     end
%                 end
%                 s=s+1;
%                 
%             end % end of r loop
            htit = axes('visible','off');
            %title(['Freq', num2str(freqclick)],'parent',htit,'visible','on');
            
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









