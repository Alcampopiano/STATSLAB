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
    options.savesvg='no';
    
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
    options.savesvg='no';
    
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%   
            if strcmp(STATS.sample_results.(fields{1}).factor_A.FWE,'benhoch');
                
                % stats
                for q=1:STATS.freqbins;
                    pvals(q,:)=(STATS.sample_results.(fields{q}).factor_A.sig_pvals(options.FactorA(i),:)==0);
                end
                
                
            else 
                % stats
                for q=1:STATS.freqbins;
                    pvals(q,:)=(STATS.sample_results.(fields{q}).factor_A.pval(options.FactorA(i),:)>STATS.alpha);
                    %pvals(q,:)=(STATS.sample_results.(fields{q}).factor_A.pval(options.FactorA(i),:)>.05);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
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
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%5
            % begin plotting
            figure;
            hsub(1)=subplot(4,1,1);
            h(1)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plot1st(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
            
            %svg
            if strcmp(options.savesvg,'yes'); buildsvg(h(1),conds_ax,STATS,'ersp_group_subplot_1.svg','jet',options.timeplot); end
            
            if ersp==1; colormap(jet); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
            if itc==1; colormap(flipud(hot)); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
            
            hsub(2)=subplot(4,1,2);
            h(2)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plot2nd(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
            
            %svg
            if strcmp(options.savesvg,'yes'); buildsvg(h(2),conds_ax,STATS,'ersp_group_subplot_2.svg','jet',options.timeplot); end
            
            if ersp==1; colormap(jet); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
            if itc==1; colormap(flipud(hot)); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
            
            hsub(3)=subplot(4,1,3);
            h(3)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plotdiff(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
            
            %svg
            if strcmp(options.savesvg,'yes'); buildsvg(h(3),conds_axdiff,STATS,'ersp_group_subplot_3.svg','jet',options.timeplot); end
            
            colormap(jet); caxis(conds_axdiff); cbfreeze(colorbar); freezeColors;
            
            set(allchild(gca),'buttondownfcn',{@mouseclick_callback, STATS, [leg1st,'-',leg2nd], options.timeplot, options.FactorA(i), 'A'});
            
            hsub(4)=subplot(4,1,4);
            h(4)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pvals(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
            
            %svg
            if strcmp(options.savesvg,'yes'); buildsvg(h(4),[0 1],STATS,'ersp_group_subplot_4.svg','bone',options.timeplot); end
            
            colormap(hsub(4),bone(2)); caxis([0 1]); hb=cbfreeze(colorbar); set(hb,'YTick',[0 1]); %freezeColors;
            
            % add legend and font
            set(gca,'FontSize',10)
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
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%   
            if strcmp(STATS.sample_results.(fields{1}).factor_A.FWE,'benhoch');
                
                % stats
                for q=1:STATS.freqbins;
                    pvals(q,:)=(STATS.sample_results.(fields{q}).factor_A.sig_pvals(options.FactorA(i),:)==0);
                end
                
                
            else 
                % stats
                for q=1:STATS.freqbins;
                    pvals(q,:)=(STATS.sample_results.(fields{q}).factor_A.pval(options.FactorA(i),:)>STATS.alpha);
                    %pvals(q,:)=(STATS.sample_results.(fields{q}).factor_A.pval(options.FactorA(i),:)>.05);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
             
                
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
                
                set(allchild(gca),'buttondownfcn',{@mouseclick_callback, STATS, [leg1st,'-',leg2nd], options.timeplot, options.FactorA(i), 'A'});
                
                hsub(4)=subplot(4,1,4);
                h(4)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pvals(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                colormap(hsub(4),bone(2)); caxis([0 1]); hb=cbfreeze(colorbar); set(hb,'YTick',[0 1]); %freezeColors;
                
                % add legend and font
                set(gca,'FontSize',10);
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
                
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%   
            if strcmp(STATS.sample_results.(fields{1}).factor_B.FWE,'benhoch');
                
                % stats
                for q=1:STATS.freqbins;
                    pvals(q,:)=(STATS.sample_results.(fields{q}).factor_B.sig_pvals(options.FactorB(i),:)==0);
                end
                
                
            else 
                % stats
                for q=1:STATS.freqbins;
                    pvals(q,:)=(STATS.sample_results.(fields{q}).factor_B.pval(options.FactorB(i),:)>STATS.alpha);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
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
                
                set(allchild(gca),'buttondownfcn',{@mouseclick_callback, STATS, [leg1st,'-',leg2nd], options.timeplot, options.FactorB(i), 'B'});
                
                hsub(4)=subplot(4,1,4);
                h(4)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pvals(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                colormap(hsub(4),bone(2)); caxis([0 1]); hb=cbfreeze(colorbar); set(hb,'YTick',[0 1]); %freezeColors;
                
                % add legend and font
                set(gca,'FontSize',10)
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
                
             %%%%%%%%%%%%%%%%%%%%%%%%%%   
            if strcmp(STATS.sample_results.(fields{1}).factor_AxB.FWE,'benhoch');
                
                % stats
                for q=1:STATS.freqbins;
                    pvals(q,:)=(STATS.sample_results.(fields{q}).factor_AxB.sig_pvals(options.FactorAB(i),:)==0);
                    plotdiff(q,:)=STATS.sample_results.(fields{q}).factor_AxB.test_stat(options.FactorAB(i),:);
                end
                
                
            else 
                % stats
                for q=1:STATS.freqbins;
                    pvals(q,:)=(STATS.sample_results.(fields{q}).factor_AxB.pval(options.FactorAB(i),:)>STATS.alpha);
                    plotdiff(q,:)=STATS.sample_results.(fields{q}).factor_AxB.test_stat(options.FactorAB(i),:);

                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
                
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
      
                set(allchild(gca),'buttondownfcn',{@mouseclick_callback, STATS, ['(',leg1st,')','-','(',leg2nd,')'], options.timeplot, options.FactorAB(i), 'AxB'});
                
                hsub(2)=subplot(2,1,2);
                h(2)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pvals(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                colormap(hsub(2),bone(2)); caxis([0 1]); hb=cbfreeze(colorbar); set(hb,'YTick',[0 1]); %freezeColors;
                
                % add legend and font
                set(gca,'FontSize',10)
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

    function mouseclick_callback(gcbo, eventdata, STATS, titl, timeplot, oncon, onfact)
        
        diffcol=[0 0 0];
        CIcol=[.5 .5 .5];
        
        flds=fieldnames(STATS.sample_results);
        %get the point that was clicked on
        cP = get(gca,'Currentpoint');
        % ms_plot = cP(1,1);
        freqclick = cP(1,2);
        [V I]=min(abs(STATS.TF_freqs-freqclick));
        fact=['factor_', onfact];
        
        % extract CI
        CI_low(1,:)=STATS.sample_results.(flds{I}).(fact).CI{oncon,1}(1,timeplot);
        CI_up(1,:)=STATS.sample_results.(flds{I}).(fact).CI{oncon,1}(2,timeplot);
        
        % extract difference wave
        pdiff=STATS.sample_results.(flds{I}).(fact).test_stat(oncon,timeplot);
        
        % begin plotting
        figure;
        
        % plot zeroline
        plot(STATS.TF_times(timeplot),zeros(1,length(STATS.TF_times(timeplot))),'r','LineWidth',1);
        hold on
        
        % plot CI
        jbfill(STATS.TF_times(timeplot),CI_up,CI_low,CIcol,CIcol,1,1);
        hold on
        
        % plot diff wave
        plot(STATS.TF_times(timeplot),pdiff,'Color',diffcol);
        axis tight
        grid on
        
        tit=title([titl, ' at ', num2str(STATS.TF_freqs(I)), ' Hz']);
        set(tit,'Interpreter', 'none');
    end


end









