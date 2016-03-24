function [STATS]=pbsubjectfigTF(STATS,infodisplay,varargin)

try
    fields=fieldnames(STATS.sample_results);
catch
    error('you must run the same analysis at the group level, to plot single-subject time frequency results')
end

% freq band fields
sub_fields=fieldnames(STATS.subject_results);
band_fields=fieldnames(STATS.subject_results.(sub_fields{1}));

ersp=0;
itc=0;

% special case where using 'all' plots all possible contrasts
if any(strcmp(varargin,'all'));
    
    % this finds out if the design was factorial or not and sets default
    % options accordingly
    %     if size(fieldnames(STATS.subject_results.subject_1.(band_fields{1})),1)==3;
    %         isfactorial=1;
    %
    %         options = struct('FactorA', 1:size(STATS.subject_results.subject_1.(band_fields{1}).factor_A.contrasts,2), ...
    %             'FactorB', 1:size(STATS.subject_results.subject_1.(band_fields{1}).factor_B.contrasts,2), ...
    %             'FactorAB', 1:size(STATS.subject_results.subject_1.(band_fields{1}).factor_AxB.contrasts,2));
    %
    %         if infodisplay
    %             disp('Condition names'); disp(STATS.condnames)
    %             disp('FactorA'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_A.contrasts)
    %             disp('FactorB'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_B.contrasts)
    %             disp('FactorAB'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_AxB.contrasts)
    %         end
    %
    %     elseif size(fieldnames(STATS.subject_results.subject_1.(band_fields{1})),1)==1;
    %         isfactorial=0;
    %         options = struct('FactorA', 1:size(STATS.subject_results.subject_1.(band_fields{1}).factor_A.contrasts,2));
    %
    %         if infodisplay
    %             disp('Condition names'); disp(STATS.condnames)
    %             disp('FactorA'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_A.contrasts)
    %         end
    %
    %     end
    
    
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
    
    %     % this finds out if the design was factorial or not and sets default
    %     % options accordingly
    %     if size(fieldnames(STATS.subject_results.subject_1.(band_fields{1})),1)==3;
    %         isfactorial=1;
    %         options = struct('FactorA', [],'FactorB', [], 'FactorAB', []);
    %
    %         if infodisplay
    %             disp('Condition names'); disp(STATS.condnames)
    %             disp('FactorA'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_A.contrasts)
    %             disp('FactorB'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_B.contrasts)
    %             disp('FactorAB'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_AxB.contrasts)
    %         end
    %
    %     elseif size(fieldnames(STATS.subject_results.subject_1.(band_fields{1})),1)==1;
    %         isfactorial=0;
    %         options = struct('FactorA', []);
    %
    %         if infodisplay
    %             disp('Condition names'); disp(STATS.condnames)
    %             disp('FactorA'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_A.contrasts)
    %         end
    %
    %
    %     end
    %
    
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

% preallocate
pval_gather=zeros(length(band_fields), length(STATS.TF_times));
pval_tmp=zeros(length(band_fields), length(STATS.TF_times));

% switch cases for factorial vs. single-factor
switch STATS.design
    case 'w'
        
        numfigs=size(options.FactorA,2);
        
        for i=1:numfigs
            k=1;
            m=1;
            
            % sub loop
            for s=1:length(sub_fields);
                
                % stats
                for q=1:STATS.freqbins;
                    pval_tmp(q,:)=(STATS.subject_results.(sub_fields{s}).(band_fields{q}).factor_A.pval(options.FactorA(i),:)<=.05);
                end
                
                pval_gather=pval_gather+pval_tmp;
            end
            
            % pval proportions
            pval_prop=pval_gather./length(sub_fields);
            
            
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
            figure(i);
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
            h(4)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pval_prop(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
            colormap(flipud(hot)); caxis([0 1]); hb=cbfreeze(colorbar);  freezeColors; %set(hb,'YTick',[0 1]);
            
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
        
        
    case 'ww'
        
        if ~isempty(options.FactorA);
            numfigs=size(options.FactorA,2);
            
            for i=1:numfigs
                k=1;
                m=1;
                
                for s=1:length(sub_fields);
                    
                    % stats
                    for q=1:STATS.freqbins;
                        pval_tmp(q,:)=(STATS.subject_results.(sub_fields{s}).(band_fields{q}).factor_A.pval(options.FactorA(i),:)<=.05);
                    end
                    
                    pval_gather=pval_gather+pval_tmp;
                end
                
                % pval proportions
                pval_prop=pval_gather./length(sub_fields);
                
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
                h(4)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pval_prop(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                colormap(flipud(hot)); caxis([0 1]); hb=cbfreeze(colorbar);  freezeColors; %set(hb,'YTick',[0 1]);
                
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
                
                
                for s=1:length(sub_fields);
                    
                    % stats
                    for q=1:STATS.freqbins;
                        pval_tmp(q,:)=(STATS.subject_results.(sub_fields{s}).(band_fields{q}).factor_B.pval(options.FactorB(i),:)<=.05);
                    end
                    
                    pval_gather=pval_gather+pval_tmp;
                end
                
                % pval proportions
                pval_prop=pval_gather./length(sub_fields);
                
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
                h(4)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pval_prop(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                colormap(flipud(hot)); caxis([0 1]); hb=cbfreeze(colorbar);  freezeColors; %set(hb,'YTick',[0 1]);
                
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
                
                for s=1:length(sub_fields);
                    
                    % stats
                    for q=1:STATS.freqbins;
                        pval_tmp(q,:)=(STATS.subject_results.(sub_fields{s}).(band_fields{q}).factor_AxB.pval(options.FactorAB(i),:)<=.05);
                    end
                    
                    pval_gather=pval_gather+pval_tmp;
                end
                
                % pval proportions
                pval_prop=pval_gather./length(sub_fields);
                
                % group level interaction difference data
                for q=1:STATS.freqbins;
                    plotdiff(q,:)=STATS.sample_results.(fields{q}).factor_AxB.test_stat(options.FactorAB(i),:);
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
                
                set(allchild(gca),'buttondownfcn',{@mouseclick_callback, STATS, ['(',leg1st,')','-','(',leg2nd,')'], options.timeplot, options.FactorAB(i), 'AxB'});
                
                hsub(2)=subplot(2,1,2);
                h(2)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pval_prop(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                colormap(flipud(hot)); caxis([0 1]); hb=cbfreeze(colorbar);  freezeColors; %set(hb,'YTick',[0 1]);
                
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
        
        
        % work from here
    case 'bw' % series of 1-ways
        
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



    function mouseclick_callback(gcbo, eventdata, STATS, titl, timeplot, oncon, onfact)
        
        diffcol=[0 0 0];
        CIcol=[.5 .5 .5];
        
        flds=fieldnames(STATS.sample_results);
        subflds=fieldnames(STATS.subject_results);
        %get the point that was clicked on
        cP = get(gca,'Currentpoint');
        % ms_plot = cP(1,1);
        freqclick = cP(1,2);
        [V I]=min(abs(STATS.TF_freqs-freqclick));
        fact=['factor_', onfact];
        
        figure;
        for qq=1:length(subflds);
            
            % extract CI
            CI_lower(1,:)=STATS.subject_results.(subflds{qq}).(flds{I}).(fact).CI{oncon,1}(1,timeplot);
            CI_upper(1,:)=STATS.subject_results.(subflds{qq}).(flds{I}).(fact).CI{oncon,1}(2,timeplot);
            
            % extract difference wave
            pdiff=STATS.subject_results.(subflds{qq}).(flds{I}).(fact).test_stat(oncon,timeplot);
            
            % begin plotting
            
            subplot(length(subflds),1,qq);
            % plot zeroline
            plot(STATS.TF_times(timeplot),zeros(1,length(STATS.TF_times(timeplot))),'r','LineWidth',1);
            hold on
            
            % plot CI
            jbfill(STATS.TF_times(timeplot),CI_upper,CI_lower,CIcol,CIcol,1,1);
            hold on
            
            % plot diff wave
            plot(STATS.TF_times(timeplot),pdiff,'Color',diffcol);
            axis tight
            grid on
        end
        
        htit = axes('visible','off');
        title([titl, ' at ', num2str(STATS.TF_freqs(I)), ' Hz'], 'parent', htit, 'visible', 'on');
    end

end









