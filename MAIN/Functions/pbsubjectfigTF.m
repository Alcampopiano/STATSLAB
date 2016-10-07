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
        
        options = struct('FactorA', 1:size(STATS.subject_results.subject_1.(band_fields{1}).factor_A.contrasts,2), ...
            'FactorB', 1:size(STATS.subject_results.subject_1.(band_fields{1}).factor_B.contrasts,2), ...
            'FactorAB', 1:size(STATS.subject_results.subject_1.(band_fields{1}).factor_AxB.contrasts,2));
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_A.contrasts)
            disp('FactorB'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_B.contrasts)
            disp('FactorAB'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_AxB.contrasts)
        end
        
    elseif strcmp(STATS.design,'w');
        
        options = struct('FactorA', 1:size(STATS.subject_results.subject_1.(band_fields{1}).factor_A.contrasts,2));
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_A.contrasts)
        end
        
    elseif strcmp(STATS.design,'bw'); % considered not factorial in single-subject cases
        band_fields=fieldnames(STATS.subject_results.Factor_A1.subject_1);
        
        options = struct('FactorA', 1:size(STATS.subject_results.Factor_A1.subject_1.(band_fields{1}).factor_A.contrasts,2));
        
        if infodisplay
            disp('j level labels'); disp(STATS.jlabels)
            disp('k level labels'); disp(STATS.klabels)
            disp('contrasts'); disp(STATS.subject_results.Factor_A1.subject_1.(band_fields{1}).factor_A.contrasts)
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
            disp('FactorA'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_A.contrasts)
            disp('FactorB'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_B.contrasts)
            disp('FactorAB'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_AxB.contrasts)
        end
        
    elseif strcmp(STATS.design,'w');
        
        options = struct('FactorA', []);
        
        if infodisplay
            disp('Condition names'); disp(STATS.condnames)
            disp('FactorA'); disp(STATS.subject_results.subject_1.(band_fields{1}).factor_A.contrasts)
        end
        
    elseif strcmp(STATS.design,'bw'); % considered not factorial in single-subject cases
        band_fields=fieldnames(STATS.subject_results.Factor_A1.subject_1);
        
        options = struct('FactorA', []);
        
        if infodisplay
            disp('j level labels'); disp(STATS.jlabels)
            disp('k level labels'); disp(STATS.klabels)
            disp('contrasts'); disp(STATS.subject_results.Factor_A1.subject_1.(band_fields{1}).factor_A.contrasts)
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
            %svg
            if strcmp(options.savesvg,'yes'); buildsvg(h(1),conds_ax,STATS,'ersp_subject_subplot_1.svg','jet',options.timeplot); end
            if ersp==1; colormap(jet); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
            if itc==1; colormap(flipud(hot)); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
            
            hsub(2)=subplot(4,1,2);
            h(2)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plot2nd(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
            %svg
            if strcmp(options.savesvg,'yes'); buildsvg(h(2),conds_ax,STATS,'ersp_subject_subplot_2.svg','jet',options.timeplot); end
            if ersp==1; colormap(jet); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
            if itc==1; colormap(flipud(hot)); caxis(conds_ax); cbfreeze(colorbar); freezeColors; end
            
            hsub(3)=subplot(4,1,3);
            h(3)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(plotdiff(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
            %svg
            if strcmp(options.savesvg,'yes'); buildsvg(h(3),conds_axdiff,STATS,'ersp_subject_subplot_3.svg','jet',options.timeplot); end
            colormap(jet); caxis(conds_axdiff); cbfreeze(colorbar); freezeColors;
            
            set(allchild(gca),'buttondownfcn',{@mouseclick_callback, STATS, [leg1st,'-',leg2nd], options.timeplot, options.FactorA(i), 'A'});
            
            hsub(4)=subplot(4,1,4);
            h(4)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pval_prop(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
            %svg
            if strcmp(options.savesvg,'yes'); buildsvg(h(4),[0 1],STATS,'ersp_subject_subplot_4.svg','bonesub',options.timeplot); end
            colormap(hsub(4),flipud(hot)); caxis([0 1]); hb=cbfreeze(colorbar);  %freezeColors; %set(hb,'YTick',[0 1]);
            
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
            
            % reset
            pval_gather=zeros(length(band_fields), length(STATS.TF_times));
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
                colormap(hsub(4),flipud(hot)); caxis([0 1]); hb=cbfreeze(colorbar);  %freezeColors; %set(hb,'YTick',[0 1]);
                
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
                
                % reset
                pval_gather=zeros(length(band_fields), length(STATS.TF_times));

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
                colormap(hsub(4),flipud(hot)); caxis([0 1]); hb=cbfreeze(colorbar);  %freezeColors; %set(hb,'YTick',[0 1]);
                
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
                
                % reset
                pval_gather=zeros(length(band_fields), length(STATS.TF_times));
                
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
                colormap(hsub(2),flipud(hot)); caxis([0 1]); hb=cbfreeze(colorbar);  %freezeColors; %set(hb,'YTick',[0 1]);
                
                % add legend and font
                set(gca,'FontSize',10);
                h=axes('visible','off');
                lh=title(['(',leg1st,')','-','(',leg2nd,')'],'parent',h,'visible','on','Position',[.5 1.05 0]);
                % get rid of subscripts that occur when there are underscores
                set(lh,'Interpreter', 'none');
                axis tight
                grid on
                
                % clear variables before next iteration
                clear plot1st plot2nd plotdiff leg1st leg2nd leg2ndlist leg1stlist lh
                
                % reset
                pval_gather=zeros(length(band_fields), length(STATS.TF_times));
                
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
                    
                    % sub loop
                    for s=1:length(subnames);
                        
                        % stats
                        for q=1:STATS.freqbins;
                            pval_tmp(q,:)=(STATS.subject_results.(factnames{v}).(subnames{s}).(band_fields{q}).factor_A.pval(options.FactorA(i),:)<=.05);
                        end
                        
                        pval_gather=pval_gather+pval_tmp;
                    end
                    
                    % pval proportions
                    pval_prop=pval_gather./length(subnames);
                    
                    for j=1:STATS.levels(2);
                        
                        % get the condition waveforms
                        if STATS.subject_results.(factnames{1}).subject_1.(fields{1}).factor_A.contrasts(j,options.FactorA(i))==1
                            %plot1st(k,:)=STATS.condwaves_trim(j,:);
                            c=load(['group_TFwaves_',STATS.savestring,'_',STATS.condnames{y}]);
                            plot1st(:,:,k)=c.condition;
                            leg1stlist{k}=STATS.klabels{j};
                            %k=k+1;
                            %lab1stlist{k}=STATS.klabels{j};
                            k=k+1;
                        elseif STATS.subject_results.(factnames{1}).subject_1.(fields{1}).factor_A.contrasts(j,options.FactorA(i))==-1
                            %plot2nd(m,:)=STATS.condwaves_trim(j,:);
                            c=load(['group_TFwaves_',STATS.savestring,'_',STATS.condnames{y}]);
                            plot2nd(:,:,m)=c.condition;
                            leg2ndlist{m}=STATS.klabels{j};
                            %lab2ndlist{m}=STATS.klabels{j};
                            m=m+1;
                        end
                        
                        
                        
                        y=y+1;
                    end
                    
                    % reduce, if needed, the condition waveforms
                    % should this be SUM or MEAN if pooling across levels?
                    plot1st=mean(plot1st,3);
                    plot2nd=mean(plot2nd,3);
                    plotdiff=plot1st-plot2nd;
                    
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
                    
                    set(allchild(gca),'buttondownfcn',{@mouseclick_callback, STATS, [STATS.jlabels{v},' @ ',leg1st,' - ',leg2nd], options.timeplot, options.FactorA(i), ['A', num2str(v)]});
                    
                    hsub(4)=subplot(4,1,4);
                    h(4)=surf(STATS.TF_times(options.timeplot),STATS.TF_freqs,double(pval_prop(:,options.timeplot)),'facecolor','interp','linestyle','none'); axis tight; view(0,90);
                    colormap(hsub(4),flipud(hot)); caxis([0 1]); hb=cbfreeze(colorbar);  %freezeColors; %set(hb,'YTick',[0 1]);
                    
                    % concatenate, if needed, the legend lables
                    leg1st=strjoin_statslab(leg1st,'+');
                    leg2nd=strjoin_statslab(leg2nd,'+');
                    titlestr=[STATS.jlabels{v},' @ ',leg1st,' - ',leg2nd];
                    
                    % add legend and font
                    set(gca,'FontSize',10);
                    h=axes('visible','off');
                    lh=title(titlestr,'parent',h,'visible','on','Position',[.5 1.05 0]);
                    
                    % get rid of subscripts that occur when there are underscores
                    set(lh,'Interpreter', 'none');
                    axis tight
                    grid on
                    
                    % clear variables before next iteration
                    clear plot1st plot2nd plotdiff leg1st leg2nd leg2ndlist leg1stlist lh
                    
                    % reset
                    pval_gather=zeros(length(band_fields), length(STATS.TF_times));
                end
                
                %y=y+STATS.levels(2)-1;
            end
            
        end
end

disp('******* Saving STATS structure *******')
save(['STATS_',STATS.savestring,'.mat'],'STATS');



    function mouseclick_callback(gcbo, eventdata, STATS, titl, timeplot, oncon, onfact)
        
        diffcol=[0 0 0];
        CIcol=[.5 .5 .5];
        %get the point that was clicked on
        cP = get(gca,'Currentpoint');
        % ms_plot = cP(1,1);
        freqclick = cP(1,2);
        [V I]=min(abs(STATS.TF_freqs-freqclick));
        
        
        if ~strcmp(STATS.design, 'bw')
            flds=fieldnames(STATS.sample_results);
            subflds=fieldnames(STATS.subject_results);
            
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
                %axis tight
                ylim([-15 15]);
                grid on
            end
            
            htit = axes('visible','off');
            tit=title([titl, ' -> ', num2str(STATS.TF_freqs(I)), ' Hz'], 'parent', htit, 'visible', 'on');
            % get rid of subscripts that occur when there are underscores
            set(tit,'Interpreter', 'none');
            axis tight
            grid on
            
        elseif strcmp(STATS.design, 'bw')
            
            fact=['Factor_', onfact];
            flds=fieldnames(STATS.sample_results);
            subflds=fieldnames(STATS.subject_results.(fact));
            
            figure;
            for qq=1:length(subflds);
                
                % extract CI
                CI_lower(1,:)=STATS.subject_results.(fact).(subflds{qq}).(flds{I}).factor_A.CI{oncon,1}(1,timeplot);
                CI_upper(1,:)=STATS.subject_results.(fact).(subflds{qq}).(flds{I}).factor_A.CI{oncon,1}(2,timeplot);
                
                % extract difference wave
                pdiff=STATS.subject_results.(fact).(subflds{qq}).(flds{I}).factor_A.test_stat(oncon,timeplot);
                
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
            tit=title([titl, ' -> ', num2str(STATS.TF_freqs(I)), ' Hz'], 'parent', htit, 'visible', 'on');  
            % get rid of subscripts that occur when there are underscores
            set(tit,'Interpreter', 'none');
            axis tight
            grid on
            

        end
    end

end









