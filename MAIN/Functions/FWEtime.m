function [STATS]=FWEtime(STATS,alpha,nboot,wherefrom,varargin)

% Benjamin Hochberg. Described in Groppe (mass univariate)

% check for type of analysis

TF=0;
try
    if any(strcmp({'ersp' 'itc'},STATS.measure));
        TF=1;
    end
catch
    
    if strcmp(wherefrom,'subjectTF');
        TF=1;
    end
end

switch TF
    
    case 0
        
        switch wherefrom
            
            case 'group'
                
                % get factors
                fs=fieldnames(STATS.sample_results);
                
                for j=1:length(fs);
                    
                    concol=size(STATS.sample_results.(fs{j}).contrasts,2);
                    pvec=zeros(concol,STATS.numpnts);
                    
                    % get vector of pvalues
                    for i=1:concol;
                        pvec(i,:)=STATS.sample_results.(fs{j}).pval(i,:);
                    end
                    
                    %num rows of pvalues (so number of contrasts)
                    [rc cc]=size(pvec);
                    
                    % concatenate results from each contrast
                    pvec=reshape(pvec', 1, rc*cc);
                    
                    m=length(pvec);
                    
                    % sort pvals
                    [p ind]=sort(pvec);
                    
                    % find k (largest val for which pi <= (i/m)*alpha
                    im=(1:m)/m;
                    ima=im*alpha;
                    pcomp=p<=ima;
                    k=max(find(pcomp==1));
                    
                    % so, 1 through k are rejected
                    % ones given to all points that will be rejected
                    pcomp(1:k)=1;
                    padj=ind(pcomp);
                    pnew=zeros(1,m);
                    pnew(padj)=1;
                    
                    % parse
                    significant=reshape(pnew, length(STATS.xtimes), rc)';

                    pcrit=ima(ind);
                    
                    % CI inds (based on the lowest value found in the adjusted critical
                    % vector (as in Wilcox's toolbox)
                    critCI=min(pcrit);
                    CIlow=round(critCI*nboot/2)+1;
                    CIup=nboot-CIlow-1;
                    
                    % finding quantiles for CIs
                    for k=1:rc;
                        diffs=sort(STATS.sample_results.(fs{j}).diffs{k});
                        STATS.sample_results.(fs{j}).CI{k}(1,:)=diffs(CIlow,:);
                        STATS.sample_results.(fs{j}).CI{k}(2,:)=diffs(CIup,:);
                        STATS.sample_results.(fs{j}).sig_pvals=significant(k,:);
                    end
                    %pcrit=reshape(pcrit, STATS.numpnts, rc)';
                    
                    clear pvec
                end
                
                for j=1:length(fs);
                    % remove diff arrays
                    [STATS.sample_results.(fs{j})]=rmfield(STATS.sample_results.(fs{j}),'diffs');
                end
                
                
            case 'subject'
                
                %%%%%%%%%
                % STATS is not really the stats structure here, it's a subject's results structure
                %%%%%%%%%
                
                % get factors and band names
                fac=fieldnames(STATS); % factors
                time=size(STATS.(fac{1}).pval,2);
                
                % loop factors
                for fa=1:length(fac);
                    concol=size(STATS.(fac{fa}).contrasts,2);
                    pvec=zeros(concol,time);
                    
                    % get vector of pvalues
                    for i=1:concol;
                        pvec(i,:)=STATS.(fac{fa}).pval(i,:);
                    end

                    %num rows of pvalues (so number of contrasts)
                    [rc cc]=size(pvec);
                    
                    % concatenate results from each contrast
                    pvec=reshape(pvec', 1, rc*cc);
                    
                    m=length(pvec);
                    
                    % sort pvals
                    [p ind]=sort(pvec);
                    
                    % find k (largest val for which pi <= (i/m)*alpha
                    im=(1:m)/m;
                    ima=im*alpha;
                    pcomp=p<=ima;
                    k=max(find(pcomp==1));
                    
                    % ones given to all points that will be rejected
                    pcomp(1:k)=1;
                    padj=ind(pcomp);
                    pnew=zeros(1,m);
                    pnew(padj)=1;
                    
                    % parse
                    significant=reshape(pnew, cc, rc)';
                    
                    pcrit=ima(ind);
                    critCI=min(pcrit);
                    CIlow=round(critCI*nboot/2)+1;
                    CIup=nboot-CIlow-1;
                    
                    % finding quantiles for CIs
                    for k=1:rc;
                        diffs=sort(varargin{fa}{k});
                        STATS.(fac{fa}).CI{k}(1,:)=diffs(CIlow,:);
                        STATS.(fac{fa}).CI{k}(2,:)=diffs(CIup,:);
                        STATS.(fac{fa}).sig_pvals(k,:)=significant(k,:);
                    end
                    
                    clear pvec
                end
        end
        
        
    case 1
        
        switch wherefrom
            
            case 'group'
                
                % clean up structure a bit
                STATS.tmpA=STATS.sample_results.factor_A;
                [STATS.sample_results]=rmfield(STATS.sample_results,'factor_A');
                try
                    STATS.tmpB=STATS.sample_results.factor_B;
                    [STATS.sample_results]=rmfield(STATS.sample_results,'factor_B');
                    STATS.tmpAB=STATS.sample_results.factor_AxB;
                    [STATS.sample_results]=rmfield(STATS.sample_results,'factor_AxB');
                catch
                end
                
                % get factors and band names
                fs=fieldnames(STATS.sample_results); % bands
                ft=fieldnames(STATS.sample_results.(fs{1})); % factors
                tmpstr={'tmpA','tmpB','tmpAB'}; % tmp structures
                
                % loop factors
                for fa=1:length(ft);
                    concol=size(STATS.sample_results.(fs{1}).(ft{fa}).contrasts,2);
                    pvec=zeros(concol,STATS.timesout,length(fs));
                    
                    for j=1:length(fs);
                        % get vector of pvalues
                        for i=1:size(STATS.sample_results.(fs{1}).(ft{fa}).contrasts,2);
                            pvec(i,:,j)=STATS.sample_results.(fs{j}).(ft{fa}).pval(i,:);
                        end
                    end
                    %num rows of pvalues (so number of contrasts)
                    [rc cc pp]=size(pvec);
                    
                    % concatenate results from each contrast
                    pvec=reshape(pvec, 1, rc*cc*pp);
                    
                    m=length(pvec);
                    
                    % sort pvals
                    [p ind]=sort(pvec);
                    
                    % find k (largest val for which pi <= (i/m)*alpha
                    im=(1:m)/m;
                    ima=im*alpha;
                    pcomp=p<=ima;
                    k=max(find(pcomp==1));
                    
                    % ones given to all points that will be rejected
                    pcomp(1:k)=1;
                    padj=ind(pcomp);
                    pnew=zeros(1,m);
                    pnew(padj)=1;
                    
                    % parse
                    significant=reshape(pnew, rc, cc, pp);
                    pcrit=ima(ind);
                    critCI=min(pcrit);
                    CIlow=round(critCI*nboot/2)+1;
                    CIup=nboot-CIlow-1;
                    
                    % if you wnat to add crit vals to stats stucture
                    %pcrit=reshape(pcrit, rc,cc,pp);
                    
                    % finding quantiles for CIs
                    for k=1:rc;
                        diffs=sort(STATS.(tmpstr{fa}).diffs{k});
                        for b=1:pp;
                            STATS.sample_results.(fs{b}).(ft{fa}).CI{k}(1,:)=diffs(CIlow,:,b);
                            STATS.sample_results.(fs{b}).(ft{fa}).CI{k}(2,:)=diffs(CIup,:,b);
                            STATS.sample_results.(fs{b}).(ft{fa}).sig_pvals(k,:)=significant(k,:,b);
                            %STATS.sample_results.(fs{b}).factor_A.critval(k,:)=pcrit(k,:,pp);
                        end
                    end
                    
                    clear pvec
                end
                
                % clean structure
                try
                    [STATS]=rmfield(STATS,'tmpA');
                    [STATS]=rmfield(STATS,'tmpB');
                    [STATS]=rmfield(STATS,'tmpAB');
                catch
                end
                
                
            case 'subjectTF'
                
                %%%%%%%%%
                % STATS is not really the stats structure here, it's the results structure
                %%%%%%%%%
                
                % get factors and band names
                bd=fieldnames(STATS); % bands
                fac=fieldnames(STATS.(bd{1})); % factors
                
                time=size(STATS.(bd{1}).(fac{1}).pval,2);
                
                % loop factors
                for fa=1:length(fac);
                    concol=size(STATS.(bd{1}).(fac{fa}).contrasts,2);
                    pvec=zeros(concol,time,length(bd));
                    
                    for j=1:length(bd);
                        % get vector of pvalues
                        for i=1:concol;
                            pvec(i,:,j)=STATS.(bd{j}).(fac{fa}).pval(i,:);
                        end
                    end
                    %num rows of pvalues (so number of contrasts)
                    [rc cc pp]=size(pvec);
                    
                    % concatenate results from each contrast
                    pvec=reshape(pvec, 1, rc*cc*pp);
                    
                    m=length(pvec);
                    
                    % sort pvals
                    [p ind]=sort(pvec);
                    
                    % find k (largest val for which pi <= (i/m)*alpha
                    im=(1:m)/m;
                    ima=im*alpha;
                    pcomp=p<=ima;
                    k=max(find(pcomp==1));
                    
                    % ones given to all points that will be rejected
                    pcomp(1:k)=1;
                    padj=ind(pcomp);
                    pnew=zeros(1,m);
                    pnew(padj)=1;
                    
                    % parse
                    significant=reshape(pnew, rc, cc, pp);
                    pcrit=ima(ind);
                    critCI=min(pcrit);
                    CIlow=round(critCI*nboot/2)+1;
                    CIup=nboot-CIlow-1;
                    
                    % if you wnat to add crit vals to stats stucture
                    %pcrit=reshape(pcrit, rc,cc,pp);
                    
                    % finding quantiles for CIs
                    for k=1:rc;
                        diffs=sort(varargin{fa}{k});
                        for b=1:pp;
                            STATS.(bd{b}).(fac{fa}).CI{k}(1,:)=diffs(CIlow,:,b);
                            STATS.(bd{b}).(fac{fa}).CI{k}(2,:)=diffs(CIup,:,b);
                            STATS.(bd{b}).(fac{fa}).sig_pvals(k,:)=significant(k,:,b);
                        end
                    end

                    clear pvec
                end
                
        end
end


end




