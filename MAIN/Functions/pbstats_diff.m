function [psihat_stat pvalgen pcrit conflow confup] = pbstats_diff(data, con, nboot, alpha, FWE)
%bla bla bla

% uses con2way (which builds appropriate linear contrasts) to do percentile
% bootstrap tests on differences defined by the contrast coefficients. Used
% Rom/Hoch corrections for FWE.


% Expects data in the form of bootstrapped ERPs, for each subject, for each condition
% (i.e., one file for each anova cell, for each subject, containing thousands of bootstrapped ERPs.
% Differences based on contrast coefficients are calculated here using the input waveforms.


% preallocate
[rowcol, concol]=size(con);
connum=concol;
%psihat=zeros(connum,nboot);
conflow=zeros(connum,1);
confup=zeros(connum,1);
%[rowdat coldat]=size(data);

% create diferences based on the contrast coefficients
%datcon=data*con;
datcon=data;

% take the mean difference (mostly to plot it later)
psihat_stat=mean(datcon,1);

% define random indices
%bootinds=randi(N,nboot,N);

%boostrap values
%for i = 1:nboot;
%    psihat(:,i)=trimmean(datcon(bootinds(i,:),:),40,1);
%end


psihat=datcon';


% generalized p value
pval=zeros(1,connum);
for i=1:connum;
    pval(i)=(sum(psihat(i,:)>0)+.5*sum(psihat(i,:)==0))/nboot;
    pval(i)=min(pval(i),1-pval(i));
end
pvalgen=2*pval;

% FWE corrections Rom's method (only suitable for n<= 10 hypotheses
% (contrasts to test)

%if alpha==.05;
%   dvec=[.025 .025 .0169 .0127 .0102 .00851 .0073 .00639 .00568 .00511];
%    if connum>10;
%        connumvec=11:connum;
%        avec=.05./connumvec;
%        dvec=[dvec avec];
%    end
%end

%if alpha==.01;
%   dvec=[.005 .005 .00334 .00251 .00201 .00167 .00143 .00126 .00112 .00101];
%    if connum>10;
%        connumvec=11:connum;
%        avec=.01./connumvec;
%        dvec=[dvec avec];
%    end
%end

% this makes FWE methods easily expandable in the future

switch FWE
    
    case 'Rom'
        
        % FWE corrections to p values
        dvec=alpha./(2*(1:connum));
        dvec=2*dvec;
        zvec=dvec(1:connum);
        [~, pvalgen_inds]=sort(0-pvalgen);
        pcrit(pvalgen_inds)=zvec;
        
        %CIs around contrast differences
        CIlow=round(dvec(connum)*nboot/2)+1;
        CIup=nboot-CIlow-1;
        
        for i=1:connum;
            temp=sort(psihat(i,:));
            conflow(i,1)=temp(CIlow);
            confup(i,1)=temp(CIup);
        end
        
    case 'none'
        
        dvec=zeros(1:connum);
        dvec(:)=alpha; % dvec is just alpha connum times in a row
        pcrit=alpha;
        
        %CIs around contrast differences
        CIlow=round(dvec(connum)*nboot/2)+1;
        CIup=nboot-CIlow-1;
        
        for i=1:connum;
            temp=sort(psihat(i,:));
            conflow(i,1)=temp(CIlow);
            confup(i,1)=temp(CIup);
        end
end



end

