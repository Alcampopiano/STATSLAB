function [conflow confup] = motherCI(dataCI, connum, numpnts, nsamp, alpha)
%{
this function finds the CI of the CI bounds. I call this the mother CI. IN
practice I have found it to be very wide (conservative) compared to other
techniques that still account for variability accross subjects.
%}

% find CIs for input matrices (dataCI)

% Rom/Hoch corrections for FWE.

conflow=zeros(connum,numpnts);
confup=zeros(connum,numpnts);

%{
% generalized p value
pval=zeros(1,connum);
for i=1:connum;
    pval(i)=(sum(psihat(i,:)>0)+.5*sum(psihat(i,:)==0))/nboot;
    pval(i)=min(pval(i),1-pval(i));
end
pvalgen=2*pval;
%}

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


% SHOULD THE MOTHER CIS BE CORRECTED FOR FWE, BEACUSE EACH BOOTRUN THEY
% WERE ALREADY CORRECTED...??
% FWE corrections to p values
dvec=alpha./(2*(1:connum));
dvec=2*dvec;

%{
zvec=dvec(1:connum);
[~, pvalgen_inds]=sort(0-pvalgen);
pcrit(pvalgen_inds)=zvec;
%}

%CIs around contrast differences
CIlow=round(dvec(connum)*nsamp/2)+1;
CIup=nsamp-CIlow-1;


for i=1:connum;
    templow=sort(dataCI{1}{i});
    conflow(i,:)=templow(CIlow,:);
    tempup=sort(dataCI{2}{i});
    confup(i,:)=tempup(CIup,:);
    clear templow tempup
end


end

