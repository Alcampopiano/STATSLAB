function [con] = con1way(numconds)
% create all possible pairwise contrast coefficients for 1-way ANOVA


numcols=(numconds^2-numconds)/2;
con=zeros(numconds,numcols);
q=0;
for i=1:numconds;
    for j=1:numconds;
        if i<j;
            q=q+1;
            con(i,q)=1;
            con(j,q)=0-1;
        end
    end
end

end

