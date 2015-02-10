function [conA conB conAB] = con2way(jlvls, klvls);

% create linear contrasts for two-way

% contast for factor A
Jsize=(jlvls^2-jlvls)/2;
Ksize=(klvls^2-klvls)/2;

conA=zeros(jlvls*klvls,Jsize);
j=0;
for i=1:jlvls;
    for k=1:jlvls;
        if i<k;
            j=j+1;
            tempmat=zeros(jlvls,klvls);
            tempmat(i,:)=1;
            tempmat(k,:)=0-1;
            conA(:,j)=reshape(tempmat',jlvls*klvls,1);
        end
    end
end


% contast for factor B
conB=zeros(jlvls*klvls,Ksize);
j=0;
for i=1:klvls;
    for k=1:klvls;
        if i<k;
            j=j+1;
            tempmat=zeros(jlvls,klvls);
            tempmat(:,i)=1;
            tempmat(:,k)=0-1;
            conB(:,j)=reshape(tempmat',jlvls*klvls,1);
        end
    end
end

% contast for interaction AxB
conAB=zeros(jlvls*klvls,Jsize*Ksize);
j=0;
for i=1:jlvls;
    for k=1:jlvls;
        if i<k;
            for ii=1:klvls;
                for kk=1:klvls;
                    if ii<kk;
                        j=j+1;
                        tempmat=zeros(jlvls,klvls);
                        tempmat(i,ii)=1;
                        tempmat(i,kk)=0-1;
                        tempmat(k,ii)=0-1;
                        tempmat(k,kk)=1;
                        conAB(:,j)=reshape(tempmat',jlvls*klvls,1);
                    end
                end
            end
        end
    end
end
end

