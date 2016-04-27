function X=gh_simulator(g,h)

Z=normrnd(0,1,1,100000);


% desity of pops distribution
%[bandwidth,f,x,cdf]=kde(Z,5000);

%figure;
%plot(x,f,'k', 'LineWidth' ,2);
if g~=0;
    
    t1=exp(g*Z)-1;
    t2=g;
    t3=exp(h*Z.^2/2);
    X=(t1/t2).*t3;
    
elseif g==0;
    
    X=Z.*exp(h*Z.^2/2);
end

%hold on
% desity of pops distribution
%[bandwidth,f,x,cdf]=kde(X,5000);
%plot(x,f,'r', 'LineWidth' ,2);

%xlim([-5 5]);
end

