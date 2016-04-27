%[sim]=simulation_realdata(nboot)
%% get some erp data


% [fnames]=uigetfile('.mat','get','MultiSelect','On');
% 
% for i=1:length(fnames);
% load(fnames{i});
% ERP(i,:)=mean(data,1);
% end
% save('ERPdata.mat','ERP')

load ERPdata.mat
nboot=5000;
%% one sample ttest based on null

u=mean(ERP,1);
utrim=trimmean(ERP,40);
[rows cols]=size(ERP);
T=zeros(nboot,cols);
count_bv=zeros(12,1536);
sum_bv=zeros(12,1536);
for i=1:nboot;
    
    
    %for k=1:1536
    % take two samples
    bv=randi(rows,1,rows);
    grp1=ERP(bv,:);
    
    % do a regular ttest
    E=mean(grp1)-u;
    SE=std(grp1)/(sqrt(rows));
    T(i,:)=E./SE;
    
    
    % do a trimmed ttest w winsorizing
    Et=trimmean(grp1,40)-utrim;
    sw=std(winsorize(grp1,.2));
    SEt=sw./(.6*sqrt(rows));
    Ttrim(i,:)=Et./SEt;
    
      [C] = countmember(1:12,bv);

     if i==1;
         sum_bv=C;
         
     else 
         sum_bv=C+sum_bv;
     end

end

% create normal curve;
R=normrnd(0,1,1,100000);

% desity of norm distribution
[bandwidth,f,x,cdf]=kde(R,1000);

for i=1:1536;
    % desity of estimated distribution
    %[bandwidth,ft(:,i),xt(i,:),cdf]=kde(T(:,i),2^12); %2^12
    
    [ft(:,i),xt(i,:)] = ksdensity(T(:,i));
end

for i=1:1536;
    % desity of trimmed distribution
    %[bandwidth,ft_trim(:,i),xt_trim(i,:),cdf]=kde(Ttrim(:,i),2^12);
    
    [ft_trim(:,i),xt_trim(i,:)] = ksdensity(Ttrim(:,i));
end
%% CIs

%CIs around contrast differences
CIlow=round(.05*nboot/2)+1;
CIup=nboot-CIlow-1;

SE=std(ERP)/(sqrt(rows));
% assumed ttest
for i=1:1536;
    %temp=sort(T(:,i));
    %conflow=temp(CIlow);
    %confup=temp(CIup);
    
   
   conflow2_trad(i)=u(i)-(-1.96*SE(i));
   confup2_trad(i)=u(i)-(1.96*SE(i));
end

SE=std(ERP)/(sqrt(rows));
% boot ttest no trimming
for i=1:1536;
    temp=sort(T(:,i));
    conflow(i)=temp(CIlow);
    confup(i)=temp(CIup);
    
   conflow2(i)=u(i)-(conflow(i)*SE(i));
   confup2(i)=u(i)-(confup(i)*SE(i));
end

% trimmed t-test
    sw=std(winsorize(ERP,.2),1);
    SEt=sw./(.6*sqrt(rows));
for i=1:1536;
    temp=sort(Ttrim(:,i));
    conflow_trim(i)=temp(CIlow);
    confup_trim(i)=temp(CIup);
    

   conflow2_trim(i)=utrim(i)-(conflow_trim(i)*SEt(i));
   confup2_trim(i)=utrim(i)-(confup_trim(i)*SEt(i));
    
end

%% figure
colors = hsv(12);
v = VideoWriter('real_data_sim_stats.avi');
open(v);

fHand = figure;
set(fHand,'units','normalized','position', [ 0.1300    0.1100    0.7750    0.8150]);
for k=1:3:1536;
    
    h4=subplot(4,1,4);
    h4d=plot(STATS.xtimes,u);
    hold on
    h4ds=scatter(STATS.xtimes(k),u(k),50, 'k', 'Filled');
    ylabel('GRAND ERP','fontsize', 18)
    
    h3=subplot(4,1,3);
    for i = 1:12;
        h3d(i)=bar(i,ERP(i,k), 'facecolor', colors(i,:));
        hold on
    end
    ylim([-5 5]);
    set(gca, 'XTick', 1:12);
    ylabel('SUBJECT UV')
    
    h2=subplot(4,1,2);
    
    h2d=plot(linspace(-1.96,1.96,2),zeros(1,2),'-ko','MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on
    h2d2=plot(linspace(conflow(k),confup(k),2),zeros(1,2)+.05,'-ro','MarkerFaceColor','r','MarkerEdgeColor','r');
    hold on
    h2d2_trim=plot(linspace(conflow_trim(k),confup_trim(k),2),zeros(1,2)+.1,'-bo','MarkerFaceColor','b','MarkerEdgeColor','b');
    %hold on
    %h2d2_trad=plot(linspace(confup_trad(k),conflow_trad(k),2),zeros(1,2)+.05,'-ko','MarkerFaceColor','k','MarkerEdgeColor','k');
    xlim([-4 4])
    ylim([-.05 .15])
    set(gca, 'YTick', []);
    
    ylabel('ASSUMED & ESTIMATED CIs')
    
    h1=subplot(4,1,1);
    plot(x,f,'k','linewidth',2)
    xlim([-4 4])
    ylim([0 .5])
    
    hold(gca,'on');
    h1b=plot(xt(k,:),ft(:,k),'r', 'linewidth',2);
    hold(gca,'on');
    h1b_trim=plot(xt_trim(k,:),ft_trim(:,k),'b', 'linewidth',2);
    ylabel('ASSUMED & ESTIMATED T')
    ht=title([num2str(STATS.xtimes(k)),' ms, ', num2str(k), ' TF']);
    
    %waitforbuttonpress;
    frame = getframe(gcf);
    writeVideo(v,frame);
    delete([ht ht h3d  h1b h4d h4ds h2d h2d2 h2d2_trim h1b_trim]);
end
close(v);
% 
% for i=1:1536;
% %figure; plot(x,f,'k','linewidth',2)
% hold(gca,'on');
% h=plot(xt(i,:),ft(:,i),'r', 'linewidth',2);
% ht=title([num2str(STATS.xtimes(i)), ' ms']);
% waitforbuttonpress;
% delete([h ht]);
% end
% close(gcf);
% %saveas(gcf,['real_data_sim_', num2str(i),'.png']);
% %close(gcf);
% %end
% 
% v = VideoWriter('real_data_sim.avi');
% open(v);
% fid=figure;
% xlim([-20 20])
% for i=1:3:1536;
% figure(fid);
% hold(gca,'on');
% plot(x,f,'k','linewidth',2)
% hold on; plot(xt(i,:),ft(:,i),'r', 'linewidth',2);
% frame = getframe(gcf);
% writeVideo(v,frame);
% hold(gca,'on');
% plot(xt(i,:),ft(:,i),'w', 'linewidth',2);
% %saveas(gcf,['real_data_sim_', num2str(i),'.png']);
% %close(gcf);
% end
% close(v);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
