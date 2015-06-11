function [] = splash_fade()

% load images
lowess_hot=imread('lowess_hotsc.png');
slopeCI_example=imread('SlopeCI_example.png');
Subject_BS_128=imread('Subject_BS_128.png');
textstats=imread('statstext.png');
alpha_lowhi=linspace(0,1,10);
alpha_hilow=linspace(1,0,10);

% if ims are too big, supress warning
warning off

%%%%%%%%
% create the gui figure
% sz = [891 1068]; % figure size
% screensize = get(0,'ScreenSize');
% xpos = ceil((screensize(3)-sz(2))/2);
% ypos = ceil((screensize(4)-sz(1))/2);
%
% h = figure(...
%     'visible','off',...
%     'position',[xpos, ypos, sz(2), sz(1)],...
%     'units','pixels',...
%     'renderer','OpenGL',...
%     'MenuBar','none',...
%     'PaperPositionMode','auto',...
%     'Name','Data Input Dialog',...
%     'NumberTitle','off',...
%     'Tag','gui',...
%     'Resize','off');

%%%%%%%
% this aviods a weird positioning bug
figure('visible','off');
imshow(Subject_BS_128);
set(gcf,'menubar','none')
%%%%%%%

h=imshow(Subject_BS_128);
set(gcf,'visible','on');
set(gcf, 'CloseRequestFcn', '')
set(gcf,'menubar','none')
set(gcf,'NumberTitle','off');
set(h, 'AlphaData', alpha_lowhi(2))
for i=1:length(alpha_hilow)-2;
    set(h, 'AlphaData', alpha_lowhi(i+1))
    pause(.0005)
end
pause(.5)

set(h, 'AlphaData', alpha_hilow(2))
for i=1:length(alpha_hilow)-1;
    set(h, 'AlphaData', alpha_hilow(i+1))
    pause(.0005)
end
pause(.5)


h=imshow(slopeCI_example);
set(gcf,'menubar','none')
set(gcf,'NumberTitle','off');
set(h, 'AlphaData', alpha_lowhi(2))
for i=1:length(alpha_hilow)-2;
    set(h, 'AlphaData', alpha_lowhi(i+1))
    pause(.0005)
end
pause(.5)
set(h, 'AlphaData', alpha_hilow(2))
for i=1:length(alpha_hilow)-1;
    set(h, 'AlphaData', alpha_hilow(i+1))
    pause(.0005)
end
pause(.5)
h=imshow(lowess_hot);

set(gcf,'menubar','none')
set(gcf,'NumberTitle','off');
set(h, 'AlphaData', alpha_lowhi(2))
for i=1:length(alpha_hilow)-2;
    set(h, 'AlphaData', alpha_lowhi(i+1))
    pause(.0005)
end
pause(.5)
set(h, 'AlphaData', alpha_hilow(2))
for i=1:length(alpha_hilow)-1;
    set(h, 'AlphaData', alpha_hilow(i+1))
    pause(.0005)
end

pause(.5)
h=imshow(textstats);
set(gcf,'menubar','none')
set(gcf,'NumberTitle','off');
set(h, 'AlphaData', alpha_lowhi(2))
for i=1:length(alpha_hilow)-1;
    set(h, 'AlphaData', alpha_lowhi(i))
    pause(.005)
end

warning on
waitforbuttonpress
delete(gcf);



