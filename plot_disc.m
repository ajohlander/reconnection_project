function [fLMN] = plot_disc(b1,b2,b3,b4,t)
%PLOT_DISC Summary of this function goes here
%   Detailed explanation goes here



%---------------Figure-------------------
fLMN = irf_plot(3,'newfigure');
set(gcf,'PaperUnits','centimeters')
xSize = 15; ySize = 8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop

h = zeros(1,3);

h(1) = irf_panel('L');
h(2) = irf_panel('M');
h(3) = irf_panel('N');
hold(h(1))
hold(h(2))
hold(h(3))

%t = zeros(1,4); %for no time-shift

% Maximum
plot(h(1),b1(:,1)-t(1),b1(:,2),'k')   %C1
plot(h(1),b2(:,1)-t(2),b2(:,2),'r')   %C2
plot(h(1),b3(:,1)-t(3),b3(:,2),'g')   %C3
plot(h(1),b4(:,1)-t(4),b4(:,2),'b')   %C4

% Intermediate
plot(h(2),b1(:,1)-t(1),b1(:,3),'k')   %C1
plot(h(2),b2(:,1)-t(2),b2(:,3),'r')   %C2
plot(h(2),b3(:,1)-t(3),b3(:,3),'g')   %C3
plot(h(2),b4(:,1)-t(4),b4(:,3),'b')   %C4

% Minimum
plot(h(3),b1(:,1)-t(1),b1(:,4),'k')   %C1
plot(h(3),b2(:,1)-t(2),b2(:,4),'r')   %C2
plot(h(3),b3(:,1)-t(3),b3(:,4),'g')   %C3
plot(h(3),b4(:,1)-t(4),b4(:,4),'b')   %C4


% Plot zero-lines
plot(h(1),[b1(1,1)-1,b1(end,1)+1],[0,0],'k--')
plot(h(2),[b1(1,1)-1,b1(end,1)+1],[0,0],'k--')
plot(h(3),[b1(1,1)-1,b1(end,1)+1],[0,0],'k--')


xlim(h(1),[b1(1,1)-.5,b1(end,1)+.5])
xlim(h(2),[b1(1,1)-.5,b1(end,1)+.5])
xlim(h(3),[b1(1,1)-.5,b1(end,1)+.5])

irf_timeaxis(h(3))

clear h fn

end

