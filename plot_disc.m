function [f,h] = plot_disc(b1,b2,b3,b4,t)
%PLOT_DISC Plots time shifted magnetic field in a LMN-frame.
%   PLOT_DISC(b1,b2,b3,b4,t) plots the magnetic field in LMN-frame given
%   magnetic field in LMN-frame b1,...b4 and a time vector with time
%   differences t.
%   [f,h] = PLOT_DISC(b1,b2,b3,b4,t) returns th figure handle f and a 
%   vector of axes handles h=[h(1),h(2),h(3)] for the panels.


%---------------Figure-------------------
f = irf_plot(3,'newfigure');
set(gcf,'PaperUnits','centimeters')
xSize = 15; ySize = 20;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop

h = zeros(1,4);

h(1) = irf_panel('L');
h(2) = irf_panel('M');
h(3) = irf_panel('N');
hold(h(1))
hold(h(2))
hold(h(3))

title(h(1),{' ';' ';' '});
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

%Labels
ylabel(h(1),'B_{L} [nT]','FontSize',16)
ylabel(h(2),'B_{M} [nT]','FontSize',16)
ylabel(h(3),'B_{N} [nT]','FontSize',16)

%No tick labels for the first two panels
set(h(1),'XTickLabel',[])
set(h(2),'XTickLabel',[])

%Set time interval of axes
xlim(h(1),[b1(1,1),b1(end,1)])
xlim(h(2),[b1(1,1),b1(end,1)])
xlim(h(3),[b1(1,1),b1(end,1)])

irf_timeaxis(h(3))


end

