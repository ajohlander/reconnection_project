function out = c_4_v_timing_mva(x1,x2,x3,x4,R,column)
%C_4_V_TIMING_MVA Performs timing and minimum variance analysis on four s/c
%field data with a graphical user interface.
%   C_4_V_TIMING_MVA(b1,b2,b3,b4,R,column) interactive discontinuity
%   analysis analyzer on magnetic field b1,...b4 with position R using
%   column number 'column'. R has the form R.R1,...R.R4.
%   C_4_V_TIMING_MVA('B?',R) uses B1, B2, B3 and B4 from workspace
%   C_4_V_TIMING_MVA('B?',R, column) also uses column for selecting the
%   component of the magnetic field. 2=x, 3=y, 4=z. If
%   column is not passed as an argument, the default value is 2.
%
%   The user is asked to define an interval in which the discontinuity
%   analysis should be made. Two new figures open when this is done.
%   In the first window, the normal vectors obtained from MVA and timing is
%   shown in 3D space. Normal vectors that do not fulfill l2/l3>5 is shown
%   as dotted arrows. The dotted arrows are not used in determining the
%   maximum angle.
%   In the second figure, the magnetic field is plotted in the
%   LMN-system, where L=maximum, M=intermiediate and N = minimum.
%
%   A prompt asks the user to save the variables: velocity of discontinuity V,
%   uncertainty in velocity dV, all 5 normal vectors in n, where
%   n.nTiming is the normal vector from timing and n.n1,...n4 is normal
%   vectors obtained from MVA, all vectors from MVA v, and eigenvalues from
%   MVA. The variables do not appear until the main figure is closed.



column = 0;

% Getting parameters
if(nargin == 2 || nargin == 3 ) % Input ('B?', 'R?',...) or ('B?', R,...)
    
    if(ischar(x1) && strfind(x1,'?')) % 'B?'
        var_strB = x1;
        evalin('base',['if ~exist(''' irf_ssub(var_strB,1) '''), c_load(''' var_strB ''');end' ]);
        c_eval('b?=evalin(''base'',irf_ssub(var_strB,?));');
    end
    
    if(ischar(x2))
        if(strfind(x2,'?'))
            var_strR = x2;
            evalin('base',['if ~exist(''' irf_ssub(var_strR,1) '''), c_load(''' var_strR ''');end' ]);
            c_eval('R.R?=evalin(''base'',irf_ssub(var_strR,?));');
        end
    else % R
        R = x2;
    end
    
    if(nargin == 3)
        column = x3;
    end
        
elseif(nargin ==5 || nargin ==6)    %Input (b1,b2,b3,b4,R/'R?',...)
    b1 = x1;
    b2 = x2;
    b3 = x3;
    b4 = x4;
    if(ischar(R) && strfind(R,'?'))
        var_strR = R;
        evalin('base',['if ~exist(''' irf_ssub(var_strR,1) '''), c_load(''' var_strR ''');end' ]);
        c_eval('R.R?=evalin(''base'',irf_ssub(var_strR,?));');
    end % R
end

if(column == 0) % Not in input
    irf.log('notice','Using second column');
    column = 2; % x is the default column.
end

tint = [b1(1,1),b1(end,1)];


%---------------GUI-------------------
fGUI = irf_plot(4,'newfigure');
set(gcf,'PaperUnits','centimeters')
xSize = 15; ySize = 15;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop

h = zeros(1,4);
h(1) = irf_panel('C1');
h(2) = irf_panel('C2');
h(3) = irf_panel('C3');
h(4) = irf_panel('C4');
hold(h(1))
hold(h(2))
hold(h(3))
hold(h(4))

%Bcolumn
pb = zeros(1,4);
pb(1) = plot(h(1),b1(:,1),b1(:,column));
pb(2) = plot(h(2),b2(:,1),b2(:,column));
pb(3) = plot(h(3),b3(:,1),b3(:,column));
pb(4) = plot(h(4),b4(:,1),b4(:,column));

%line at 0
fLine = plot(h(1),tint,[0,0]);
set(fLine,'Color',[.5 .5 .5])
fLine = plot(h(2),tint,[0,0]);
set(fLine,'Color',[.5 .5 .5])
fLine = plot(h(3),tint,[0,0]);
set(fLine,'Color',[.5 .5 .5])
fLine = plot(h(4),tint,[0,0]);
set(fLine,'Color',[.5 .5 .5])

%Meta
xlim(h(1),tint)
xlim(h(2),tint)
xlim(h(3),tint)
xlim(h(4),tint)


ystr = '';

switch column
    case 2
        ystr = 'B_{x} [nT]';
    case 3
        ystr = 'B_{y} [nT]';
    case 4
        ystr = 'B_{z} [nT]';
end

ylabel(h(1),ystr,'FontSize',16)
ylabel(h(2),ystr,'FontSize',16)
ylabel(h(3),ystr,'FontSize',16)
ylabel(h(4),ystr,'FontSize',16)

irf_timeaxis(h)

% Loop for the GUI
while true
    % Makes the GUI-window the current figure
    nFig = get(fGUI,'Parent');
    figure(nFig{1})
    
    clear x % For good measure
    try
        [x,~] = ginput(2);
    catch % If user closes window, the program closes
        disp('Figure closed. Any saved variables can be found in the Workspace');
        return
    end
    
    % Remove any existing colored area
    if (exist('a','var')==1)
        delete(a)
    end
    % The s/c that was clicked is the reference s/c
    clickedAx = gca;
    nCluster = find(h == clickedAx);
    
    b = 0;
    switch nCluster
        case 1
            b = b1;
        case 2
            b = b2;
        case 3
            b = b3;
        case 4
            b = b4;
    end
    
    % Finds the time interval from user data
    clickInd1 = find_closest_index(x(1),b(:,1));
    clickInd2 = find_closest_index(x(2),b(:,1));
    
    M = [clickInd1, clickInd2];
    
    if(M(1)>M(2))
        M = fliplr(M);
    end
    
    absB = sqrt(sum(b(:,2:4)'.^2));
    maxB = max(absB);
    
    %Draws a colored area that marks the time interval chosen
    a = zeros(1,4);
    for i = 1:4
        a(i) = area(h(i),[b(M(1,1),1) b(M(1,2),1)], [maxB maxB],-maxB, 'FaceColor', [0.5,1,0.5]);
        set(a(i),'EdgeColor','none')
        uistack(pb(i),'top')
    end
    
    
    % Calls for the timing function
    [V,dV,n,v,l] = disc_timing(b1,b2,b3,b4,R,M,nCluster);
    
    % Asks the user to save the variables to the Workspace.
    % The function does not return any variables due to the loop structure.
    checkLabels = {'Save velocity to variable named:' ...
        'Save dV to variable named:' 'Save normal vectors to variable named:'...
        'Save MVA matrices to variable named:' 'Save MVA eigenvalues to variable named:'};
    varNames = {'V','dV','n','v','l'};
    items = {V,dV,n,v,l};
    export2wsdlg(checkLabels,varNames,items,...
        'Save outputs to Workspace');
    
    out = V;
end
end



function [V,dV,n,v,L] = disc_timing(b1,b2,b3,b4,R,M,nCluster)
%DISC_TIMING Performes timing and MVA for magnetic field data.
%   [V,dV,n,v,L] = DISC_TIMING(b1,b2,b3,b4,R,M,nCluster) 

switch nCluster
    case 1
        b = b1;
    case 2
        b = b3;
    case 3
        b = b3;
    case 4
        b = b4;
end

bCut = b(M(1):M(2),:);
bCut1 = b1(M(1):M(2),:);
bCut2 = b2(M(1):M(2),:);
bCut3 = b3(M(1):M(2),:);
bCut4 = b4(M(1):M(2),:);

%Possible bug: c_4_v_xcorr interpolates R, might be too few points.

% Correlation function for velocity calculation
[V,dV] = c_4_v_xcorr([bCut(1,1),bCut(end,1)],b1,b2,b3,b4,R.R1,R.R2,R.R3,R.R4);

% Calculating time difference from velocity and position
rInd = find_closest_index(bCut1(1,1),R.R1(:,1));

R1 = R.R1(rInd,2:4);
R2 = R.R2(rInd,2:4);
R3 = R.R3(rInd,2:4);
R4 = R.R4(rInd,2:4);

switch nCluster
    case 1
        tDiff = [(R1-R1)*V'/(V*V'),(R2-R1)*V'/(V*V'),(R3-R1)*V'/(V*V'),(R4-R1)*V'/(V*V')];
    case 2
        tDiff = [(R1-R2)*V'/(V*V'),(R2-R2)*V'/(V*V'),(R3-R2)*V'/(V*V'),(R4-R2)*V'/(V*V')];
    case 3
        tDiff = [(R1-R3)*V'/(V*V'),(R2-R3)*V'/(V*V'),(R3-R3)*V'/(V*V'),(R4-R3)*V'/(V*V')];
    case 4
        tDiff = [(R1-R4)*V'/(V*V'),(R2-R4)*V'/(V*V'),(R3-R4)*V'/(V*V'),(R4-R4)*V'/(V*V')];
end

%
%-------------MVA-------------

%Time-shifting the four intervals
M1 = [find_closest_index(bCut1(1,1)+tDiff(1),b1(:,1)),find_closest_index(bCut1(end,1)+tDiff(1),b1(:,1))];
M2 = [find_closest_index(bCut2(1,1)+tDiff(2),b2(:,1)),find_closest_index(bCut2(end,1)+tDiff(2),b2(:,1))];
M3 = [find_closest_index(bCut3(1,1)+tDiff(3),b3(:,1)),find_closest_index(bCut3(end,1)+tDiff(3),b3(:,1))];
M4 = [find_closest_index(bCut4(1,1)+tDiff(4),b4(:,1)),find_closest_index(bCut4(end,1)+tDiff(4),b4(:,1))];

% bShift is the magnetic field centered around the discontinuity for each
% s/c
bShift = b(M(1):M(2),:);
bShift1 = b1(M1(1):M1(2),:);
bShift2 = b2(M2(1):M2(2),:);
bShift3 = b3(M3(1):M3(2),:);
bShift4 = b4(M4(1):M4(2),:);

% Minimum variance analysis. Possibly do nested?
[~,l,v_minvar] = irf_minvar(bShift);
[~,l1,v_minvar1] = irf_minvar(bShift1);
[~,l2,v_minvar2] = irf_minvar(bShift2);
[~,l3,v_minvar3] = irf_minvar(bShift3);
[~,l4,v_minvar4] = irf_minvar(bShift4);

%LMN coordinate system of the clicked s/c
bLMN1 = xyz2lmn(bShift1,v_minvar);
bLMN2 = xyz2lmn(bShift2,v_minvar);
bLMN3 = xyz2lmn(bShift3,v_minvar);
bLMN4 = xyz2lmn(bShift4,v_minvar);

nTiming = V/norm(V);
nMVA1 = v_minvar1(3,:);
nMVA2 = v_minvar2(3,:);
nMVA3 = v_minvar3(3,:);
nMVA4 = v_minvar4(3,:);


theta1 = acosd(dot(nTiming,nMVA1));
theta2 = acosd(dot(nTiming,nMVA2));
theta3 = acosd(dot(nTiming,nMVA3));
theta4 = acosd(dot(nTiming,nMVA4));

if theta1>90
    nMVA1 = -nMVA1;
end
if theta2>90
    nMVA2 = -nMVA2;
end
if theta3>90
    nMVA3 = -nMVA3;
end
if theta4>90
    nMVA4 = -nMVA4;
end

% Outputs
% normal vector
n = [];
n.nTiming = nTiming;
n.n1 = nMVA1;
n.n2 = nMVA2;
n.n3 = nMVA3;
n.n4 = nMVA4;
% eigenvalues
L = [];
L.l1 = l1;
L.l2 = l2;
L.l3 = l3;
L.l4 = l4;
% all vectors
v = [];
v.v1 = v_minvar1;
v.v2 = v_minvar2;
v.v3 = v_minvar3;
v.v4 = v_minvar4;


% ----------Plotting the arrows-----------
R = mean([R1;R2;R3;R4]); %Center of the tetrahedron
dist = [norm(R-R1),norm(R-R2),norm(R-R3),norm(R-R4)];
arrLength = max(dist);
Re = 6371; %Earth radius

arrColor = [[0 0 0];[1 0 0];[0 1 0];[0 0 1];[0 .5 .5]];
fArrow = irf_plot(1,'newfigure');
hold(fArrow);
axis equal

ar = zeros(1,5);

%-----OLD AND NOW BROKEN CODE-------

% firstArrow = arrow(R./Re,(R+nTiming*arrLength)./Re,'FaceColor',arrColor(5,:),'EdgeColor',arrColor(5,:));
% ar(1) = arrow(R1./Re,(R1+nMVA1*arrLength)./Re,'FaceColor',arrColor(1,:),'EdgeColor',arrColor(1,:));
% ar(2) = arrow(R2./Re,(R2+nMVA2*arrLength)./Re,'FaceColor',arrColor(2,:),'EdgeColor',arrColor(2,:));
% ar(3) = arrow(R3./Re,(R3+nMVA3*arrLength)./Re,'FaceColor',arrColor(3,:),'EdgeColor',arrColor(3,:));
% ar(4) = arrow(R4./Re,(R4+nMVA4*arrLength)./Re,'FaceColor',arrColor(4,:),'EdgeColor',arrColor(4,:));


% delete(firstArrow) %First arrow becomes wierd
% arrow(R./Re,(R+nTiming*arrLength)./Re,'FaceColor',arrColor(5,:),'EdgeColor',arrColor(5,:))

%--------NEW CODE FOR MATLAB 2014b AND LATER--------------
ar(1) = line([R1(1), R1(1)+nMVA1(1)*arrLength]/Re,...
    [R1(2), R1(2)+nMVA1(2)*arrLength]/Re,...
    [R1(3), R1(3)+nMVA1(3)*arrLength]/Re);
ar(2) = line([R2(1), R2(1)+nMVA2(1)*arrLength]/Re,...
    [R2(2), R2(2)+nMVA2(2)*arrLength]/Re,...
    [R2(3), R2(3)+nMVA2(3)*arrLength]/Re);
ar(3) = line([R3(1), R3(1)+nMVA3(1)*arrLength]/Re,...
    [R3(2), R3(2)+nMVA3(2)*arrLength]/Re,...
    [R3(3), R3(3)+nMVA3(3)*arrLength]/Re);
ar(4) = line([R4(1), R4(1)+nMVA4(1)*arrLength]/Re,...
    [R4(2), R4(2)+nMVA4(2)*arrLength]/Re,...
    [R4(3), R4(3)+nMVA4(3)*arrLength]/Re);

ar(5) = line([R(1), R(1)+nTiming(1)*arrLength]/Re,...
    [R(2), R(2)+nTiming(2)*arrLength]/Re,...
    [R(3), R(3)+nTiming(3)*arrLength]/Re);


set(ar(1),'Color',arrColor(1,:))
set(ar(2),'Color',arrColor(2,:))
set(ar(3),'Color',arrColor(3,:))
set(ar(4),'Color',arrColor(4,:))
set(ar(5),'Color',arrColor(5,:))


%ar(1).Color = arrColor(1,:);


% check if l2/l3 is ok for each s/c.
eigLim = 5; %Minimum value for l2/l3
mvaOK = zeros(1,5);
mvaOK([l1(2)/l1(3),l2(2)/l2(3),l3(2)/l3(3),l4(2)/l4(3)]>eigLim) = 1;
mvaOK(5) = 1; %Timing is always OK
for i = 1:4
    if mvaOK(i)
        set(ar(i),'LineWidth',1)
    else
        set(ar(i),'LineStyle',':')
    end
end


%plot the positions of the s/c
plot3(R1(1)/Re,R1(2)/Re,R1(3)/Re,'ko')
plot3(R2(1)/Re,R2(2)/Re,R2(3)/Re,'ro')
plot3(R3(1)/Re,R3(2)/Re,R3(3)/Re,'go')
plot3(R4(1)/Re,R4(2)/Re,R4(3)/Re,'bo')
plot3(R(1)/Re,R(2)/Re,R(3)/Re,'Color',arrColor(5,:),'Marker','x')

view(-37.5,15) %3D-view

xlabel('X  [R_{E}]','FontSize',16)
ylabel('Y  [R_{E}]','FontSize',16)
zlabel('Z  [R_{E}]','FontSize',16)


set(gca,'ColorOrder',arrColor)
irf_legend(gca,{'C1','C2','C3','C4','Timing',},[0.00, 1.00]);

%Maximum angle between the normal vectors
thetaMax = 0;
nMatrix = [nMVA1;nMVA2;nMVA3;nMVA4;nTiming];

for i = 1:5
    for j = 1:5
        if (i~=j && mvaOK(i) && mvaOK(j))
            theta = acosd(dot(nMatrix(i,:),nMatrix(j,:)));
            if theta>thetaMax
                thetaMax = theta;
            end
        end
    end
end
arrStr1 = ['V = [', num2str(V),'] \pm ','[', num2str(dV),']'];
arrStr2 =  ['Maximum normal deviation = ', num2str(thetaMax), '^{o}'];
title({arrStr1;arrStr2})

%eigStr = [num2str(l1(2)/l1(3)),num2str(l2(2)/l2(3)),' ',...
%    num2str(l3(2)/l3(3)),' ',num2str(l4(2)/l4(3))];

irf_legend({'l_2/l_3'},[0.90 1.00])
%set(gca,'ColorOrder',arrColor)
irf_legend({num2str(l1(2)/l1(3)),num2str(l2(2)/l2(3)),...
    num2str(l3(2)/l3(3)),num2str(l4(2)/l4(3))},[1.00, 0.95])

%---------------------------------------------------------------

%Plots the magnetic field in the LMN coordinate system
[~,h] = plot_disc(bLMN1,bLMN2,bLMN3,bLMN4,tDiff);

LMNstr1 = ['s/c ', num2str(nCluster), '     ', 'l2/l3 = ', num2str(l(2)/l(3))];
LMNstr2 = ['n_{minvar} = [', num2str(v_minvar(3,:)),'] '];
LMNstr3 = ['n_{timing} = [', num2str(V/norm(V)),']'];


title(h(1),{LMNstr1;[LMNstr2,LMNstr3]})

end


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


function [index] = find_closest_index(num,vec)
%FIND_CLOSEST_INDEX Finds relevant index of a vector
%   index = ANJO.FINDCLOSESTINDEX(num, vec) returns the index of the vector
%   vec where vec is the closest to the number num.

index = find(num>=vec, 1, 'last' );

if index == length(vec)
    return;
elseif(num-vec(index) >= vec(index+1)-num)
    index = index+1;
end

if isempty(index)
    index = 1;
end
end


function [bLMN] = xyz2lmn(b,v)
%XYZ2LMN Coordinate transformation from GSE or GSM to a LMN-frame.
%   [bLMN] = XYZ2LMN(b,v) returns the magnetic field in a minimum variance
%   frame bLMN given magnetic field data 'b' and a matrix 'v' containing the
%   eigenvectors from minimum variance analysis.

bLMN = zeros(size(b));

bLMN(:,1) = b(:,1);
bLMN(:,2) = b(:,2:4)*v(1,:)';
bLMN(:,3) = b(:,2:4)*v(2,:)';
bLMN(:,4) = b(:,2:4)*v(3,:)';

end

