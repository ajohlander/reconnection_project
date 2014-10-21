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

ar = zeros(1,4);
firstArrow = arrow(R./Re,(R+nTiming*arrLength)./Re,'FaceColor',arrColor(5,:),'EdgeColor',arrColor(5,:));
ar(1) = arrow(R1./Re,(R1+nMVA1*arrLength)./Re,'FaceColor',arrColor(1,:),'EdgeColor',arrColor(1,:));
ar(2) = arrow(R2./Re,(R2+nMVA2*arrLength)./Re,'FaceColor',arrColor(2,:),'EdgeColor',arrColor(2,:));
ar(3) = arrow(R3./Re,(R3+nMVA3*arrLength)./Re,'FaceColor',arrColor(3,:),'EdgeColor',arrColor(3,:));
ar(4) = arrow(R4./Re,(R4+nMVA4*arrLength)./Re,'FaceColor',arrColor(4,:),'EdgeColor',arrColor(4,:));


delete(firstArrow) %First arrow becomes wierd
arrow(R./Re,(R+nTiming*arrLength)./Re,'FaceColor',arrColor(5,:),'EdgeColor',arrColor(5,:))

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
LMNstr3 = ['n_{timing} = [' num2str(V/norm(V)),']'];


title(h(1),{LMNstr1;[LMNstr2,LMNstr3]})

%suptitle(legstr);
%v = c_4_v(rC1,rC2,rC3,rC4,tVec);
%v = 1;




end

