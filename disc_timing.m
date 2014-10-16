function [V] = disc_timing(b1,b2,b3,b4,r1,r2,r3,r4,M,nCluster)
%DISC_TIMING Summary of this function goes here
%   Detailed explanation goes here


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

M
bCut = b(M(1):M(2),:);
bCut1 = b1(M(1):M(2),:);
bCut2 = b2(M(1):M(2),:);
bCut3 = b3(M(1):M(2),:);
bCut4 = b4(M(1):M(2),:);

% Correlation function for velocity calculation
[V,dV] = c_4_v_xcorr([bCut(1,1),bCut(end,1)],b1,b2,b3,b4,r1,r2,r3,r4);

% Calculating time difference from velocity and position
rInd = find_closest_index(bCut1(1,1),r1(:,1));

R1 = r1(rInd,2:4);
R2 = r2(rInd,2:4);
R3 = r3(rInd,2:4);
R4 = r4(rInd,2:4);

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

fArrow = figure;
arrow([0,0,0],nTiming)



plot_disc(bLMN1,bLMN2,bLMN3,bLMN4,tDiff);

titstr = ['s/c ', num2str(nCluster), '     ', 'l2/l3 = ', num2str(l(2)/l(3)),...
    '   ', 'n_{mivar} = ' num2str(v_minvar(3,:)),'    n_{timing} = ' num2str(V/norm(V)) ];
suptitle(titstr);
%v = c_4_v(rC1,rC2,rC3,rC4,tVec);
%v = 1;




end

