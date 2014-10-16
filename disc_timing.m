function [V] = disc_timing(b1,b2,b3,b4,r1,r2,r3,r4,nCluster)
%DISC_TIMING Summary of this function goes here
%   Detailed explanation goes here


coord = 2; %2=max, 3=inter, 4=min
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


%[bLMN1,l1,v1] = irf_minvar(b1);
[~,l,v_minvar] = irf_minvar(b);
bLMN1 = xyz2lmn(b1,v_minvar);
bLMN2 = xyz2lmn(b2,v_minvar);
bLMN3 = xyz2lmn(b3,v_minvar);
bLMN4 = xyz2lmn(b4,v_minvar);


% 
% zInd1 = get_crossing(bLMN1,'mean',coord);
% zInd2 = get_crossing(bLMN2,'mean',coord);
% zInd3 = get_crossing(bLMN3,'mean',coord);
% zInd4 = get_crossing(bLMN4,'mean',coord);
% 
% if(length(zInd1)>1)
%     disp('Warning: more than one zero crossing')
%     zInd1 = zInd1(1);
% end
% if(length(zInd2)>1)
%     disp('Warning: more than one zero crossing')
%     zInd2 = zInd2(1);
% end
% if(length(zInd3)>1)
%     disp('Warning: more than one zero crossing')
%     zInd3 = zInd3(1);
% end
% if(length(zInd4)>1)
%     disp('Warning: more than one zero crossing')
%     zInd4 = zInd4(1);
% end
% 
% 
% tD1 = bLMN1(zInd1,1);
% tD2 = bLMN2(zInd2,1);
% tD3 = bLMN3(zInd3,1);
% tD4 = bLMN4(zInd4,1);
% tVec = [tD1,tD2,tD3,tD4];
% tDiff = tVec-tD1;
% %tDiff = [tD1-tD1,tD2-tD1,tD3-tD1,tD4-tD1];



[V,dV] = c_4_v_xcorr([b(1,1),b(end,1)],b1,b2,b3,b4,r1,r2,r3,r4);

R1 = r1(1,2:4);
R2 = r2(1,2:4);
R3 = r3(1,2:4);
R4 = r4(1,2:4);
tDiff = [0,(R2-R1)*V'/(V*V'),(R3-R1)*V'/(V*V'),(R4-R1)*V'/(V*V')];

fn = plot_disc(bLMN1,bLMN2,bLMN3,bLMN4,tDiff);

titstr = ['s/c ', num2str(nCluster), '     ', 'l2/l3 = ', num2str(l(2)/l(3)),...
    '   ', 'n_{mivar} = ' num2str(v_minvar(3,:)),'    n_{timing} = ' num2str(V/norm(V)) ];
suptitle(titstr);
%v = c_4_v(rC1,rC2,rC3,rC4,tVec);
%v = 1;




end

