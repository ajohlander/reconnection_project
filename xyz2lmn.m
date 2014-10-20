function [bLMN] = xyz2lmn(b,v)
%XYZ2LMN Summary of this function goes here
%   Detailed explanation goes here


bLMN = zeros(size(b));

bLMN(:,1) = b(:,1);
bLMN(:,2) = b(:,2:4)*v(1,:)';
bLMN(:,3) = b(:,2:4)*v(2,:)';
bLMN(:,4) = b(:,2:4)*v(3,:)';


end

