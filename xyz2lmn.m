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

