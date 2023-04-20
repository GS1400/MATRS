
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fitGrad.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fitGrad estimates the gradient by fitting
%
function [g]=fitGrad(n,Y,F,xb,fb,maxY,tune)
Y         = adjustY(Y,tune); % adjust Y
Y         = Y';
xb        = xb';
K         = min(size(Y,1)-1,maxY);
distance  = sum((Y-ones(size(Y,1),1)*xb).^2,2);
[~,ind]   = sort(distance);
ind       = ind(2:K+1);
S         = Y(ind,:) -ones(K,1)*xb;
R         = triu(qr(S,0));
R         = R(1:n,:);
warning off
L = inv(R)';
warning on
sc     = sum((S*L').^2,2).^(3/2);
sc     = adjustVec(sc,tune);
b      = (F(ind)-fb)'./sc;
b      = adjustVec(b,tune);
A      = [Y(ind,:) - ones(K,1)*xb];
A      = A./(sc*ones(1,size(A,2)));
warning off
y = A\b;
warning on
y  = adjustVec(y,tune);
g = y(1:n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%