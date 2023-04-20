%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% iusequence.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function D=iusequence(N,d);
% generates a sequence of N vectors in the d-dimensional unit cube 
% such that for each leading subsequence, arbitrary vectos are 
% close to one of the vectors in the subsequence.
%
% MA     % data struture for Matrix adaptation parameters
% N      % number of sample points wanted 
% d      % dimesnion of sample points wanted
% m      % maximum number of points wanted

% D(:,k) % k-th point of the sequence
% 
function D=iusequence(MA,ilambda,d,m)
% for k=1 we pick the zero vector 
D=zeros(d,ilambda);  % storage for the sequence to be constructed
m=max(m,ilambda);
R = randi([-MA.int.dir,MA.int.dir],d,m);
if isfield(MA.int,'pinit')
    R(:,1)=MA.int.pinit;
end       
zerm=zeros(1,m);     % for later vectorization        
for k=1:ilambda
  % find the reservoir vector with largest minimum distance
  % from the vectors already chosen
   % pick a vector from the reservoir
  if k==1
     j=1;
  else 
     % find the reservoir vector with largest minimum distance
     % from the vectors already chosen
     [~,j]=max(u); 
     j=j(1);             % break ties
  end
  D(:,k)=R(:,j);  
  ii=0;
  while ii<10
      ii=ii+1;
      r= randi([-MA.int.dir,MA.int.dir],d,1);
      if norm(r)~=0, break;end
  end
  R(:,j)=r;  % update the reservoir
  % update minimum squared distance vector 
  onk = ones(1,k); u(j)=min(sum((r(:,onk)-D(:,1:k)).^2,1));
  s=sum((R-D(:,k+zerm)).^2,1);
  if k==1, u=s; else, u=min(u,s); end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




