%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% usequence.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function D=usequence(n,m,z,p)
% generates a sequence of m vectors in the n-dimensional unit cube 
% such that for each leading subsequence, arbitrary vectos are 
% close to one of the vectors in the subsequence.
% 2*D-1 can be used as an efficient replacement for Halton sequences
% or other low discrepancy sequences.  
%
% input
%
% m      % maximum number of points wanted
% n      % dimesnion of sample points wanted
% z      % determine components type (z_i>0: integer, z_i=0: real)
% p      % first point of the sequence
% det    % generator setting (1: deterministic, 0: non-deterministic)
%
% output
%
% D(:,k) % k-th point of the sequence
% 
function D=usequence(n,m,z,p,det)
if det, % for deterministic case
   rsaved = rand("state");
   % assign rdet 
   rdet   = readrdet;
   rand ("state",rdet);
end;
% for k=1 we pick the zero vector 
D=zeros(n,m); % storage for the sequence to be constructed
% for k=2 we pick the all-one vector
for i=1:n
   if z(i)>0, R(i,:) = randi([-z(i),z(i)],1,m); 
   else, R(i,:)=2*rand(1,m)-1; % random reservoir of m vectors
   end
end
R(:,1)=p; zerM=zeros(1,m);% for later vectorization        
for k=1:m
   if k==1, j=1;
  else 
     % find the reservoir vector with largest minimum distance
     % from the vectors already chosen
     [~,j]=max(u); 
     j=j(1);             % break ties
   end
   D(:,k)=R(:,j); 
   for i=1:n
      if z(i)>0, r(:,i) = randi([-z(i),z(i)],1,1); 
      else, r(:,i)=2*rand-1; % random reservoir of m vectors
      end
   end
   R(:,j)=r; % update the reservoir
   % update minimum squared distance vector 
   onk = ones(1,k);
   u(j)=min(sum((r(:,onk)-D(:,1:k)).^2,1));
   s=sum((R-D(:,k+zerM)).^2,1);
   if k==1, u=s; else, u=min(u,s); end
end
if det, rand("state", rsaved); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%