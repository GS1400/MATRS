%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% obils_search.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% zhat = obils_search(R,y,l,u,beta) produces the optimal solution to 
% the upper triangular box-constrained integer least squares problem
% min_{z}||y-Rz|| s.t. z in [l, u] by a search algorithm.

% Inputs:
%    R - n by n real nonsingular upper triangular matrix
%    y - n-dimensional real vector
%    l - n-dimensional integer vector, lower bound 
%    u - n-dimensional integer vector, upper bound    
%
% Outputs:
%    zhat - n-dimensional integer vector (in double precision). 

% Subfunctions: init, update (included in this file)

% Main Reference: 
% X.-W. Chang and Q. Han. Solving Box-Constrained Integer Least 
% Squares Problems, IEEE Transactions on Wireless Communications,  
% 7 (2008), pp. 277-287.


% Authors: Xiao-Wen Chang, www.cs.mcgill.ca/~chang
%          Xiangyu Ren, Jiqqun Shen
% Copyright (c) 2015. Scientific Computing Lab, McGill University.
% April 2015. Last revision: December 2015


function zhat = obils_search(R,y,l,u)


% ------------------------------------------------------------------
% --------  Initialization  ----------------------------------------
% ------------------------------------------------------------------

n = size(R,1);

% Current point
z = zeros(n,1);

% c(k)=(y(k)-R(k,k+1:n)*z(k+1:n))/R(k,k)
c = zeros(n,1);

% d(k): left or right search direction at level k   
d = zeros(n,1);

% lflag(k) = 1 if the lower bound is reached at level k
lflag = zeros(size(l));
% uflag(k) = 1 if the upper bound is reached at level k
uflag = lflag;

% Partial squared residual norm for z
% prsd(k) = (norm(y(k+1:n)-R(k+1:n,k+1:n)*z(k+1:n)))^2
prsd = zeros(n,1);

% Store some quantities for efficiently calculating c
% S(k,n) = y(k),
% S(k,j-1) = y(k) - R(k,j:n)*z(j:n) = S(k,j) - R(k,j)*z(j), j=k+1:n
S = zeros(n,n);
S(:,n) = y;

% path(k): record information for updating S(k,k:path(k)-1) 
path = n*ones(n,1); 

% The level at which search starts to move up to a higher level
ulevel = 0;  

% Initial search radius
beta = Inf; 

% ------------------------------------------------------------------
% --------  Search process  ----------------------------------------
% ------------------------------------------------------------------

k = n;

% Find the initial integer in [l(n), u(n)] 
c(k) = S(k,k) / R(k,k);
[z(k),d(k),lflag(k),uflag(k)] = init(c(k),l(k),u(k));
gamma = R(k,k) * (c(k) - z(k));

% dflag for down or up search direction 
dflag = 1; % Intend to move down to a lower level

ii=0;
while 1
    ii=ii+1;
    if ii>=100, zhat = z; break; end
    if dflag == 1
       % Temporary partial squared residual norm at level k
       newprsd = prsd(k) + gamma * gamma;
       
       if newprsd < beta  % Inside the ellispoid
           if k ~= 1 % Move to level k-1
               % Update path  
               if ulevel ~= 0 
                   path(ulevel:k-1) = k;
                   for j = ulevel-1 : -1 : 1
                       if path(j) < k
                           path(j) = k;
                       else
                           break;  % Note path(1:j-1) >= path(j)
                       end
                   end
               end
               
               % Update S
               k = k - 1;
               for j = path(k):-1:k+1
                   S(k,j-1) = S(k,j) - R(k,j)*z(j);
               end
               
               % Update the partial squared residual norm
               prsd(k) = newprsd;
               
               % Find the initial integer in [l(k), u(k)]
               c(k) = S(k,k) / R(k,k);
               [z(k),d(k),lflag(k),uflag(k)] = init(c(k),l(k),u(k));
               gamma = R(k,k) * (c(k) - z(k));
               
           else    % A valid point is found, update search radius
               zhat = z;
               beta = newprsd;
           end          
           ulevel = 0;           
       else           
           dflag = 0;  % Will move back to a higher level           
       end       
    else
       % Outside the ellispoid
       if k == n % The optimal solution has been found, terminate
           break;
       else
       % Move back to level k+1
           if ulevel == 0
               ulevel = k;
           end
           k = k + 1;
           if lflag(k) ~= 1 || uflag(k) ~= 1
               % Find a new integer at level k  
               [z(k),d(k),lflag(k),uflag(k)] = ...
                      update(z(k),d(k),lflag(k),uflag(k),l(k),u(k));
               gamma = R(k,k) * (c(k) - z(k)); 
               dflag = 1;
           end
       end
    end
end


% ------------------------------------------------------------------
% --------  Subfunctions  ------------------------------------------
% ------------------------------------------------------------------

function [z_k,d_k,lflag_k,uflag_k] = init(c_k,l_k,u_k)
%
% Find the initial integer and the search direction at level _k 

z_k = round(c_k);
if z_k <= l_k
    z_k = l_k;
    lflag_k = 1;  % The lower bound is reached
    uflag_k = 0;
    d_k = 1;
elseif z_k >= u_k
    z_k = u_k;
    uflag_k = 1;  % The upper bound is reached
    lflag_k = 0;
    d_k = -1;
else
    lflag_k = 0; 
    uflag_k = 0;
    if c_k > z_k
        d_k = 1;
    else
        d_k = -1;
    end
end


function [z_k,d_k,lflag_k,uflag_k] = update(z_k,d_k,lflag_k,uflag_k,l_k,u_k)
%
% Find a new integer at level k and record it if it hits a boundary.

z_k = z_k + d_k; 
if z_k == l_k
    lflag_k = 1;
    d_k = -d_k + 1;
elseif z_k == u_k
    uflag_k = 1;
    d_k = -d_k - 1;
elseif lflag_k == 1
    d_k = 1;
elseif uflag_k == 1
    d_k = -1;
else
    if d_k > 0
        d_k = -d_k - 1;
    else
        d_k = -d_k + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


