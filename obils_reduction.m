%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% obils_reduction.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [R,y,l,u,p] = obils_reduction(B,y,l,u) reduces the general overdetermined 
% box-constrained integer least squares problem to an upper triangular one 
% by the QR factorization with column permutations 
% Q'*B*P = [R; 0]. The orthogonal matrix Q is not produced.  
%
% Inputs:
%    B - m-by-n real matrix with full column rank
%    y - m-dimensional real vector
%    l - n-dimensional integer vector, lower bound 
%    u - n-dimensional integer vector, upper bound
%
% Outputs:
%    R - n-by-n real nonsingular upper triangular matrix
%    y - n-dimensional vector transformed from the input y, y:=(Q'*y)(1:n)
%    l - permuted input lower bound l, i.e., l := P'*l 
%    u - permuted input upper bound u, i.e., u := P'*u
%    p - n-dimensional permutation vector representing P

% Main Reference: 
% S. Breen and X.-W. Chang. Column Reordering for 
% Box-Constrained Integer Least Squares Problems, 
% Proceedings of IEEE GLOBECOM 2011, 6 pages.

% Authors: Xiao-Wen Chang, www.cs.mcgill.ca/~chang
%          Xiangyu Ren
% Copyright (c) 2015-2018. Scientific Computing Lab, McGill University.
% Last revision: December 2018

function [R,y,l,u,p] = obils_reduction(B,y,l,u)


[~,n] = size(B);

% Transform B and y by the QR factorization 
U = qr([B,y]);
R = triu(U(1:n,1:n));
y = U(1:n,n+1);
R0 = R;
y0 = y;

% Permutation vector 
p = 1:n;

% Inverse transpose of R
G = inv(R)';

% Determine the column permutatons 
for k = n : -1 : 2
    maxDist = -1;
    
    % Determine the k-th column
    for i = 1:k
        alpha = y(i:k)' * G(i:k,i);
        x_i = max(min(round(alpha),u(i)),l(i));
        if (alpha < l(i) || alpha > u(i) || alpha == x_i)
            dist = 1 + abs(alpha - x_i); 
        else 
            dist = 1 - abs(alpha - x_i);
        end
        dist_i = dist / norm(G(i:k,i));
        if dist_i > maxDist
            maxDist = dist_i;
            j = i;
            x_j = x_i;
        end
    end

    % Perform permutations
    p(j:k) = p([j+1:k,j]);
    l(j:k) = l([j+1:k,j]);
    u(j:k) = u([j+1:k,j]);

    % Update y, R and G for the new dimension-reduced problem
    y(1:k-1) = y(1:k-1) - R(1:k-1,j) * x_j;
    R(:,j) = [];
    G(:,j) = [];

    for t = j : k - 1
        % Triangularize R and G by Givens rotation        
        [W, R([t,t+1],t)] = planerot(R([t,t+1],t));
        R([t,t+1],t+1:k-1) = W * R([t,t+1],t+1:k-1);
        G([t,t+1],1:t) = W * G([t,t+1],1:t);
        % Apply the Givens rotation W to y
        y(t:t+1) = W * y(t:t+1);
    end
end

% Reorder the columns of R0 according to p  
R0 = R0(:,p);

% Transform R0 and y0 by the QR factorization
U = qr([R0,y0]);
R = triu(U(:,1:n));
y = U(:,n+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

