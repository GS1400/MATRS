%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% minq8sep.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,y,nit,ier,acc,itref] = minq8sep(c,d,A,b,eq,maxit,prt)
% solves the problem min c'*x+0.5*x'*diag(d)*x s.t. A*x>= b, 
% (A*x)(eq)=b(eq) by solving the dual problem with MINQ8 and making one 
% step of iterative refinement if 
% abs(r)<=nnz(A)*eps*(abs(A)*abs(x)+abs(b)) (with componentwise 
% inequalities) is violated, where r is the vector of violations of the 
% constraints
% 
% Input:
% c      vector of length n
% d      vector of length n with positive entries
% A      m times n matrix
% b      vector of length b
% eq     index set of at most length n containing indices in 1:m
%        indicates the equality constraints
% maxit  limit on the number of iterations in each call to MINQ8
%        (default m)
% prt    = 0  no printing (default)
%        ~= 0 prints intermediate function values in MINQ8
%
% Output
% x      optimal primal solution  (if vector of length n)
%        no feasible point exists (if empty)
% y      (vector of length m) optimal dual solution (if x nonempty)
%                      certificate of infeasibility (if x empty)
% nit    sum of numbers of iterations of the two calls to MINQ8
% ier    error flag
%        = 0 accuracy satisfactory
%        = 1 feasible set probably empty
% acc    maximum of the absolute values of the constraint violations
%        (if x nonemtpy)
%        norm(A'*y,inf) if x empty 
% itref  = 0 if no iterative refinement was needed
%        = 1 if one step of iterative refinement was made
%
% Calls minq8.m and its subprograms
%
function [x,y,nit,ier,acc,itref] = minq8sep(c,d,A,b,eq,maxit,prt)
if min(d)<=0, error('diagonal must be positive'), end
small = 1.e-8;
tol = eps; % tolerance for the stopping criterion of MINQ8
[m,n] = size(A);
data.gam = 0;
data.c = -b;
data.A = A';
data.b = c;
data.D = 1./d;
ineq = 1:m;
ineq(eq) = [];  % index set of the inequality constraints
yu = Inf*ones(m,1);
yl = zeros(m,1);yl(eq) = -yu(eq);
ys = zeros(m,1);
if nargin < 6
  maxit = m;
end
if nargin < 7
  prt = 0;
end
[y,f,g,nit,ier] = minq8(data,yl,yu,ys,maxit,tol,prt);
if norm(ier)>0  % check whether ier is a certificate of inreasibility
  z = A'*ier;
  zeta = b'*ier;
  if norm(z)<small && zeta > 0 && length(find(ier>=0)) == length(ier)
  % ier is a certificate of infeasibility
    x = [];
    y = ier;
    ier = 1;
    itref = 0;
    acc = max(abs(z));
    return 
  end
end
x=(A'*y-c)./d;
% check for accuracy
res=A*x-b;ressmall=nnz(A)*eps*(abs(A)*abs(x)+abs(b))
res(ineq)=min(res(ineq),0);
acc = max(abs(res));
if isempty(find(abs(res)>ressmall))
  % accuracy satisfactory
  ier=0;
  itref=0;
  return
end

% one step of iterative refinement
itref=1;
data.c=res;
data.b = zeros(n,1);
[dy,f,g,nit0,ier]=minq8(data,yl-y,yu,zeros(m,1),maxit,tol,prt);
x=x+(A'*dy)./d;
y=y+dy;
nit = nit+nit0;

% check for accuracy
res=A*x-b;ressmall=nnz(A)*eps*(abs(A)*abs(x)+abs(b));
res(ineq)=min(res(ineq),0);
acc = max(abs(res));
if isempty(find(abs(res)>sqrt(nnz(A))*ressmall))
  % accuracy satisfactory
  ier=0;
else
  % feasible set probably empty
  ier=1;
end;






