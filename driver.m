
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% driver.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% driver for reentrant bound-constrained optimization with mixed-integer
% variables
%
% replace in the driver the starting point and the Rosenbrock function 
% by your own problem data.
%
clear
clc
%rng('default')
%%%%%%%%%%%%%%%%%%%%%%
% problem definition %
%%%%%%%%%%%%%%%%%%%%%%
cas=1;  
switch cas
    case 1, st.m=2; st.n=3;   % problem size, [m,n]=size(x)
           x=[-1 -1 -1; -2,-2,-2];  % starting point 
           xx = [x(1,2);x(2,1)];
           indInt=[0 1 1;0 1 0];
    case 2, st.m=1; st.n=2; x=[-1 -1]; xx = x;
        indInt=[0 1];  
end    
dim = st.m*st.n;
% function handle
fun=@(xx)(xx(1)-1)^2+100*(xx(2)-xx(1)^2)^2; % Rosenbrock function
% indices for integer and non-integer 
% lower and upper bound for continuous case
low = x-10*ones(st.m,st.n);
upp = x+10*ones(st.m,st.n);
% lower and upper bounds for mixted-integer Rosenbrock function
low(indInt==1)=0; upp(indInt==0)=100;
% initial point for mixted-integer Rosenbrock function
x(indInt==1)=50; 
% noise will add to the objective function f inside MATRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for stopping test %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tune.cont=[]; tune.int=[]; 
st.prt=2; % print level
st.nfmax=10000; % maximum number of function evaluations
st.secmax=inf; % maximum time in seconds
st.accf=1e-4;  % target accuracy
st.level=0.01; % noise level
if cas==1, st.ftarget=-0.0095; else, st.ftarget=65.1289; end
% MATRS repeatedly calls MATRSstep
[xbest,fbest,info] = MATRS(fun,x,low,upp,st,indInt,tune)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    