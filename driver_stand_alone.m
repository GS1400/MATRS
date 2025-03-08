%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% driver_stand_alone.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% driver for stand-alone bound-constrained optimization with 
% mixed-integer variables
%

%%%%%%%%%%%%%%%%%%%%
% problem definition to be provided by the user
%  fun         % function handle
%  x           % initial point
%  low         % lower bound
%  upp         % upper bound
%  indMixInt   % mixed-integer indices

%%%%%%%%%%%%%%%%%%%%%
% other information to be provided by the user
% st         % structure with start, stop, read and print information
%   .x       % initial point
%   .n       % problem dimension
%   .nfmax   % maximum number of function evaluations
%   .secmax  % maximum time in seconds  
%   .ftarget % stop once a feasible point x with f(x)<=ftarget 
%               is found. ftarget=-inf disables this criteria
%   .prt     % print level
%            %  -1: nothing, 
%            %   0: only improved f 
%            %   1: little, 
%            % >=1: more and more
%
% soldet     % random generator setting
%               1: deterministic (for debugging) 
%               0: non-deterministic 
% tune       % structure including tuning parameters (default if tune=[]) 
%               .cont  % tuning parameters for continuous searches
%               .int   % tuning parameters for integer searches
%               .mint  % tuning parameters for mixed-integer searches
%            % The possible choices for the tuning parameters are 
%            % described in initTune
% solverPath % solver path

clear all
clc
format long

% To solve your problem, adapt all choices in the driver to your own
% problem. In particular, replace in the driver the starting point and
% the FLETCHCR function by your own objective function.
%
% In the examples below you can switch problems by setting solcas=1,2,3,
% types of variable by setting pbcas=1,2,3, and random generator setting
% by setting soldet=0,1
%
% 

solcas=3; pbcas=3; soldet=1; 

% initialize random number generator
if soldet 
   rdet   = readrdet; % assign rdet 
   rand ("state",rdet);
end;

%%%%%%%%%%%%%%%%%%%%%%
% problem definition %
%%%%%%%%%%%%%%%%%%%%%%
n = 2;  % problem dimension

x = 5*ones(n,1); % initial point

low = -100*ones(n,1); % lower bound

upp = 100*ones(n,1); % upper bound

% function handle
switch solcas
    case 1, fun = @(x)FLETCHCR(x,n);
    case 2, fun = @(x)EG2(x,n);
    case 3, fun = @(x)BandedTrigonometric(x,n);
    otherwise
end

% mixed-integer indices
switch pbcas
    case 1 % only continuous variables
     indMixInt = zeros(n,1); 
    case 2 % only integer variable
     indMixInt = ones(n,1); 
    case 3 % mixed-integer variable
     indMixInt = randi([0 1],n,1); % 0 (real), 1 (integer)  
    otherwise
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% structure with start, stop, %
% and print criteria          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stopping test %
st.ftarget = -10^(12); st.nfmax = 100*n; st.secmax = inf;  

st.prt =0; % print level 

st.n = n; % problem dimension

st.x = x; % initial point

st.soldet=soldet;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other information to be % 
% provided by the user    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tuning choice 
tune=[]; % default choice

% solver path
solverPath = '/users/kimiaei/MPCsoftware/MATRS-main';

% path for minq8, used to solve the trust-region subproblem in cMATRS
eval(['addpath ',solverPath,'/minq8']) 
% path for obils, used to solve the trust-region subproblem in iMATRS
eval(['addpath ',solverPath,'/obils'])
% path for test problems
eval(['addpath ',solverPath,'/problems'])

if st.prt>0, rng, end

%%%%%%%%%%%%%%%%%%%%%%%%%
% call the MATRS solver %
%%%%%%%%%%%%%%%%%%%%%%%%%
[xbest,fbest,info] = MATRS(fun,x,low,upp,st,indMixInt,tune);
if st.prt>0 % display results
    disp(['index set indMixInt for mixed-integer variable',...
           ' (1:integer, 0:continuous): '])
    indMixInt',
    disp('xbest found by MATRS:')
    xbest',
    disp(['fbest=f(xbest)=',num2str(fbest)]),
    info
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  