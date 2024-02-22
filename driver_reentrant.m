%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% driver_reentrant.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% driver for reentrant bound-constrained optimization 
% with mixed-integer variables
%
% This deriver shows the use of the reentrant implementation of 
% the MATRS solver. MATRSstep takes a history of the previous feasible
% points and their function values and returns only a new feasible point
% for evaluation. Thus  the execution of a solver on a single processing 
% system is stopped when a function value is requested and resumes with
% a new call once the function value is available. This provides maximal
% flexibility for its use in tunning software 
%

%%%%%%%%%%%%%%%%%%%%
% problem definition to be provided by the user
%
%  fun         % function handle
%  x           % initial point
%  low         % lower bound
%  upp         % upper bound
%  indMixInt   % mixed-integer indices

%%%%%%%%%%%%%%%%%%%%%
% other information to be provided by the user
%
%  soldet     % random generator setting
%                1: deterministic (for debugging) 
%                0: non-deterministic 
%  tune       % structure including tuning parameters (default if tune=[])
%                .cont  % tuning parameters for continuous searches
%                .int   % tuning parameters for integer searches
%                .mint  % tuning parameters for mixed-integer searches
%             % The possible choices for the tuning parameters are 
%             % described in initTune
%  solverPath % solver path
%  nfmax      % maximum number of function evaluations
%  prt        % print level
%                -1: nothing, 
%                 0: only improved f 
%                 1: little, 
%               >=1: more and more
%

clear all
clc
format long

% To solve your problem, adapt all choices in the driver to your own
% problem. In particular, replace in the driver the starting point and
% the Rosenbrock function by your own objective function.
%
% You can switch random generator by setting soldet=0,1
%
soldet=0; 

%%%%%%%%%%%%%%%%%%%%%%
% problem definition %
%%%%%%%%%%%%%%%%%%%%%%

x = [-1 -1];  % initial point 

low = -100*ones(2,1); % lower bound

upp = 100*ones(2,1); % upper bound

indMixInt=[0 1]; % mixed-integer indices 

nfmax=1000;  % stopping test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other information provided by the user %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print level
prt =0;    

% tuning choice
tune=[]; % default choice

% solver path %
solverPath = ...
          input('e.g.: /users/kimiaei/MPCsoftware/MATRS-main >> ','s');

% path for minq8, used to solve the trust-region subproblem in cMATRS
eval(['addpath ',solverPath,'/minq8']) 
% path for obils, used to solve the trust-region subproblem in iMATRS
eval(['addpath ',solverPath,'/obils']) 

% initialize random number generator
if soldet, 
   % assign rdet for deterministic case
   rdet   = readrdet;
   rand ("state",rdet);
end;
if prt>0, rng, end

%%%%%%%%%%%%%%%%%%%%%%%%%
% call the MATRS solver %
%%%%%%%%%%%%%%%%%%%%%%%%%
% MATRSinit in the paper 
MATRSstep(tune,[],prt,indMixInt,low,upp,nfmax); 

% loop with reentrant call 
for nf=1:nfmax
    f=(x(1)-1)^2+100*(x(2)-x(1)^2)^2; % compute function value
    x= MATRSstep(x(:),f,prt);
    if prt>=0,
       if nf==1, fb=f;
       elseif f<fb
          fb=f; 
          disp(['f = ',num2str(fb),' improved at nf = ',num2str(nf)])
        end
    end
end

% generate output by calling MATRSread in the paper
[xbest,fbest,infostep]= MATRSstep 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%