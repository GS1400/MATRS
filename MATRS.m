
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,info] = MATRS(fun,x,low,upp,st,indInt,tune);
%
% solves the bound-constrained noisy derivative-free optimization problem 
%    min f(x) 
% with mixed-integer variables
%
% fun          % function handle for objective function
% x            % starting point (must be specified)
% low          % lower bound on x
% upp          % upper bound on x
% st           % structure with stop and print criteria
%              %   (indefinite run if no stopping criterion is given)
%  .secmax     %   stop if sec>=secmax (default: inf)
%  .nfmax      %   stop if nf>=nfmax   (default: inf)
%  .qs         %   stop if qs<=accf    (default: 1e-4)
%  .ftarget    %   function value accepted as optimal (default: 0)
%  .prt        %   printlevel (default: -1)
%              %   -2: nothing, -1: litte, >=0: more and more
% indInt       % indices for type of variables
%              %  0: continuous, 1: integer
% tune         % optional structure specifying tuning parameters
%              %   for details see initTune.m
%
% x            % best point found 
% f            % function value at best point found 
% info         % performance information for MADFOrun
%  .finit      %   initial function value
%  .ftarget    %   target function value (to compute qf)
%  .qf         %   (ftarget-f)/(finit-f)
%  .initTime   %   inital cputime
%  .done       %   done with the search?
%  .acc        %   stop when qf<=acc
%  .secmax     %   stop if sec>=secmax 
%  .nfmax      %   stop if nf>=nfmax 
%  .finit      %   the initial f
%  .prt        %   printlevel 
%              %     -1: nothing, 0: litte, >=1: more and more
% 
function [xbest,fbest,info] = MATRS(fun,x,low,upp,st,indInt,tune);
persistent m n dim finfo
if nargin<7
   message = 'MATRS needs fun, x, low, upp, st, indInt, tune as input';
   disp(message)  
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% initial checks %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(st,'prt'), prt=st.prt; else, prt=-1; end
if prt>=0,
  disp(' ')
  disp('==============================================================')
  disp('start MATRS')
  disp('==============================================================')
end;
% check function handle
if isempty(fun)
  message = 'MATRS needs the function handle fun to be defined';
  disp(message)  
  return
elseif ~isa(fun,'function_handle')
  message = 'fun should be a function handle';
  disp(message)
  return
end;
% starting point
if isempty(x)
  message = 'starting point must be defined';
  disp(message) 
  return      
elseif ~isa(x,'numeric')
  message = 'x should be a numeric vector'; 
  disp(message ) 
  return       
end
m=st.m; n=st.n; dim=m*n;
[mx,nx]=size(x);
if mx~=m && nx~=n,
  sizex=[mx,nx], sizeNeeded=[m,n] 
  error('dimension mismatch');
end
if dim==0
  % no variables
  xbest=zeros(0);fbest=inf;
  finfo.error='no variables';
  return;
end
if isempty(low)
  message = 'low must be defined';
  disp(message) 
  return      
elseif ~isa(low,'numeric')
  message = 'low should be a numeric vector'; 
  disp(message ) 
  return       
end
if isempty(upp)
  message = 'upp must be defined';
  disp(message) 
  return      
elseif ~isa(upp,'numeric')
  message = 'upp should be a numeric vector'; 
  disp(message ) 
  return       
end
if ~all(low<=upp) 
  % feasible domain empty
  xbest=NaN;fbest=NaN;
  finfo.error='feasible domain empty';
  return;
end
if length(upp)~=dim
  nupp=dim,nupp=length(upp)
  finfo.error='dimension mismatch'; 
end
if length(low)~=dim
  nlow=dim,nlow=length(low)
  finfo.error='dimension mismatch'; 
end
% add info on stopping criteria
if isfield(st,'secmax'), info.secmax=st.secmax;
else, info.secmax=inf;
end;
if isfield(st,'nfmax'), info.nfmax=st.nfmax;
else, info.nfmax=inf;
end;
if isfield(st,'ftarget'), info.ftarget=st.ftarget;
else, info.ftarget=0;
end;
if isfield(st,'accf'), info.accf=st.accf;
else, info.acc=1e-4;
end;
% noise level
if isfield(st,'level'), level=st.level;
else, level=0;
end;
info.prt = prt; % print level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% initialize solver environment %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MATRSstep(tune,fun,prt,indInt,low,upp,info.nfmax)
initTime=cputime;
%%%%%%%%%%%%%%%%%%%
%%%% main loop %%%%
%%%%%%%%%%%%%%%%%%%
nf=0;  
while 1
  % get function value
  f  = fun(x);
  %%%%%%%%%%%%%%%%%%
  % add noise to f %
  %%%%%%%%%%%%%%%%%%
  epsilon = (2*rand-1)*level; % uniform
  fs = f+epsilon; % absolute
  nf = nf+1; % count the number of function evaluations
  % get new point for function evaluation
  x    = MATRSstep(x(:),fs,prt);
  x    = reshape(x,m,n); 
  % restore original format of x
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % checking stopping test %
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  sec       = (cputime-initTime);
  info.done = (sec>info.secmax)|(nf>=info.nfmax);
  if nf>1
     qs        = (f-info.ftarget)/(finit-info.ftarget);
     qs        = max(qs,0);
     info.qs   = qs;
     info.done = (info.done|info.qs<=info.accf);
     if qs<qsb,
         if prt>=0
            disp(['relative function value qs improved at nf=',...
                num2str(nf),' to f=',num2str(qs)]) 
         end
         qsb=qs; 
     end
  else
      finit = f; qsb=inf;
  end
  info.sec  = sec;
  if info.done, break; end
end
%%%%%%%%%%%%%%%%%%%%%
%%%% return info %%%%
%%%%%%%%%%%%%%%%%%%%%
[xbest,fbest,infostep]= MATRSstep;
% update info
info.niter_iTRS          = infostep.number_iTRS; 
info.niter_cTRS          = infostep.number_cTRS; 
info.niter_mixed_integer = infostep.number_mixed_integer;
info.nf                  = infostep.nf;
info.qs_ConvergeQuality  = qsb;
info.times_in_second     = sec;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% solution status  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if info.qs_ConvergeQuality<=st.accf
  info.status_of_converge = 'accuracy reached';
elseif nf>=st.nfmax
   info.status_of_converge = 'nfmax reached';
elseif sec>=st.secmax
  info.status_of_converge = 'secmax reached';
elseif fbest<=-1e+12
  info.status_of_converge ='function -1e+40 reached';
else
  info.status_of_converge = 'unknown';
end;
if prt>=0,
  disp('==============================================================')
  disp('end MATRS')
  disp('==============================================================')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

