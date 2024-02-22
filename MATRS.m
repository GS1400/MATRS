%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,info] = MATRS(fun,x,low,upp,st,indMixInt,tune);
%
% solves the bound-constrained derivative-free optimization problem 
%    min f(x), s.t. low <= x <= upp 
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
%  .ftarget    %   target function value (default: -inf)
%                  it is explained in driver
%  .prt        %   printlevel (default: -1)
%              %    -1: nothing, 
%              %     0: only improved f 
%              %     1: litte, 
%              %   >=1: more and more
% indMixInt    % indices for type of variables
%              %  0: continuous, 1: integer
% tune         % optional structure specifying tuning parameters
%              %   for details see initTune.m
%
% x            % best point found 
% f            % function value at best point found 
% info         % performance information for MATRS
%  .ftarget    %   target function value 
%  .initTime   %   inital cputime
%  .done       %   done with the search?
%  .secmax     %   stop if sec>=secmax 
%  .nfmax      %   stop if nf>=nfmax 
%  .prt        %   print level
%              %     -1: nothing, 
%              %      0: only improved f 
%              %      1: little, 
%              %    >=1: more and more
% 
function [xbest,fbest,info] = MATRS(fun,x,low,upp,st,indMixInt,tune);
if nargin<7
   message = 'MATRS needs fun, x, low, upp, st, indMixInt, tune as input';
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
info.error='';
% check function handle
if isempty(fun)
  message = 'MATRS needs the function handle fun to be defined';
  info.error=message;
  xbest=NaN;fbest=NaN;
  return
elseif ~isa(fun,'function_handle')
  message = 'fun should be a function handle';
  info.error=message;
  xbest=NaN;fbest=NaN;
  return
end;
% starting point
if isempty(x)
  message = 'starting point must be defined';
  info.error=message;
  xbest=NaN;fbest=NaN;
  return      
elseif ~isa(x,'numeric')
  message   = 'x should be a numeric vector'; 
  info.error= message;
  xbest=NaN;fbest=NaN;
return       
end
dim=length(x);
if isempty(low)
  message = 'low must be defined';
  info.error=message;
  xbest=NaN;fbest=NaN;
return      
elseif ~isa(low,'numeric')
  message = 'low should be a numeric vector'; 
  info.error=message;
  xbest=NaN;fbest=NaN;
  return       
end
if isempty(upp)
  message = 'upp must be defined';
  info.error=message; 
  xbest=NaN;fbest=NaN;
  return      
elseif ~isa(upp,'numeric')
  message = 'upp should be a numeric vector'; 
  info.error=message; 
  xbest=NaN;fbest=NaN;
  return       
end
if ~all(low<=upp) 
  % feasible domain empty
  xbest=NaN;fbest=NaN;
  info.error='feasible domain empty';
  xbest=NaN;fbest=NaN;
  return;
end
if length(upp)~=dim
  nupp=dim,nupp=length(upp)
  info.error='dimension mismatch'; 
  return;
end
if length(low)~=dim
  nlow=dim,nlow=length(low)
  info.error='dimension mismatch';
  return;
end
% add stopping criteria to info
if isfield(st,'secmax'), info.secmax=st.secmax;
else, info.secmax=inf;
end;
if isfield(st,'nfmax'), info.nfmax=st.nfmax;
else, info.nfmax=inf;
end;
if isfield(st,'ftarget'), info.ftarget=st.ftarget;
else, info.ftarget=-inf;
end;
info.prt = prt; % print level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% initialize solver environment %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initTime=cputime;
info.initTime=initTime;
MATRSstep(tune,[],prt,indMixInt,low,upp,info); % MATRSinit of the paper
%%%%%%%%%%%%%%%%%%%
%%%% main loop %%%%
%%%%%%%%%%%%%%%%%%%
nf=0; Ii = find(indMixInt==1);
while 1
  % get function value
  f  = fun(x);
  nf = nf+1; % count the number of function evaluations
  x  = MATRSstep(x(:),f,prt);
  if nf==1, fb=f; end
  if ~isempty(x)
      if any(floor(x(Ii))-x(Ii)~=0)
         input(['bug in MATRSstep; there is a least one real component',...
                ' in the integer components of x'])
         x(Ii) = round(x(Ii));
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%
      % checking stopping test %
      %%%%%%%%%%%%%%%%%%%%%%%%%%
      sec       = (cputime-initTime);
      info.done = (sec>info.secmax)|(nf>=info.nfmax);
      if nf>1, info.done = (info.done|f<=info.ftarget);end
      info.sec  = sec;
      if f<fb, 
          fb=f; 
          if prt>=0
             disp(['f = ',num2str(fb),...
                   ' improved at nf = ',num2str(nf)]) 
          end
      end
  else
      %%%%%%%%%%%%%%%%%%%%%%%%%%
      % checking stopping test %
      %%%%%%%%%%%%%%%%%%%%%%%%%%
      sec       = (cputime-initTime);
      info.done = 1;
      info.sec  = sec;
  end
  if info.done, break; end
end
%%%%%%%%%%%%%%%%%%%%%%%
%%%% generate info %%%%
%%%%%%%%%%%%%%%%%%%%%%%
[xbest,fbest,infostep]= MATRSstep; % MATRSread of the paper
% update info
info.niter_iTRS          = infostep.number_iTRS;
info.niter_cTRS          = infostep.number_cTRS; 
info.niter_miMATRS       = infostep.number_miMATRS;
info.nf                  = infostep.nf;
info.times_in_second     = sec;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% solution status  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if fbest<=st.ftarget, info.status_of_converge = 'accuracy reached';
elseif nf>=st.nfmax, info.status_of_converge = 'nfmax reached';
elseif sec>=st.secmax, info.status_of_converge = 'secmax reached';
elseif fbest<=-1e+12, info.status_of_converge ='f is unbounded below';
else, info.status_of_converge = 'unknown';
end;
if prt>=0,
  disp('==============================================================')
  disp('end MATRS')
  disp('==============================================================')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%