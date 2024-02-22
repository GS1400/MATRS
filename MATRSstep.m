%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRSstep.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performs one step of MATRS
% 
%
% function MATRSstep(itune,ctune,fun,prt,ni,nc); % initialize solver
%
% itune      % structure contains those tuning integer parameters that  
%            % are explicitly specified; the others take default values
%            %   (for choices and defaults see initTune.m)
%
% ctune      % structure contains those tuning continuous parameters that  
%            % are explicitly specified; the others take default values
%            %   (for choices and defaults see initTune.m)
%
% mitune     % structure contains those tuning mixed-integer parameters 
%            % that are explicitly specified; the others take default
%            % values (for choices and defaults see initTune.m)
%
% fun          % function handle (empty if used in reentrant mode)
% prt          % print level
% ni           % number of integer vaiables
% nc           % number continous vaiables
%
%
% function x=MATRSstep(x,f,prt); % suggest new point
%
% x            % point evaluated (input) or to be evaluated (output)
% f            % function value at input x
% prt          % change print level (optional)
%
% function [xbest,fbest,infos]=MATRSstep();  % read results
%
% xbest        % current best point
% fbest        % function value found at xbest
% infos        % performance information for MATRS
%  .finit      %   initial function value
%  .ftarget    %   target function value, explained in driver
%  .initTime   %   inital cputime
%  .done       %   done with the search?
% 
function [x,f,infos] = MATRSstep(x,f,p,indMixInt,low,upp,budget)

persistent MA state prt itune ctune mitune info ni nc dim 

persistent xbest fbest  xfinal ffinal ordering nfmax secmax initTime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MATRSread of the paper %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0 % return info 
  if ffinal<fbest, f=ffinal; x=xfinal;
  else, f=fbest; x=xbest;
  end
  infos.nf             = info.nf;
  infos.number_iTRS    = info.niTR; 
  infos.number_cTRS    = info.ncTR; 
  infos.number_miMATRS = info.number_miMATRS;
  infos.number_succ_cTRS  = MA.cont.succ_cTRS;
  infos.number_succ_iTRS  = MA.int.succ_iTRS;
  infos.state = state;
  return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MATRSinit of the paper %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout==0 % initialize solver environment
  info = []; dim = size(indMixInt,2)*size(indMixInt,1);
  indMixInt=indMixInt(:); low=low(:); upp=upp(:);
  % problem dimension
  % subspace dimension of continuous variables
  info.cI = find(indMixInt==0); 
  % subspace dimension of integer variables
  info.iI = setdiff(1:dim,info.cI);
  nc = length(info.cI); ni = length(info.iI); 
  if nc>0 % lower and upper bound on continuous variables
    info.clow = low(info.cI); info.cupp = upp(info.cI);
  end
  if ni>0 % lower and upper bound on continuous variables
     info.ilow = low(info.iI); info.iupp = upp(info.iI);
  end
  if isempty(x)
     xc = []; yc = []; zc=[];
  else
     xc = x.cont; yc = x.int; zc=x.mint;
  end
  % tuning parameters
  [ctune,itune,mitune] = initTune(xc,yc,zc,nc,ni);
  % budgets for stopping tests
  if isfield(budget,'nfmax'), nfmax = budget.nfmax; 
  else, nfmax=budget;
  end
  if isfield(budget,'secmax'), secmax = budget.secmax; 
  else, secmax=inf;
  end
  if isfield(budget,'initTime'),initTime = budget.initTime; end
  info.nc=nc; info.ni=ni;  prt = p; infos=info; state=1; 
  return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% initialization %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if state==1 
    xbest = x; fbest=f; info.nf=1; ordering = 'ascend';
    if ni>0 % initialization for cMATRS
        MA.int          = [];
        MA.int.sigma    = itune.sigmai;
        MA.int.ar       = itune.ainit;
        MA.int.Parent.y = x;
        MA.int.s        = zeros(ni, 1);
        MA.int.M        = itune.M;
        MA.int.echi     = sqrt(ni)*(1-1/ni/4-1/ni/ni/21);
        MA.int.a        = eye(itune.ilambda,1); 
        MA.int.kappa    = 1;
        MA.int.D        = ones(ni,itune.ilambda);
        MA.int.good     = 1;
        MA.int.succ_iTRS = 0;
        MA.int.unsucc_iTRS =0;
    end
    if nc>0 % initialization for iMATRS
        MA.cont          = [];
        MA.cont.sigma    = ctune.sigmac;
        MA.cont.ar       = ctune.ainit;
        MA.cont.echi     = sqrt(nc)*(1-1/nc/4-1/nc/nc/21);
        MA.cont.Parent.y = x;
        MA.cont.s        = zeros(nc, 1);
        MA.cont.M        = ctune.M;
        MA.cont.kappa    = 1;
        MA.cont.D        = ones(nc,ctune.clambda); 
        MA.cont.a        = ones(ctune.clambda,1);
        MA.cont.good     = 1;
        MA.cont.succ_cTRS = 0;
        MA.cont.unsucc_cTRS =0;
    end
    % create the first list for saving all evaluated points and their f
    info.XF=Inf*ones(nfmax,dim+1);
    info.XF(1,1:dim)=x'; info.XF(1,dim+1)=fbest;
    % create the second list, used to compute combination directions
    info.X = xbest; info.F = fbest; info.dim=dim;
    MA.int.iit=0; MA.cont.cit=0; 
    info.ncTR =0; % number of calls to cTRS
    info.niTR =0; % number of calls to iTRS
    info.number_miMATRS=0; % number of calls to miMATRS
    MA.int.succ_iTRS =0; MA.int.unsucc_iTRS=0;
    MA.cont.succ_cTRS =0; MA.cont.unsucc_cTRS=0;
    if nc>0 
       state=2; % go to cMATRS
    elseif nc==0 && ni>0
        state=30; % go to iMATRS
    end
    if ni>0&& nc>0, MA.mixed.ae = mitune.ainit;end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATRS alternates calls to cMATRS, iMATRS, and miMATRS in this order %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save new point and its function value
xfinal = x; ffinal = f; info.f=f; 
while 1 % main loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% cMATRS: contniuous MATRS %%%%%%%%%%%%%%%%%%%
    %%%% cMATRS includes cMutation, selection, cRecom, and cTRS %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nc>0  % start of cMATRS
       if state==2 % initialization for cMutation
            MA.cont.cit = MA.cont.cit+1;  
           if prt>=1
               disp('========================================')
               fprintf([num2str(MA.cont.cit),'th cMATRS\n'])
               disp('========================================')
           end
           if prt>=1 
             fprintf([num2str(MA.cont.cit),'th cMutation is done \n'])
           end
           if ~MA.cont.good, MA.cont.OffspringPop0=MA.cont.MutInfoList;
           else, MA.cont.OffspringPop0 =[];
           end
           info.xbest0  = xbest; 
           MA.cont.goodm  = 0; MA.cont.MutInfoList = []; MA.cont.cm = 1;
           state  = 3;
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % cMutation: continuous mutation phase %
       % cMutation calls cLSS                 %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       while state<14 && state>=3  % cMutation is performing`
            if state==3 
               % compute continuous dustribution and mutation directions
               MA.cont.MutInfo.sd = MA.cont.D(:,MA.cont.cm); 
               MA.cont.MutInfo.p  = MA.cont.M*MA.cont.MutInfo.sd;
               MA.cont.cxbest     = xbest(info.cI);  
               % check feasibility
               % find largest real step size 
               % find  mutation step size
               info.ls ='c'; info.step='m';
               [MA,info]=getalpha(MA,info,ctune);
               % check whether or not cMutation ends
               if ~MA.cont.feasible
                    state=8; % go to try opposite direction
                    info.ftrial=[]; info.ftrial0=[];
               else % cMutation does not stop
                    info.ytrial = xbest;
                    % compute the first trial point and project it 
                    % into [low upp]
                    info.statep='c'; info.beta=MA.cont.alp0;
                    info.p=MA.cont.MutInfo.p;
                    [info]=projectStep(xbest,info);
                    state=4;
                    x=info.ytrial;
                    MA.cont.SSext=info.beta;
                    return; % return MATRSstep to compute f at x
               end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % cLSS: continuous line search %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if state==4 % f at x was computed; initialize extrapolation
                % save new point and its function value in the history
                info=UpdatePoints(info,mitune);
                MA.cont.extrapDone=0;
                MA.cont.SSext=[MA.cont.SSext info.beta];
                % check whether the first trial point is 
                % the first trial point of extrapolation
                if info.ftrial<fbest 
                    info.enf = info.nf-1;
                    MA.cont.extrapDone=1;
                    MA.cont.alpha = MA.cont.alp0; 
                    % calculate the second trial point
                    % along the mutation direction
                    info.ytrial  = xbest;
                    % compute the step size
                    info.beta  = ...
                          min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                    % compute trial point and project it into [low upp]
                    info.statep='c';info.p=MA.cont.MutInfo.p;
                    [info]=projectStep(xbest,info);
                    state=5;
                    x=info.ytrial;
                    return; % return MATRSstep to compute f at x
                else % no decrease in f 
                    state=8; % go to try opposite direction
                end
            end
            if state==5 
                % save new point and its function value in the history
                info=UpdatePoints(info,mitune);
                MA.cont.SSext=[MA.cont.SSext info.beta];
                state=6; % go to do extrapolation 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % extrapolation along    %
            % the mutation direction %
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            while state==6||state==7 % extrapolation is doing
                if state==6
                    % increase step size
                    if (MA.cont.alpha<MA.cont.alpmax) && ...
                                                (info.ftrial < fbest)
                        % compute new feasible point
                        MA.cont.alpha = ... 
                        min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                        % chekc whether alpmax reaches or not
                        if(MA.cont.alpha < MA.cont.alpmax)
                           info.ytrial = xbest; 
                           % compute new step size
                           info.beta   = ...
                           min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                           % compute trial point and project it 
                           % into [low upp]
                           info.statep='c';info.p=MA.cont.MutInfo.p;
                           [info]=projectStep(xbest,info);
                           state=7;
                           x=info.ytrial;
                           return; % return MATRSstep to compute f at x
                        end
                    end
                    % go to update information
                    state=13; 
                    break
                end
               if state==7 
                   % f at x was computed; save x and f in the history
                   info=UpdatePoints(info,mitune);
                   state=6; % go back to continue extrapolation
                   MA.cont.SSext=[MA.cont.SSext info.beta];
               end
            end % end of extrapolation
            if state==8 % opposite direction is tried
                MA.cont.MutInfo.p   = -MA.cont.MutInfo.p; 
                MA.cont.cxbest      = xbest(info.cI);       
                info.ls             = 'c'; info.step='m';
                [MA,info]           = getalpha(MA,info,ctune);
                if ~MA.cont.feasible % check termination of cMutation
                    state=13; % go to update information
                else % cMutation is peroming
                     % compute new trial point
                     if~isempty(info.ftrial),info.ftrial0=info.ftrial;end
                     info.ytrial = xbest;
                     % compute trial point and project it into [low upp]
                     info.statep='c'; info.beta=MA.cont.alp0;
                     info.p=MA.cont.MutInfo.p;
                     [info]=projectStep(xbest,info);
                     state=9;
                     x=info.ytrial;
                     MA.cont.SSext=info.beta; 
                    return; % return MATRSstep to compute f at x
                end
            end
            % check whether extrapolation can be done or not
            if state==9 
                  % f at x was computed; save x and f in the history
                  info=UpdatePoints(info,mitune);
                  MA.cont.SSext=[MA.cont.SSext info.beta];
                  % check whether the first trial point is 
                  % the first trial point of extrapolation
                 if info.ftrial<fbest
                    MA.cont.extrapDone=1; 
                    info.enf = info.nf-1;
                    % initialize step size and best point
                    MA.cont.alpha = MA.cont.alp0;
                    % calculate a new feasible trial point
                    info.ytrial  = xbest;
                    info.beta    = ... 
                           min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                    % compute trial point and project it into [low upp]
                    info.statep='c'; info.p=MA.cont.MutInfo.p;
                    [info]=projectStep(xbest,info);
                    state=10;
                    x=info.ytrial;
                    return; % return MATRSstep to compute f at x
                 else % no decrease in f along opposite direction
                     MA.cont.MutInfo.sd = -MA.cont.MutInfo.sd;
                     state=13; % go to update information
                      if ~isempty(info.ftrial0)
                         if info.ftrial0<info.ftrial
                             % function value at mutation point
                             info.ftrial=info.ftrial0; 
                         end
                      end
                 end
            end
            if state==10 % save x and f in the history
                info=UpdatePoints(info,mitune);
                state=11; % go to do extrapolation
                 MA.cont.SSext=[MA.cont.SSext info.beta];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extrapoltion along opposite %
            % mutation direction          %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            while state==11||state==12 % loop for extrapolation
                if state==11
                    % expansion step (increase stepsize)
                    if (MA.cont.alpha<MA.cont.alpmax) && ...
                                              (info.ftrial < fbest)
                        % step size calulation and best point updating
                        MA.cont.alpha = ...
                          min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                        % next point to be tested
                        if(MA.cont.alpha < MA.cont.alpmax)
                           info.ytrial = xbest;  
                           % compute new step size
                           info.beta = min(MA.cont.alpmax,...
                                         ctune.cnu*MA.cont.alpha);
                           % compute trial point and project it
                           % into [low upp]
                           info.statep='c'; info.p=MA.cont.MutInfo.p;
                           [info]=projectStep(xbest,info);
                           state=12;
                           x=info.ytrial;
                           return; % return MATRSstep to compute f at x
                        end
                    end
                    % no decrease in f along both p and -p
                    state=13; % go to update information
                    break
                end
               if state==12 % save x and f in the history
                   info=UpdatePoints(info,mitune);
                   state=11; % go back to do more extrapolation
                   MA.cont.SSext=[MA.cont.SSext info.beta];
               end
            end % end of ectrapolation along opposite mutation direction
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update information after     %
            % termination of extrapolation %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if state==13 
                 if MA.cont.extrapDone 
                      % extrapolation was done
                      % extrapolation with amost two trial points
                      % one of which with the lowest inexact f
                      % is accepted as a new best point
                      X = info.XF(info.enf:info.nf-1,1:dim)';
                      F = info.XF(info.enf:info.nf-1,dim+1);
                      [ftrial,ib]=min(F);
                      ib=ib(1); ytrial = X(:,ib);
                      if ftrial < fbest, fbest=ftrial; xbest=ytrial; end
                      if prt>=1 
                        fprintf(['cMutation: fbest=',...
                                                 num2str(fbest),'\n'])
                      end
                      MA.cont.a(MA.cont.cm) = MA.cont.SSext(ib);
                      MA.cont.goodm=1;
                      if prt>=2
                         fprintf('best point was updated\n')
                      end
                 else % extrapolation cannot find a reduction of f
                      % update list of extrapolation step sizes
                      if MA.cont.feasible
                           MA.cont.a(MA.cont.cm) = info.beta;
                      else
                          MA.cont.a(MA.cont.cm) = MA.cont.sigma;
                      end
                      MA.cont.MutInfo.f   = info.ftrial;
                      MA.cont.MutInfoList = [MA.cont.MutInfoList ...
                                              MA.cont.MutInfo];
                 end
                 % check whether cMutation ends or not
                  MA.cont.cm=MA.cont.cm+1;
                  if MA.cont.cm>size(MA.cont.D,2)
                      state=14; 
                      break;% end of cMutation
                  else % go back to continue cMutation
                      state=3; 
                  end
            end
       end % end of cMutation           
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % computation of the set of    %
       % distribution directions by   %
       % usequnce or random generator %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if state==14  % change the set D of distribution directions
          pinit       = xbest(info.cI)-info.xbest0(info.cI);
          normpinit   = norm(pinit,inf);
          changDc = (normpinit == 0);
          if ~changDc % use usequence to generate D
               MA.cont.pinit = pinit;
               MA.cont.pinit = ctune.sc*MA.cont.pinit/normpinit;
               MA.cont.D     =  MA.cont.kappa*usequence(nc,...
                               ctune.clambda,zeros(nc,1),MA.cont.pinit,0);
               if MA.cont.kappa<ctune.ckmax
                   MA.cont.kappa=MA.cont.kappa+1;
               else, MA.cont.kappa=1;
               end
          else % use the normal distribution to generate D
              MA.cont.D  = randn(nc,ctune.clambda);
              MA.cont.kappa = 1;    
          end
          %%%%%%%%%%%%%%%%%%%
          % selection phase %
          %%%%%%%%%%%%%%%%%%%
          MA.cont.lambda = length(MA.cont.MutInfoList);
          if ~MA.cont.good && MA.cont.lambda<6
             MA.cont.MutInfoList = ...
                        [MA.cont.OffspringPop0,MA.cont.MutInfoList];
          end
          clambda = length(MA.cont.MutInfoList);
          MA.cont.good = MA.cont.goodm;
          MA.cont.goodr=0; MA.cont.okTR=~MA.cont.goodm;
          MA.cont.mu    = ceil(clambda/2);
          if clambda>=1 % initialization before selection
             cwi_raw       = log(clambda/2 + 0.5) - log((1:MA.cont.mu));
             MA.cont.wi    = cwi_raw/sum(cwi_raw);
             cmu_eff       = max(1+eps,1/sum(MA.cont.wi .^2));
             MA.cont.c_s   = min(1.999,(cmu_eff + 2)/(nc + cmu_eff + 5));
             MA.cont.c_1   = 2/((nc+1.3)^2 + cmu_eff);
             MA.cont.c_mu  = .... 
                    min( 1-MA.cont.c_1, 2*(cmu_eff - 2 + 1/cmu_eff)/...
                    ((nc + 2)^2 + cmu_eff));
             MA.cont.ds    = 1 + MA.cont.c_s + ...
                                 2*max(0,sqrt((cmu_eff-1)/(nc+1))-1);
             MA.cont.sqrt_s = sqrt(MA.cont.c_s*(2-MA.cont.c_s)*cmu_eff);
             if prt>=1
                 fprintf([num2str(MA.cont.cit),...
                                   'th selection is done \n'])
             end
             MA.cont.permut = sortData(MA.cont.MutInfoList, ordering);
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % cRecom: continuous recombination phase %
          % cRecom calls cLSS                      %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if MA.cont.mu>=2
            % compute recombination distribution direction MA.cont.sumz 
            %         recombination mutation directions MA.cont.sumd 
            MA.cont.sumz = zeros(nc,1); MA.cont.sumd = zeros(nc,1); 
            for m = 1:MA.cont.mu 
                MA.cont.sumz=MA.cont.sumz+MA.cont.wi(m)*...
                          MA.cont.MutInfoList(MA.cont.permut(m)).sd;
                MA.cont.sumd=MA.cont.sumd+MA.cont.wi(m)*...
                          MA.cont.MutInfoList(MA.cont.permut(m)).p;
             end
             MA.cont.okRec = (norm(MA.cont.sumd)~=0 && ~MA.cont.goodm);
             state=15; % go to do 
           else
             state=26; % skip cRecom
          end
       end
       if state==15 
         if MA.cont.okRec % start of cRecom
              if prt>=2
                 fprintf([num2str(MA.cont.cit),...
                          'th cRecom is done\n'])
              end
             MA.cont.cxbest = xbest(info.cI);  nc = length(info.cI);
            okRecombonation=1;  MA.cont.goodr=0; cr=0; beta=1;
            % find alpmax and alp0
            while okRecombonation % scale recombination direction at most
                                  % nscale times so that some feasible
                                  % trial points can be found
               cr=cr+1; 
               if cr>1,beta = 2*rand(nc,1)-1; beta=beta/norm(beta);end
               MA.cont.RecomInfo.p = beta.*MA.cont.sumd;
               okRecombonation = (norm(MA.cont.RecomInfo.p)~=0);
               if okRecombonation % find alpmax and alp0
                   info.ls='c'; info.step = 'r';
                   [MA,info]=getalpha(MA,info,ctune);
                   if MA.cont.alpmax>0
                       state=16;
                       break;
                   end
               end
               if cr>=ctune.nscale
                   state=26; % skip cRecom
                   break;
               end
            end
         else
             state=26; % skip cRecom
         end
       end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      % start of cLSS by cRecom %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      if state==16 % compute the trial recombination point
        % Build first point for starting linesearch
        info.ytrial  = xbest;
        % project trial point into [low upp]
        info.statep='c'; info.beta=MA.cont.alp0;
        info.p = MA.cont.RecomInfo.p;
        [info] = projectStep(xbest,info);
        state=17; x=info.ytrial;
        MA.cont.SSext=info.beta;
        return; % return MATRSstep to compute f at x 
      end
      if state==17 % save x and f in the history
        info=UpdatePoints(info,mitune);
        state=18;
        MA.cont.SSext=[MA.cont.SSext info.beta];
      end
      if state==18
             % check whether the first trial point can be 
             % the first trial point of extrapolation
             MA.cont.extrapDone=0;
            if info.ftrial<fbest
                info.enf = info.nf-1;
                MA.cont.extrapDone=1;
                % initialize MA.cont.alpha and best point
                MA.cont.alpha = MA.cont.alp0; 
                % calculate trial point
                info.ytrial = xbest;
                info.beta=min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                % project trial point into [low upp]
                info.statep='c'; info.p=MA.cont.RecomInfo.p;
                [info]=projectStep(xbest,info);
                state=19;
                x=info.ytrial;
                return; % return MATRSstep to compute f at x
            else
                state=199; % go to do opposite direction
            end
      end
      if state==19 % save x and f in the history
         state=20; info=UpdatePoints(info,mitune);
         MA.cont.SSext=[MA.cont.SSext info.beta];
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % start of extrapolation along     %
      % recombination mutation direction %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       while state==20||state==21 % extrapolation is doing
            if state==20
                % increase stepsize
                if (MA.cont.alpha<MA.cont.alpmax)&&(info.ftrial < fbest)
                    % step size calulation and best point updating
                    MA.cont.alpha = ...
                        min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                    % next point to be tested
                    if(MA.cont.alpha < MA.cont.alpmax)
                       info.ytrial     = xbest;  
                       % compute a new step size
                       info.beta = ...
                          min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                       % project trial point into [low upp]
                       info.statep='c'; info.p=MA.cont.RecomInfo.p;
                       [info]=projectStep(xbest,info);
                       state=21;
                       x=info.ytrial;
                       return; % return MATRSstep to compute f at x
                    else
                        state=25; break;
                    end
                else % go to update information
                    state=25; break;
                end
            end
           if state==21 % save x and f in the history
               info=UpdatePoints(info,mitune);
               state=20;
                MA.cont.SSext=[MA.cont.SSext info.beta];
           end
       end % end of extrapolation  
       if state==199 % opposite direction is tried
            MA.cont.RecomInfo.p=-MA.cont.RecomInfo.p;
            % compute alpmax and alp0
            info.ls='c'; info.step = 'r';
            [MA,info]=getalpha(MA,info,ctune);
            if MA.cont.alpmax>0 
               % compute a new step size
               info.ytrial = xbest;
               % compute trial point and project it into [low upp]
               info.statep='c'; info.beta=MA.cont.alp0;
               info.p=MA.cont.RecomInfo.p;
               [info]=projectStep(xbest,info);
               state=221;
               x=info.ytrial;
               MA.cont.SSext=info.beta;
               return; % return MATRSstep to compute f at x
           else
               state=25; % go to update information
           end
       end
       if state==221 % save x and f in the history
            info=UpdatePoints(info,mitune);
            MA.cont.SSext=[MA.cont.SSext info.beta];
            % check whether the first trial point can be 
            % the first trial point of extrapolation
            if info.ftrial<fbest
                MA.cont.extrapDone=1;
                info.enf=info.nf;
                % initialize step size 
                MA.cont.alpha = MA.cont.alp0; 
                info.ytrial     = xbest;
                % compute a new step size
                info.beta = min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                 % calculate a new trial recombination point
                 % compute trial point and project it into [low upp]
                 info.statep='c';info.p=MA.cont.RecomInfo.p;
                [info]=projectStep(xbest,info);
                state=22;
                x=info.ytrial;
                return; % return MATRSstep to compute f at x
            else
                state=25; % go to update information
            end
       end
       if state==22 % save x and f in the history
           info=UpdatePoints(info,mitune);
           state=23;
           MA.cont.SSext=[MA.cont.SSext info.beta];
       end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % start of extrapolation along opposite %
      % recombination mutation direction      %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       while state==23||state==24 % extrapolation is doing
            if state==23
                % increase stepsize
                if (MA.cont.alpha<MA.cont.alpmax) && (info.ftrial<fbest)
                    % MA.cont.alpha calulation and best point updating
                    MA.cont.alpha = ...
                         min(MA.cont.alpmax, ctune.cnu*MA.cont.alpha);
                    % next point to be tested
                    if(MA.cont.alpha < MA.cont.alpmax)
                       info.ytrial  = xbest; 
                       info.beta = ...
                           min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                       % compute trial point and project it into [low,upp]
                       info.statep='c'; info.p=MA.cont.RecomInfo.p;
                       [info]=projectStep(xbest,info);
                       state=24;
                       x=info.ytrial;
                       return; 
                    else  % go to upfdate information
                        state=25; break;
                    end
                else % go to upfdate information
                    state=25; break;
                end
            end
           if state==24 % save x and f in the history
              info=UpdatePoints(info,mitune);
              state=23; % go back to do more extrapolation
              MA.cont.SSext=[MA.cont.SSext info.beta];
           end
       end  % end of extrapolation
       %%%%%%%%%%%%%%%%%%%%%%%
       % update information  %
       %%%%%%%%%%%%%%%%%%%%%%%
      if state==25  % update information                         
          if MA.cont.extrapDone
                % choose a trial point with the lowest inexact
                % function value among all evaluated points by
                % extrapolation
                X = info.XF(info.enf:info.nf-1,1:dim)';
                F = info.XF(info.enf:info.nf-1,dim+1);
                [ftrial,ib]=min(F);
                ib=ib(1); ytrial = X(:,ib);
                if ftrial < fbest, fbest=ftrial; xbest=ytrial; end
                if prt>=2 
                   fprintf(['cRecom: fbest=',num2str(fbest),'\n'])
                end
                MA.cont.ar = MA.cont.SSext(ib);
                MA.cont.goodr=1;
          else
              MA.cont.ar = info.beta;
          end % end of cLSS
          if MA.cont.lambda>=1 % update M and sigma
             info.var='c'; 
             [MA]=updateM(MA,info,ctune);
             % update sigma
             pow = norm(MA.cont.s)/MA.cont.echi - 1; 
             MA.cont.sigma = ...
                    MA.cont.sigma*exp((MA.cont.c_s/MA.cont.ds)*pow);
             MA.cont.sigma = max(ctune.csigmamin,...
                             min(ctune.csigmamax,MA.cont.sigma));
          end
          MA.cont.good = MA.cont.good ||MA.cont.goodr;
          MA.cont.okTR = ~MA.cont.okRec||~MA.cont.goodr;
          state=26; % go to perform cTRS
      end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % cTRS: continuous trust-region strategy %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if state==26 % initialization for cTRS
           MA.cont.Del  = min(ctune.Delmax,MA.cont.sigma);
           MA.cont.okTR = MA.cont.okTR && info.nf>nc+2; info.cTRSdone=1;
           if MA.cont.okTR % perform integer trust-region algorithm
               if prt>=2
                 fprintf([num2str(MA.cont.cit),'th cTRS is done\n'])
               end
               MA.cont.goodt=0;  
               MA.cont.M = adjustY(MA.cont.M,ctune);
               MA.cont.G = MA.cont.M'*MA.cont.M; 
               state = 27; % go to perform cTRS
           else
               state=29; % cTRS is skipped
           end
       end
       while ismember(state,[27 28 277]) % loop for cTRS
            if state==27
                % choose the sample points
                XX = info.XF(info.nf-nc-2:info.nf-1,1:dim)';
                FF = info.XF(info.nf-nc-2:info.nf-1,dim+1)';
                MA.cont.cxbest = xbest(info.cI);
                % approximate gradient
                MA.cont.g  = fitGrad(nc,XX(info.cI,:),...
                             FF,MA.cont.cxbest,fbest,nc+1,ctune);
                if norm(MA.cont.g)==0
                   I = randperm(info.nf-1,nc+1);
                   XX = info.XF(I,1:dim)';
                   FF = info.XF(I,dim+1)';
                    MA.cont.g  = fitGrad(nc,XX(info.cI,:),...
                             FF,MA.cont.cxbest,fbest,nc+1,ctune);
                        
                    if norm(MA.cont.g)==0, state=29; break ;end    
                end 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % solve the trust-region subproblem by minq8 %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               [A,D]    = eig(MA.cont.G);
               diagD     = diag(D);
               ind = diagD==0;
               if ~isempty(ind), [Dmax]=max(diagD);diagD(ind)=Dmax;end
               data.gam = fbest;
               data.c   = MA.cont.g;
               data.D   = diagD; 
               data.b   = zeros(nc,1);
               data.A   = A;
               xl  = max(info.clow-MA.cont.cxbest,-MA.cont.Del);
               xu  = min(info.cupp-MA.cont.cxbest,MA.cont.Del);
               % call minq8
               warning off
               [MA.cont.p,~,~,~,~] = minq8 ...
                    (data,xl,xu,xbest(info.cI),10000,1e-8,0);
               warning on
               if norm(MA.cont.p)==0, 
                   ptr = xbest(info.cI)-XX(info.cI,end);
                   if norm(ptr)~=0
                      MA.cont.p = MA.cont.Del*(ptr/norm(ptr));
                   else
                       state=29; % skip cTRS
                       break;
                   end
               end
               state=277;
            end
            if state==277  % evaluate the trial point 
                info.ytrial     = xbest;
                info.statep='c'; info.beta=1;
                info.p = MA.cont.p;
                [info] = projectStep(xbest,info);
                state=28;
                x=info.ytrial;
                return; % return MATRSstep to compute f at x
            end
            if state==28 % save x and f in the history
                info=UpdatePoints(info,mitune);
                gcsc = max(-MA.cont.g'*MA.cont.p,...
                       sqrt(eps)*abs(MA.cont.g)'*abs(MA.cont.p));
                if gcsc==0
                    if fbest>info.ftrial, succcTR=0;
                    else, succcTR=1;
                    end
                else
                    mueff = (fbest-info.ftrial)/gcsc;
                    succcTR  = (mueff*abs(mueff-1) > ctune.czeta);
                end
                if info.cTRSdone
                  info.ncTR =info.ncTR+1; info.cTRSdone=0;
                end
                % decrease condition
                % Updating iterate and TR radius.
                if succcTR  % successful iteration
                   fbest = info.ftrial; xbest = info.ytrial; 
                    if prt>=2 
                       fprintf(['cTRS: fbest=',num2str(fbest),'\n'])
                    end
                   MA.cont.goodt=1;
                   % expand TR radius
                   MA.cont.Del = max(norm(MA.cont.p,inf),...
                                 ctune.ctheta)*MA.cont.Del;
                  MA.cont.succ_cTRS = MA.cont.succ_cTRS+1; 
                else % unsuccessful iteration
                    MA.cont.unsucc_cTRS = MA.cont.unsucc_cTRS+1;  
                    % reduce TR radius
                    MA.cont.Del = min(norm(MA.cont.p,inf),...
                                  MA.cont.Del)/ctune.ctheta;
                    if MA.cont.Del<=ctune.cDeltamin
                       MA.cont.good = MA.cont.good ||MA.cont.goodt;
                        state=29; 
                        break; 
                    end
                end
                % go back to form the trust-region subproblem, solve it
                state=27;
            end
       end % end of cTRS
       if state==29  
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % change affine scaling matrix %
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if norm(MA.cont.M,inf)>=ctune.mmax
              MA.cont.s = zeros(nc, 1); MA.cont.M = eye(nc); 
          end
          if ni>0, state=30; % go to perform iMATRS
          else, state=2; % go to next iteration of cMATRS
          end
       end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% iMATRS: integer MATRS %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ni>0  
        if state==30  % initialization for iMutation
             MA.int.iit = MA.int.iit+1;  
             if prt>=1
               disp('========================================')
               fprintf([num2str(MA.int.iit),'th iMATRS\n'])
               disp('========================================')
            end
            if prt>=2
                 fprintf([num2str(MA.int.iit),'th iMutation is done \n'])
             end
             if ~MA.int.good, MA.int.OffspringPop0=MA.int.MutInfoList;
             else, MA.int.OffspringPop0 =[];
             end
             info.xbest0i=xbest; 
             MA.int.goodm=0; MA.int.MutInfoList=[];
             state=31; 
             MA.int.im=1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% iMutation: integer mutation phase %%%%%%%%
        %%%%%%%%%%%%%%% iMutation calls iLSS %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while state<42 && state>=31 % loop for iMutation
            if state==31 
                % compute distribution and mutation directions
                MA.int.MutInfo.sd = MA.int.D(:,MA.int.im); 
                MA.int.MutInfo.p  = ...
                                round(MA.int.M*MA.int.MutInfo.sd);
                MA.int.ixbest = xbest(info.iI);    
                % check feasibility and find largest allowed step size
                info.ls ='i'; info.step='m';
                [MA,info]=getalpha(MA,info,itune);
                if ~MA.int.feasible
                     % go to try opposite direction
                     state=36;
                     info.ftrial=[];  info.ftrial0=[];
                else % compute the first trial point
                     info.ytrial  = xbest;
                     % compute trial point and project it into [low upp]
                     info.statep='i'; info.beta=MA.int.alp0;
                     info.p=MA.int.MutInfo.p;
                     [info]=projectStep(xbest,info);
                     x=info.ytrial;
                     sxf  = size(info.XF,1);
                     diff = ... 
                     (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                     [mn,ind]=min(sum(diff,2));
                     if (mn<=10^-16) % x is not a new feasible point
                        fval=info.XF(ind(1),info.dim+1); 
                        info.ftrial = fval;
                        state = 322;
                        if cputime-initTime>=secmax  % secmax reached
                           x=''; % this is needed to avoid ifinity loop
                                 % when no feasible point exists
                           return
                        end
                    else % x is a new point
                        state=32;
                        MA.int.SSext=info.beta;
                        return; % return MATRSstep to compute f at x
                    end
                end
            end
           if state==32 % save x and f in the history
               info=UpdatePoints(info,mitune);
               state=322;
               MA.int.SSext=[MA.int.SSext info.beta];
           end
           if state==322
                % check whether the first trial point can be 
                % the first trial point of extrapolation
                MA.int.extrapDone=0;
                if info.ftrial<fbest
                    MA.int.extrapDone=1;
                    info.enf = info.nf-1;
                    % initialize alpha 
                    MA.int.alpha = MA.int.alp0;
                    % calculate the second trial point
                    info.ytrial = xbest;
                    info.beta   = ...
                          min(MA.int.alpmax,itune.inu*MA.int.alpha);
                    % compute trial point and project it into [low upp]
                    info.statep='i'; info.p=MA.int.MutInfo.p;
                    [info]=projectStep(xbest,info);
                    x=info.ytrial;
                    sxf  = size(info.XF,1);
                    diff = ... 
                    (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                    [mn,ind]=min(sum(diff,2));
                    if (mn<=10^-16) % x is not a new feasible point
                       fval=info.XF(ind(1),info.dim+1); 
                       info.ftrial = fval;
                       state = 34;
                        if cputime-initTime>=secmax% secmax reached
                            x=''; % this is needed to avoid ifinity loop
                                 % when no feasible point exists
                            return
                        end
                    else % x is a new point
                        state = 33; 
                        return; % return MATRSstep to compute f at x
                    end
                else
                    state=36; % go to try opposite direction
                end
           end
           if state==33 % check whether x is a new point
                info=UpdatePoints(info,mitune);
                state=34; 
                MA.int.SSext=[MA.int.SSext info.beta];
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % start of extrapolation along %
           % integer mutation direction   %
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
           while state==34||state==35 % loop for extrapolation
                if state==34 
                    % increase stepsize
                    if (MA.int.alpha<MA.int.alpmax)&&(info.ftrial<fbest)
                        % step size calulation and best point updating
                        MA.int.alpha = ...
                           min(MA.int.alpmax,itune.inu*MA.int.alpha);
                        % next point to be tested
                        if(MA.int.alpha < MA.int.alpmax)
                           info.ytrial  = xbest; 
                           % compute a new step size
                           info.beta    = ...
                              min(MA.int.alpmax,itune.inu*MA.int.alpha);
                            % compute trial point and project it 
                            % into [low upp]
                           info.statep='i'; info.p=MA.int.MutInfo.p;
                           [info]=projectStep(xbest,info);
                            sxf  = size(info.XF,1);
                            diff = ... 
                            (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                            [mn,ind]=min(sum(diff,2));
                            if (mn<=10^-16) 
                              % x is not a new feasible point
                              fval=info.XF(ind(1),info.dim+1); 
                              info.ftrial = fval;  state=34;
                               if cputime-initTime>=secmax
                                   x=''; 
                                   % secmax reached
                                   % this is needed to avoid ifinity 
                                   % loop when no feasible point exists
                                  return
                               end
                            else % x is a new point
                                state=35;
                                % return MATRSstep to compute f at x
                                return; 
                            end
                        else
                             state=41; break; % go to update information
                        end
                    else
                         state=41; break; % go to update information
                    end 
                end
               if state==35 % check whether x is a new feasible point
                    info=UpdatePoints(info,mitune);
                    state=34; % go back to do more extrapolation
                    MA.int.SSext=[MA.int.SSext info.beta];
               end
           end
           if state==36 % opposite direction is tried
                MA.int.MutInfo.p = -MA.int.MutInfo.p; 
                MA.int.ixbest      = xbest(info.iI);       
                info.ls            = 'i'; info.step='m';
                [MA,info]          = getalpha(MA,info,itune);
                if ~MA.int.feasible
                    % go to update information
                    state=41;
                else % compute the first trial point along opposite dir
                    if ~isempty(info.ftrial)
                       info.ftrial0   = info.ftrial;
                    end
                     info.ytrial    = xbest;
                     % compute trial point and project it into [low upp]
                     info.statep='i'; info.beta=MA.int.alp0;
                     info.p=MA.int.MutInfo.p;
                     [info]=projectStep(xbest,info);
                     x=info.ytrial;
                     sxf  = size(info.XF,1);
                     diff = ... 
                     (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                     [mn,ind]=min(sum(diff,2));
                     if (mn<=10^-16) % x is not a new feasible point
                        fval=info.XF(ind(1),info.dim+1);
                        info.ftrial = fval;
                        state = 377;
                        if cputime-initTime>=secmax % secmax reached
                           x=''; % this is needed to avoid ifinity 
                                 % loop when no feasible point exists
                           return
                        end
                    else % x is a new point
                        state=37;
                        MA.int.SSext=info.beta;
                        return; % return MATRSstep to compute f at x
                    end
                end
           end
           if state==37 % check whether x is a new feasible point
                info=UpdatePoints(info,mitune);
                state=377;
                MA.int.SSext=[MA.int.SSext info.beta];
           end
           if state==377
                % check whether the first trial point can be 
                % the first trial point of extrapolation
                MA.int.extrapDone=0;
                if info.ftrial<fbest
                    MA.int.extrapDone=1;
                    info.enf = info.nf-1;
                    
                    MA.int.alpha = MA.int.alp0; 
                    info.ytrial = xbest;
                    info.beta   = ...
                            min(MA.int.alpmax,itune.inu*MA.int.alpha);
                    % compute trial point and project it into [low upp]
                    info.statep='i'; 
                    info.p=MA.int.MutInfo.p;
                    [info]=projectStep(xbest,info);
                     x=info.ytrial;
                     sxf  = size(info.XF,1);
                     diff = ... 
                     (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                     [mn,ind]=min(sum(diff,2));
                     if (mn<=10^-16) % x is not a new feasible point
                        fval=info.XF(ind(1),info.dim+1);
                        info.ftrial = fval;
                        state = 39;
                        if cputime-initTime>=secmax % secmax reached
                            x=''; % this is needed to avoid ifinity 
                            % loop when no feasible point exists
                           return
                        end
                    else % x is a new point
                        state=38;
                        return; % return MATRSstep to compute f at x
                    end
                else % no decrease in f
                     MA.int.MutInfo.sd = -MA.int.MutInfo.sd;
                     state=41; % go to update information
                     if ~isempty(info.ftrial0)
                         if info.ftrial0<info.ftrial
                             info.ftrial=info.ftrial0; 
                         end
                     end
                end
            end
            if state==38 % check whether x is a new feasible point
                info=UpdatePoints(info,mitune);
                state=39; 
                MA.int.SSext=[MA.int.SSext info.beta];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % start of extrapolation along %
            % opposite direction           %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            while state==39||state==40  % loop for extrapolation
                if state==39 % 
                    % expansion step (increase stepsize)
                    okstep = (MA.int.alpha<MA.int.alpmax && ...
                              info.ftrial < fbest);
                    if okstep
                        % MA.int.alpha calulation and best point updating
                        MA.int.alpha = ...
                           min(MA.int.alpmax, itune.inu*MA.int.alpha);
                        % next point to be tested
                        if(MA.int.alpha < MA.int.alpmax)
                           info.ytrial = xbest;  
                           info.beta  = ...
                              min(MA.int.alpmax,itune.inu*MA.int.alpha);
                           % compute trial point and project it 
                           % into [low upp]
                           info.statep='i'; info.p=MA.int.MutInfo.p;
                           [info]=projectStep(xbest,info);
                           x=info.ytrial;
                           sxf  = size(info.XF,1);
                           diff = ... 
                           (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                           [mn,ind]=min(sum(diff,2));
                           if (mn<=10^-16)% x is not a new feasible point
                              fval=info.XF(ind(1),info.dim+1); 
                              info.ftrial = fval;
                              state = 39;
                               if cputime-initTime>=secmax 
                                  % secmax reached
                                  x='';% this is needed to avoid ifinity 
                                  % loop when no feasible point exists
                                  return
                               end
                           else % x is a new point
                               % return MATRSstep to compute f at x
                               state=40;
                               return; 
                           end
                        else % go to update information
                           state = 41;
                           break;
                        end
                    else % go to update information
                       state = 41;
                       break;
                    end
                end
                if state==40 % check whether x is a new feasible point
                   info=UpdatePoints(info,mitune);
                   state=39; % go back to do more extrapolation
                   MA.int.SSext=[MA.int.SSext info.beta];
                end
            end
             if state==41 
                 % update information after ending extrapolation
                 if MA.int.extrapDone 
                      X = info.XF(info.enf:info.nf-1,1:dim)';
                      F = info.XF(info.enf:info.nf-1,dim+1);
                      [ftrial,ib]=min(F);
                      ib=ib(1); ytrial = X(:,ib);
                      if ftrial < fbest,fbest=ftrial;xbest=ytrial;end
                      if prt>=2 
                         fprintf(['iMutation: fbest=',...
                                  num2str(fbest),'\n'])
                      end
                      MA.int.a(MA.int.im) = MA.int.SSext(ib);
                      MA.int.goodm=1;
                 else
                      if MA.int.feasible
                         MA.int.a(MA.int.im) = info.beta; 
                      else
                         MA.int.a(MA.int.im) = MA.int.sigma;
                      end
                      MA.int.MutInfo.f    = info.ftrial;
                      MA.int.MutInfoList  = ...
                              [MA.int.MutInfoList MA.int.MutInfo];
                 end
                  MA.int.im=MA.int.im+1;
                  if MA.int.im>size(MA.int.D,2)
                      state=42; 
                      break;% end of iMutation
                  else % go back to the next iteration of iMutation
                      state=31; 
                  end
             end
        end % end of iMutation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update of the set of distribution directions %
        % by usequnce or random generator              %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if state==42 
            % generate a set of integer distribution directions
            MA = igeneratorD(xbest,ni,MA,itune,info);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % selection phase in iMATRS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            MA.int.lambda = length(MA.int.MutInfoList);
           if ~MA.int.good && MA.int.lambda<6
               MA.int.MutInfoList = ...
                   [MA.int.OffspringPop0,MA.int.MutInfoList];
           end
           MA.int.lambda = length(MA.int.MutInfoList);
           MA.int.good = MA.int.goodm;
           MA.int.okTR = ~MA.int.goodm; MA.int.goodr=0; 
           MA.int.mu   = ceil(MA.int.lambda/2);
           if MA.int.lambda>=1 
             % initialization before selection
             iwi_raw   = log(MA.int.lambda/2 + 0.5)-log((1:MA.int.mu));
             MA.int.wi = iwi_raw/sum(iwi_raw);
             imu_eff   = max(1+eps,1/sum(MA.int.wi .^2));
             MA.int.c_s  = min(1.999,(imu_eff + 2)/(ni + imu_eff + 5));
             MA.int.c_1  = 2/((ni+1.3)^2 + imu_eff);
             MA.int.c_mu = .... 
             min(1-MA.int.c_1,2*(imu_eff-2+1/imu_eff)/((ni+2)^2+imu_eff));
             MA.int.d_s      = ...
                    1+MA.int.c_s+2*max(0,sqrt((imu_eff-1)/(ni+1))-1);
             MA.int.sqrt_s = sqrt(MA.int.c_s*(2-MA.int.c_s)*imu_eff);
              if prt>=1
                  fprintf([num2str(MA.int.iit),...
                          'th selection is done\n'])
              end
             MA.int.permut = sortData(MA.int.MutInfoList, ordering);
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % iRecom: integer recombination phase %
           %%%%%%%%%% iRecom calls iLSS %%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if MA.int.mu>=2
              % compute recombination distribution and mutation directions
              MA.int.sumz = zeros(ni,1);
              MA.int.sumd = zeros(ni,1);
              for m = 1:MA.int.mu 
                 MA.int.sumz=MA.int.sumz+MA.int.wi(m)*...
                      MA.int.MutInfoList(MA.int.permut(m)).sd;
                 MA.int.sumd=MA.int.sumd+MA.int.wi(m)*...
                      MA.int.MutInfoList(MA.int.permut(m)).p;
              end
              MA.int.okRec = norm(MA.int.sumd)~=0 && ~MA.int.goodm; 
              state=43; % go to perform iRecombination
           else
              state=54; % skip iRecombination
           end
       end
       if state==43
          if MA.int.okRec % find largest allowed step szie
              if prt>=1
                 fprintf([num2str(MA.int.iit),...
                     'th iRecombination is done\n'])
              end
             MA.int.ixbest = xbest(info.iI);  ni = length(info.iI);
             okRecombonation=1; MA.int.goodr=0; ir=0; beta=1;
            while okRecombonation % scale recombintion mutation direction
                                  % at most nscale times as long as
                                  % there is no feasible point
               ir=ir+1; 
               if ir>1,beta = 2*rand(ni,1)-1;beta=beta/norm(beta);end
               MA.int.RecomInfo.p = ceil(beta.*MA.int.sumd);
               okRecombonation = (norm(MA.int.RecomInfo.p)~=0);
               if okRecombonation
                  info.ls='i'; info.step = 'r';
                  [MA,info]=getalpha(MA,info,itune);
                   if MA.int.alpmax>=1
                       state=44; % go to perform iRecombination
                       break;
                   end
               end
               if ir>=itune.nscale
                   state=54;  % skip iRecombination
                   break;
               else
                   state=44; % go to perform iRecombination
               end
            end
         else
             state=54;  % skip iRecombination
         end
      end
      if state==44 % compute the trial point
         info.ytrial  = xbest;
         % compute trial point and project it into [low upp]
         info.statep='i'; info.beta=MA.int.alp0;
         info.p=MA.int.RecomInfo.p;
         [info]=projectStep(xbest,info);
         x=info.ytrial;
         sxf  = size(info.XF,1);
         diff = ... 
         (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
         [mn,ind]=min(sum(diff,2));
         if (mn<=10^-16) % x is not a new feasible point
            fval=info.XF(ind(1),info.dim+1); info.ftrial = fval;
            state = 46;
            if cputime-initTime>=secmax % secmax reached
                x=''; % this is needed to avoid ifinity 
                % loop when no feasible point exists
               return
            end
        else % x is a new point
            state=45;
            MA.int.SSext=info.beta;
            return; % return MATRSstep to compute f at x
        end
      end
      if state==45  % check whether x is a new feasible point
        info=UpdatePoints(info,mitune);
        state=46; 
        MA.int.SSext=[MA.int.SSext,info.beta];
      end
      if state==46 
            % check whether the first trial point can be 
            % the first trial point of extrapolation
            MA.int.extrapDone=0;
            if info.ftrial<fbest % decrease in f was found
                MA.int.extrapDone=1;
                info.enf = info.nf-1;
                MA.int.alpha = MA.int.alp0; 
                info.ytrial  = xbest;
                % compute a new step size
                info.beta = min(MA.int.alpmax,itune.inu*MA.int.alpha);
                % compute trial point and project it into [low upp]
                info.statep='i'; info.p=MA.int.RecomInfo.p;
                [info]=projectStep(xbest,info);
                x=info.ytrial;
                sxf  = size(info.XF,1);
                diff = ... 
                (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                [mn,ind]=min(sum(diff,2));
                if (mn<=10^-16) % x is not a new feasible point
                    fval=info.XF(ind(1),info.dim+1); info.ftrial = fval;
                    state = 48;
                     if cputime-initTime>=secmax % secmax reached
                        x=''; % this is needed to avoid ifinity 
                        % loop when no feasible point exists
                       return
                    end
                else % x is a new point
                    state=47;
                    return; % return MATRSstep to compute f at x
                end
            else
                state=477; % go to try opposite direction
            end
      end
      if state==47 % check whether x is a new feasible point
        info=UpdatePoints(info,mitune);
        state=48; % go to do extrapolation
        MA.int.SSext=[MA.int.SSext,info.beta];
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % start of extrapolation along     %
      % recombination mutation direction %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      while state==48||state==49 % loop for extrapolation
            if state==48 
                % increase stepsize
                if (MA.int.alpha<MA.int.alpmax) && (info.ftrial<fbest)
                    MA.int.alpha = ...
                           min(MA.int.alpmax,itune.inu*MA.int.alpha);
                    % next point to be tested
                    if(MA.int.alpha < MA.int.alpmax)
                       info.ytrial     = xbest; 
                       info.beta = ...
                             min(MA.int.alpmax,itune.inu*MA.int.alpha);
                       % compute trial point and project it into [low upp]
                       info.statep='i'; info.p=MA.int.RecomInfo.p;
                       [info]=projectStep(xbest,info);
                       x=info.ytrial;
                       sxf  = size(info.XF,1);
                       diff = ... 
                       (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                       [mn,ind]=min(sum(diff,2));
                       if (mn<=10^-16) % x is not a new feasible point
                            fval=info.XF(ind(1),info.dim+1);
                            info.ftrial = fval;
                            state = 48;
                            if cputime-initTime>=secmax % secmax reached
                               x=''; % this is needed to avoid ifinity 
                               % loop when no feasible point exists
                               return
                            end
                       else % x is a new point
                            state=49;
                            return; % return MATRSstep to compute f at x
                       end
                    else
                        state=53; % go to update information
                        break; 
                    end
                else
                    state=53; % go to update information
                    break; 
                end
           end
           if state==49 % check whether x is a new feasible point
               info=UpdatePoints(info,mitune);
               state=48; % go back to continue extrapolation
                MA.int.SSext=[MA.int.SSext,info.beta];
           end
      end  % end of extrapolation
      if state==477 % opposite direction is tried
           MA.int.RecomInfo.p  = -MA.int.RecomInfo.p;
           info.ls='i'; info.step = 'r';
           [MA,info]=getalpha(MA,info,itune);
           if MA.int.alpmax>=1
               info.ytrial   = xbest;
               % compute trial point and project it into [low upp]
               info.statep='i'; info.beta=MA.int.alp0;
               info.p=MA.int.RecomInfo.p;
               [info]=projectStep(xbest,info);
               x=info.ytrial;
               sxf  = size(info.XF,1);
               diff = ... 
               (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
               [mn,ind]=min(sum(diff,2));
               if (mn<=10^-16) % x is not a new feasible point
                   fval=info.XF(ind(1),info.dim+1); info.ftrial=fval;
                   state = 1499;
                   if cputime-initTime>=secmax % secmax reached
                        x=''; % this is needed to avoid ifinity 
                        % loop when no feasible point exists
                        return
                   end
               else % x is a new point
                   state=149;
                   MA.int.SSext=info.beta;
                   return; % return MATRSstep to compute f at x
               end
           else
               state=53; % go to update information
           end
      end
      if state==149 % check whether x is a new trial point
         info=UpdatePoints(info,mitune);
         state=1499;
         MA.int.SSext=[MA.int.SSext,info.beta];
      end
      if  state==1499
            % check whether the first trial point can be 
            % the first trial point of extrapolation
            if info.ftrial<fbest % decrese in f was found
                MA.int.extrapDone=1;
                info.enf = info.nf-1;
                % initialize MA.int.alpha and best point
                MA.int.alpha = MA.int.alp0; 
                % calculate trial point
                info.ytrial  = xbest;
                info.beta  = min(MA.int.alpmax,itune.inu*MA.int.alpha);
                % compute trial point and project it into [low upp]
                info.statep='i';  info.p=MA.int.RecomInfo.p;
                [info]=projectStep(xbest,info);
                x=info.ytrial;
                sxf  = size(info.XF,1);
                diff = ... 
                (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
               [mn,ind]=min(sum(diff,2));
                if (mn<=10^-16) % x is not a new feasible point
                    fval=info.XF(ind(1),info.dim+1);info.ftrial=fval;
                    state = 151;
                     if cputime-initTime>=secmax % secmax reached
                        x=''; % this is needed to avoid ifinity 
                        % loop when no feasible point exists
                        return
                     end
                else % x is a new point
                    state=150;
                    return; % return MATRSstep to compute f at x
                end
            else
               state  = 53; % go to update information
            end
      end
      if state==150 % check whether x is a new feasible point
        info=UpdatePoints(info,mitune);
        state  = 151;
        MA.int.SSext=[MA.int.SSext,info.beta];
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % start of extrapolation along opposite %
      % recombination mutation direction      %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      while state==151||state==152 % loop for extrapolation
            if state==151 
                % expansion step (increase stepsize)
                if (MA.int.alpha<MA.int.alpmax) && (info.ftrial < fbest)
                    MA.int.alpha = ...
                         min(MA.int.alpmax, itune.inu*MA.int.alpha);
                    % next point to be tested
                    if(MA.int.alpha < MA.int.alpmax)
                       info.ytrial  = xbest;  
                       info.beta  = ...
                          min(MA.int.alpmax,itune.inu*MA.int.alpha);
                       % compute trial point and project it into [low,upp]
                       info.statep='i'; info.p=MA.int.RecomInfo.p;
                       [info]=projectStep(xbest,info);
                       x=info.ytrial;
                       sxf  = size(info.XF,1);
                       diff = ... 
                       (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                       [mn,ind]=min(sum(diff,2));
                       if (mn<=10^-16) % x is not a new feasible point
                          info.ftrial =info.XF(ind(1),info.dim+1);
                          state = 151;
                          if cputime-initTime>=secmax % secmax reached
                             x=''; % this is needed to avoid ifinity 
                             % loop when no feasible point exists
                             return
                          end
                       else % x is a new point
                           state=152;
                           return; % return MATRSstep to compute f at x
                       end
                    else
                       % go to update information
                       state=53; 
                       break;
                    end
                else
                   % go to update information
                   state=53; 
                   break;
                end
            end
           if state==152  % check whether x is a new feasible point
               info=UpdatePoints(info,mitune);
               state=151;
               MA.int.SSext=[MA.int.SSext,info.beta];
           end
      end  % end of extrapolation
      %%%%%%%%%%%%%%%%%%%%%%
      % update information %
      %%%%%%%%%%%%%%%%%%%%%%
      if state==53                
          if MA.int.extrapDone  % extrapolation was done
              % extrapolation was generated at least two trial points
              % one of which with with the lowest inexact
              % function value among all evaluated points by
              % extrapolation
              X = info.XF(info.enf:info.nf-1,1:dim)';
              F = info.XF(info.enf:info.nf-1,dim+1);
              [ftrial,ib]=min(F);
              ib=ib(1); ytrial = X(:,ib);
              if ftrial < fbest, fbest=ftrial; xbest=ytrial; end
              if prt>=2 
                 fprintf(['iRecom: fbest=',num2str(fbest),'\n'])
               end
              MA.int.goodr=1;
              MA.int.ar = MA.int.SSext(ib);
          else
              MA.int.ar = info.beta;
          end
          MA.int.good   = MA.int.good ||MA.int.goodr;
          MA.int.okTR   = ~MA.int.okRec||~MA.int.goodr;
          state  = 54; 
          if MA.int.lambda>=1 % update M and sigma
             info.var='i'; % integer
             [MA]= updateM(MA,info,itune);
              % update sigma
              pow = norm(MA.int.s)/MA.int.echi - 1;
              MA.int.sigma = ...
                       MA.int.sigma*exp((MA.int.c_s/MA.int.d_s)*pow);
              MA.int.sigma = ...
                      round(min(itune.isigmamax,max(1,MA.int.sigma)));
          end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % iTRS: integer trust-region strategy %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if state==54 % initialization for iTRS
           MA.int.Deli= ...
                     min(itune.Delmax,max(itune.Delmin,MA.int.sigma)); 
           MA.int.okTR = MA.int.okTR && info.nf>ni+2; info.iTRSdone=1;
           if MA.int.okTR % perform integer trust-region algorithm
                MA.int.istuck = 0; 
                MA.int.M = adjustY(MA.int.M,itune);
                if rank(MA.int.M) < size(MA.int.M,2)
                   state=57; % skip iTRS
                else
                    state = 55; % go to the loop of iTRS
                    warning off
                    MA.int.invM = inv(MA.int.M);
                    warning on
                    MA.int.changePTR=0;
                    if prt>=1
                       disp('=====================')
                       fprintf('iTRS is tried \n')
                       disp('======================')
                    end
                 end
           else
               state=57; % skip iTRS
           end
     end
     while ismember(state,[55 56 555 565]) % loop iTRS
            % pick sometimes randomly n+1 sample points 
            if state==55
                if MA.int.changePTR 
                   I = randperm(info.nf-1,ni+1);
                   XXi = info.XF(I,1:dim)';
                   FFi = info.XF(I,dim+1)';
                else % n+1 last evaluated point
                   XXi = info.XF(info.nf-ni-2:info.nf-1,1:dim)';
                   FFi = info.XF(info.nf-ni-2:info.nf-1,dim+1)';
                end
               MA.int.ixbest = xbest(info.iI);
               MA.int.g = ...
                     fitGrad(ni,XXi(info.iI,:),FFi,MA.int.ixbest,...
                                                     fbest,ni+1,itune);
              if norm(MA.int.g)==0
                   I = randperm(info.nf-1,ni+1);
                   XXi = info.XF(I,1:dim)';
                   FFi = info.XF(I,dim+1)';
                   MA.int.g =  ...
                         fitGrad(ni,XXi(info.iI,:),FFi,...
                                       MA.int.ixbest,fbest,ni+1,itune);
                    if norm(MA.int.g)==0, state=57; break ;end    
              end  
               ll  = max(info.ilow-MA.int.ixbest,-MA.int.Deli);
               uu  = min(info.iupp-MA.int.ixbest,MA.int.Deli);
               igls = 0;
                for i = 1 : ni
                   if ll(i) >= uu(i), igls=1; end
                end
                if ~igls
                   r = -MA.int.invM*MA.int.g;
                   U = qr([MA.int.M,r]);
                   R = triu(U(1:ni,1:ni));
                   y = U(1:ni,ni+1);
                   MA.int.p = obils_search(R,y,ll,uu);
                   if isempty(MA.int.p)||norm(MA.int.p)==0
                       ptr = xbest(info.iI)-XXi(info.iI,end);
                       if norm(prt)~=0
                          MA.int.p = MA.int.Deli*round(ptr/norm(ptr));
                       else
                           state=57; % skip iTRS
                           break;
                       end
                       if isempty(MA.int.p)||norm(MA.int.p)==0
                           state=57; % skip iTRS
                           break
                       end
                   end
                else % skip iTRS
                    state=57;  break
                end
                state=555; % go to compute a new feasible trial point
            end
            if state==555
                info.ytrial     = xbest;
                info.statep='i'; info.beta=1;
                info.p = MA.int.p;
                [info] = projectStep(xbest,info);
                sxf  = size(info.XF,1);
                diff = (info.XF(:,1:dim)-repmat(info.ytrial',sxf,1)).^2;
                [mn,~]     = min(sum(diff,2));
                if (mn>10^-16) % go to compute a new trial point
                   state=56; MA.int.changePTR=0;
                   x=info.ytrial; 
                   return; % return MATRSstep to compute f at x
                else
                    state=565; % the feasibe trial point already
                                % has been evaluated
                    if cputime-initTime>=secmax % secmax reached
                       x=''; % this is needed to avoid ifinity 
                       % loop when no feasible point exists
                       return
                    end
                end
            end
            if  state==565 % update TR radius
                 MA.int.istuck=MA.int.istuck+1; MA.int.changePTR=1;
                 if sign(rand-0.5)
                    MA.int.Deli = ...
                    abs(MA.int.Deli+sign(rand-0.5)*...
                                 randi([1,itune.Delmin],1));
                 else
                     MA.int.Deli = round(MA.int.Deli/...
                                      randi([1,itune.Delmin],1));
                 end
                 if MA.int.Deli<=0||MA.int.istuck>itune.stuckmax
                    state=57; % skip iTRS
                    break; 
                 end
                 state=55; % go back to the next iteration of iTRS
            end
            if state==56  % save x and f in the history
               info=UpdatePoints(info,mitune);
               gisci  = max(-MA.int.g'*MA.int.p,...
                             sqrt(eps)*abs(MA.int.g)'*abs(MA.int.p));
                if gisci==0
                    if fbest>info.ftrial, suciTR=0;
                    else, suciTR=1;
                    end
                else
                    mueff   = (fbest-info.ftrial)/gisci;
                    suciTR  = (mueff*abs(mueff-1) > itune.izeta);
                end
                if info.iTRSdone
                   info.niTR =info.niTR+1; info.iTRSdone=0;
                end
                if suciTR % iteration is successful
                   fbest = info.ftrial; xbest = info.ytrial;  
                   if prt>=2 
                       fprintf(['iTRS: fbest=',num2str(fbest),'\n'])
                   end
                   MA.int.goodt=1;
                   MA.int.good = MA.int.good ||MA.int.goodt;
                  % state=57;  
                   MA.int.succ_iTRS=MA.int.succ_iTRS+1; 
                   if MA.int.Deli<=itune.iDeltabar
                        MA.int.Deli=MA.int.Deli+1;
                    else
                        MA.int.Deli = ...
                            floor(itune.itheta*max(norm(MA.int.p,inf),...
                                          MA.int.Deli));
                    end
                else % iteration is unsuccessful
                    MA.int.unsucc_iTRS = MA.int.unsucc_iTRS+1;
                    if MA.int.Deli<=itune.iDeltabar
                        MA.int.Deli=MA.int.Deli-1;
                    else
                        MA.int.Deli = ...
                            floor(min(norm(MA.int.p,inf),...
                                          MA.int.Deli)/itune.itheta);
                    end
                end
                if MA.int.Deli<=0 % iTRS ends
                    state=57; 
                    break; 
                end
                % go back to form the trust-region subproblem and solve it
                state=55; 
            end
     end % end of iTRS
     if state==57
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % change affine scaling matrix %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
       if norm(MA.int.M,inf)>=itune.mmax
         MA.int.s = zeros(ni, 1); MA.int.M = eye(ni); 
       end
       if nc>0 && ni==0, state=2; % go back to cMATRS
       elseif ni>0 && nc==0, state = 30;  % go back to iMATRS
       elseif ni>0 && nc>0, state = 58; % go to miMATRS
       end
     end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% miMATRS: mixed-integer MATRS %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if state == 58 % compute combination direction
        p0 =  xbest-info.xbest0; nF =  length(info.F);
        if nF>=2
            for i=1:nF, dX(:,i) = info.X(:,i)-xbest;end
            sc = abs(max(dX')'); sc(sc==0) = 1;
            if norm(sc)==0, sc=ones(dim,1); end
        else
            sc=ones(dim,1);
        end
        sc(info.iI) = max(1,1./sc(info.iI));
        sc(info.cI) = min(1,1./sc(info.cI));
        changDic = norm(p0)==0;
        if ~changDic 
           p0          = sc.*p0;
           MA.cont.p   = p0(info.cI); 
           MA.int.p    = ceil(p0(info.iI));
           MA.mixed.p  = p0;
           MA.mixed.p(info.iI)=MA.int.p;
           if norm(MA.int.p)~=0
              % find largest allowed step sizes
              info.step='mi'; info.ls='i';
              MA.int.ixbest=xbest(info.iI);
              [MA,info]= getalpha(MA,info);
           end
           if norm(MA.cont.p)~=0
             % find largest allowed step sizes
             info.step='mi'; info.ls='c';
             MA.cont.cxbest=xbest(info.cI);
             [MA,info]= getalpha(MA,info);
           end
           if norm(MA.int.p)~=0||norm(MA.cont.p)~=0
             state=59; % go to mixed-integer phase
           else % skip mixed-integer phase and go back to cMATRS
             state=2; 
           end
        else
            MA.int.p=[]; MA.cont.p=[];
            for kk=1:length(info.F)-2
                p0        = xbest-info.X(:,end-kk);
               p0         = sc.*p0;
                MA.cont.p = p0(info.cI); 
                MA.int.p  = ceil(p0(info.iI));
                if norm(MA.int.p)~=0 || norm(MA.cont.p)~=0,
                     MA.mixed.p  = p0; MA.mixed.p(info.iI)=MA.int.p;
                    break;
                end
            end
            if (norm(MA.int.p)==0 && norm(MA.cont.p)==0) || ...
                    isempty(MA.int.p)&&isempty(MA.cont.p)
                % skip mixed-integer phase and go back to cMATRS
                state=2; 
            else
                 if norm(MA.int.p)~=0 || norm(MA.cont.p)~=0
                   if norm(MA.int.p)~=0
                      % find largest allowed step sizes
                      info.step='mi'; info.ls='i';
                      MA.int.ixbest=xbest(info.iI);
                      [MA,info]= getalpha(MA,info);
                   end
                   if norm(MA.cont.p)~=0
                     % find largest allowed step sizes
                     info.step='mi'; info.ls='c';
                     MA.cont.cxbest=xbest(info.cI);
                     [MA,info]= getalpha(MA,info);
                   end
                    state=59; % go to continue mixed-integer phase
                 else
                     % skip mixed-integer phase and go back to cMATRS
                     state=2; 
                 end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % miLSS: mixed-integer line search %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if state==59
      info.ytrial  = xbest;
      if MA.int.feasible 
         % compute trial point and project it into [low upp]
         info.statep='i'; info.p=MA.int.p; info.beta=MA.int.alp0;
         [info]=projectStep(xbest,info);
      end
      if MA.cont.feasible
         % compute trial point and project it into [low upp]
         info.statep='c'; info.p=MA.cont.p; info.beta=MA.cont.alp0;
         [info]=projectStep(xbest,info);
      end
      if MA.int.feasible||MA.cont.feasible
           x=info.ytrial;
           sxf  = size(info.XF,1);
           diff = ... 
           (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
           [mn,ind]=min(sum(diff,2));
           if (mn<=10^-16) 
                 % x is not a new feasible point
                fval=info.XF(ind(1),info.dim+1); info.ftrial = fval;
                state = 62;
           else % x is a new point
                state=60; return
           end
      else  % skip mixed-integer phase and go back to cMATRS
           state=2;
       end
   end
   if state==60 % check whether x is a new feasible point
       info.number_miMATRS =info.number_miMATRS+1;
        if prt>=1
           disp('===========================================')
           fprintf([num2str(info.number_miMATRS),'th mixed-integer\n'])
           disp('===========================================')
        end
        info=UpdatePoints(info,mitune);
        state=62; % go to do extrapolation
   end
   if state==62
        MA.mixed.extrapDone=0;
        if info.ftrial<fbest % extrapolation may be tried
            info.enf = info.nf-1;
            MA.mixed.extrapDone=1;
            if MA.int.feasible
                MA.cont.alpha = inf; MA.int.alpha = MA.int.alp0;
                % calculate trial point
                info.ytrial= xbest;
                info.beta  = min(MA.int.alpmax,itune.inu*MA.int.alpha);
                % compute trial point and project it into [low upp]
                info.statep='i'; info.p=MA.int.p; 
                [info]=projectStep(xbest,info);
            end
            if MA.cont.feasible
                MA.cont.alpha = MA.cont.alp0; MA.int.alpha = inf;
                % calculate trial point
                info.ytrial= xbest;
                info.beta = min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                % compute trial point and project it into [low upp]
                info.statep='c'; info.p=MA.cont.p; 
                [info]=projectStep(xbest,info);
            end
            if MA.int.feasible||MA.cont.feasible
               x=info.ytrial;
               sxf  = size(info.XF,1);
               diff = ... 
               (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
               [mn,ind]=min(sum(diff,2));
               if (mn<=10^-16) 
                    % x is not a new feasible point
                    fval=info.XF(ind(1),info.dim+1); info.ftrial = fval;
                    state=64;
               else % x is a new point
                    state=63; return
               end
            else % no descrease in f
                state=81; % go to update information
            end
        else
            state=71; % go to try opposite direction
        end
   end
   if state==63  % check whether x is a new feasible point
        info=UpdatePoints(info,mitune);
        state=64; % go to do extrapolation
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%
   % extrapolation along p %
   %%%%%%%%%%%%%%%%%%%%%%%%%
   while state==64||state==65 % loop for extrapolation
        if state==64
              okLS = ((MA.int.alpha<MA.int.alpmax || ...
                    MA.cont.alpha<MA.cont.alpmax) && ...
                    info.ftrial < fbest);
             if okLS
                info.ytrial  = xbest;  
                MA.cont.alpha = ...
                         min(MA.cont.alpmax, ctune.cnu*MA.cont.alpha); 
                if(MA.cont.alpha < MA.cont.alpmax)
                    info.beta = ...
                         min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                    % compute trial point and project it into [low upp]
                    info.statep='c'; info.p=MA.cont.p; 
                    [info]=projectStep(xbest,info);
                end
                MA.int.alpha = ...
                         min(MA.int.alpmax, itune.inu*MA.int.alpha);
                if (MA.int.alpha < MA.int.alpmax)  
                    info.beta  = ...
                          min(MA.int.alpmax,itune.inu*MA.int.alpha);
                    % project trial point into [low upp]
                    info.statep='i'; info.p=MA.int.p; 
                    [info]=projectStep(xbest,info);
                end
                if (MA.cont.alpha < MA.cont.alpmax)||...
                                         (MA.int.alpha < MA.int.alpmax) 
                   x=info.ytrial;
                   sxf  = size(info.XF,1);
                   diff = ... 
                   (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                   [mn,ind]=min(sum(diff,2));
                   if (mn<=10^-16)
                       % x is not a new feasible point
                       fval=info.XF(ind(1),info.dim+1); 
                       info.ftrial = fval;
                       state=64;
                   else % x is a new point
                       state=65; return
                   end
                else
                    state=81; % go to update information
                    break;
                end
             else
                 state=81; % go to update information
                 break;
             end
        end
       if state==65 % check whether x is a new feasible point
            info=UpdatePoints(info,mitune);
            state=64; % go to continue extrapolation
       end
   end % end of extrapolatio
   if state==71 % opposite direction is tried
      MA.int.p =- MA.int.p; MA.cont.p=-MA.cont.p; 
      if norm(MA.int.p)~=0
          % find largest allowed step sizes
          info.step='mi'; info.ls='i';
          MA.int.ixbest=xbest(info.iI);
          [MA,info]= getalpha(MA,info);
       end
       if norm(MA.cont.p)~=0
         % find largest allowed step sizes
         info.step='mi'; info.ls='c';
         MA.cont.cxbest=xbest(info.cI);
         [MA,info]= getalpha(MA,info);
       end
       % compute the first triap point along opposite direction
      info.ytrial = xbest;
       if MA.int.feasible
             MA.int.alpha = 0;
             % compute trial point and project it into [low upp]
             info.statep='i'; info.p=MA.int.p; info.beta=MA.int.alp0;
             [info]=projectStep(xbest,info);
       end
       if MA.cont.feasible
             MA.cont.alpha = 0; 
             % compute trial point and project it into [low upp]
             info.statep='c'; info.p=MA.cont.p; info.beta=MA.cont.alp0;
             [info]=projectStep(xbest,info);
       end
       if MA.cont.feasible ||MA.int.feasible % trial point is feasible
           x=info.ytrial;
           sxf  = size(info.XF,1);
           diff = ... 
           (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
           [mn,ind]=min(sum(diff,2));
           if (mn<=10^-16)
                % x is not a new feasible point
                fval=info.XF(ind(1),info.dim+1); info.ftrial = fval;
                state=73; % go to do extrapolation
           else % x is a new point
                state=72; return
           end
       else
           state=2; % end miMATRS and go back to cMATRS
       end
   end
   if state==72 % check whether x is a new point
        info=UpdatePoints(info,mitune);
        state=73; % go to do extrapolation
   end
   if state==73
       MA.mixed.extrapDone=0;
    if info.ftrial<fbest % the first trial point is the first 
                         % trial point of extrapolation
        MA.mixed.extrapDone=1;   
        info.enf = info.nf-1;
        % initialize alpha 
        info.ytrial= xbest;
       if MA.int.feasible
            MA.cont.alpha = inf; MA.int.alpha = MA.int.alp0;
            % compute step size
            info.beta      = min(MA.int.alpmax,itune.inu*MA.int.alpha);
            % compute trial point and project it into [low upp]
            info.statep='i'; info.p=MA.int.p; 
            [info]=projectStep(xbest,info);
       end
        if MA.cont.feasible
            MA.cont.alpha = MA.cont.alp0; MA.int.alpha = inf;
            % calculate step size
            info.beta = min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
            % project trial point into [low upp]
            info.statep='c'; info.p=MA.cont.p; 
            [info]=projectStep(xbest,info);
        end
        if MA.int.feasible||MA.cont.feasible
           x=info.ytrial;
           sxf  = size(info.XF,1);
           diff = ... 
           (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
           [mn,ind]=min(sum(diff,2));
           if (mn<=10^-16)
                % x is not a new feasible point
                fval=info.XF(ind(1),info.dim+1); info.ftrial = fval;
                state=75; % go to do extrapolation
           else % x is a new point
                state=74;
                return
           end
        else
            state=81; % go to update information
        end
    else
        state=81; % go to update information
    end
   end
   if state==74 % check whether x is a new feasible point
        info=UpdatePoints(info,mitune);
        state=75; % go to do extrapolation
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%
   % extrapolation along -p %
   %%%%%%%%%%%%%%%%%%%%%%%%%%
   while state==75||state==76 % loop for extrapolation
        if state==75
            % expansion step (increase stepsize)
            okLS = (MA.int.alpha<MA.int.alpmax || ...
                       MA.cont.alpha<MA.cont.alpmax) && ...
                                                info.ftrial < fbest;
            if okLS
                info.ytrial  = xbest;  
                 MA.cont.alpha = ...
                         min(MA.cont.alpmax, ctune.cnu*MA.cont.alpha); 
                if(MA.cont.alpha < MA.cont.alpmax)
                   
                    info.beta = ...
                         min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                    % compute trial point and project it into [low upp]
                    info.statep='c'; info.p=MA.cont.p; 
                    [info]=projectStep(xbest,info);
                end
                MA.int.alpha = ...
                         min(MA.int.alpmax, itune.inu*MA.int.alpha);
                if (MA.int.alpha < MA.int.alpmax)  
                   info.beta  = ...
                         min(MA.int.alpmax,itune.inu*MA.int.alpha);
                   % compute trial point and project it into [low upp]
                   info.statep='i'; info.p=MA.int.p; 
                   [info]=projectStep(xbest,info);
                end
                if (MA.cont.alpha < MA.cont.alpmax)...
                                      ||(MA.int.alpha < MA.int.alpmax) 
                       x=info.ytrial;
                       sxf  = size(info.XF,1);
                       diff = ... 
                       (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                       [mn,ind]=min(sum(diff,2));
                       if (mn<=10^-16)
                          % x is not a new feasible point
                          fval=info.XF(ind(1),info.dim+1); 
                          info.ftrial = fval;
                          state=75;
                       else % x is a new point
                          state=76; return
                       end 
                else
                   state=81; break; % go to update information
                end
            else
                state=81; break; % go to update information
            end
        end
        if state==76 % check whether x is a new feasible point
            info=UpdatePoints(info,mitune);
            state=75; % go back to do more extrapolation
        end
   end % end of extrapolation
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  update information after %
   %  the end of extrapolation %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if state==81 
     if MA.mixed.extrapDone % extrapolation 
           % choose a trial point with lowest 
           % f computed by extrapolation
           X = info.XF(info.enf:info.nf-1,1:dim)';
           F = info.XF(info.enf:info.nf-1,dim+1);
           [ftrial,ib]=min(F); ib=ib(1); ytrial = X(:,ib);
           if ftrial < fbest
               % save step size
               MA.mixed.ae = norm(xbest-ytrial)/norm(MA.mixed.p);
               % update the best point
               fbest=ftrial; xbest=ytrial;
           end
           if prt>=2 
              fprintf(['Mixed-integer: fbest=',num2str(fbest),'\n'])
           end
     end
     state=2; % go back to cMATRS
   end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%