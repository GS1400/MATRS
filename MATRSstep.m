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
% 
%
% To achieve the reentrant behavior, we need to record the place where 
% the algorithm should be continued. This is done using the variable 
% state, used in MATRSstep, which encodes the possible places where 
% a function value is required. state takes a total of 81 different  
% values:
%
%    initialization               %  state='initial' 
%
%    cMATRS
%             % cMutation         %  state takes 'cMut1',...'cMut12'
%             % selection         %  state takes 'cSel'
%             % cRecom            %  state takes 'cRec1',...'cRec13' 
%             % cTRS              %  state takes 'cTR1',...'cTR4'     
%             % reset M           %  state takes 'cMres'
%
%    iMATRS   
%             % iMutation         %  state takes 'iMut1',...'cMut14' 
%             % selection         %  state takes 'iSel'
%             % iRecom            %  state takes 'iRec1',...'iRec14' 
%             % iTRS              %  state takes 'iTR1',...'iTR5'
%             % reset M           %  state takes 'iMres'
%
%    miMATRS                      %  state takes 'miLS1',...'miLS14' 
%
% A detailed description of the meaning of each state is in MATRSstep.m.  



function [x,f,infos]=MATRSstep(x,f,p,indMixInt,low,upp,budget)

persistent MA state prt itune ctune mitune info ni nc dim 

persistent xbest fbest  xfinal ffinal ordering nfmax secmax initTime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MATRSread of the paper %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0 % return info 
  if ffinal<fbest, f=ffinal; x=xfinal;
  else, f=fbest; x=xbest;
  end
  infos.nf            =info.nf;
  infos.number_iTRS   =info.niTR; 
  infos.number_cTRS   =info.ncTR; 
  infos.number_miMATRS=info.number_miMATRS;
  infos.number_succ_cTRS =MA.cont.succ_cTRS;
  infos.number_succ_iTRS =MA.int.succ_iTRS;
  infos.state=state;
  return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MATRSinit of the paper %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout==0 % initialize solver environment
  info=[]; dim=size(indMixInt,2)*size(indMixInt,1);
  indMixInt=indMixInt(:); low=low(:); upp=upp(:);
  % problem dimension
  % subspace dimension of continuous variables
  info.cI=find(indMixInt==0); 
  % subspace dimension of integer variables
  info.iI=setdiff(1:dim,info.cI);
  nc=length(info.cI); ni=length(info.iI); 
  if nc>0 % lower and upper bound on continuous variables
    info.clow=low(info.cI); info.cupp=upp(info.cI);
  end
  if ni>0 % lower and upper bound on continuous variables
     info.ilow=low(info.iI); info.iupp=upp(info.iI);
  end
  if isempty(x)
     xc=[]; yc=[]; zc=[];
  else
     xc=x.cont; yc=x.int; zc=x.mint;
  end
  % tuning parameters
  [ctune,itune,mitune]=initTune(xc,yc,zc,nc,ni);
  % budgets for stopping tests
  if isfield(budget,'nfmax'), nfmax=budget.nfmax; 
  else, nfmax=budget;
  end
  if isfield(budget,'secmax'), secmax=budget.secmax; 
  else, secmax=inf;
  end
  if isfield(budget,'initTime'),initTime=budget.initTime; end
  info.nc=nc; info.ni=ni;  prt=p; infos=info; state='initial'; 
  return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% initialization %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(state,'initial') 
    xbest=x; fbest=f; info.nf=1; ordering='ascend';
    if ni>0 % initialization for cMATRS
        MA.int         =[];
        MA.int.sigma   =itune.sigmai;
        MA.int.ar      =itune.ainit;
        MA.int.Parent.y=x;
        MA.int.s       =zeros(ni, 1);
        MA.int.M       =itune.M;
        MA.int.echi    =sqrt(ni)*(1-1/ni/4-1/ni/ni/21);
        MA.int.a       =eye(itune.ilambda,1); 
        MA.int.kappa   =1;
        MA.int.D       =ones(ni,itune.ilambda);
        MA.int.good    =1;
        MA.int.succ_iTRS=0;
        MA.int.unsucc_iTRS=0;
    end
    if nc>0 % initialization for iMATRS
        MA.cont         =[];
        MA.cont.sigma   =ctune.sigmac;
        MA.cont.ar      =ctune.ainit;
        MA.cont.echi    =sqrt(nc)*(1-1/nc/4-1/nc/nc/21);
        MA.cont.Parent.y=x;
        MA.cont.s       =zeros(nc, 1);
        MA.cont.M       =ctune.M;
        MA.cont.kappa   =1;
        MA.cont.D       =ones(nc,ctune.clambda); 
        MA.cont.a       =ones(ctune.clambda,1);
        MA.cont.good    =1;
        MA.cont.succ_cTRS=0;
        MA.cont.unsucc_cTRS=0;
    end
    % create the first list for saving all evaluated points and their f
    info.XF=Inf*ones(nfmax,dim+1);
    info.XF(1,1:dim)=x'; info.XF(1,dim+1)=fbest;
    % create the second list, used to compute combination directions
    info.X=xbest; info.F=fbest; info.dim=dim;
    MA.int.iit=0; MA.cont.cit=0; 
    info.ncTR=0; % number of calls to cTRS
    info.niTR=0; % number of calls to iTRS
    info.number_miMATRS=0; % number of calls to miMATRS
    MA.int.succ_iTRS=0; MA.int.unsucc_iTRS=0;
    MA.cont.succ_cTRS=0; MA.cont.unsucc_cTRS=0;
    if nc>0 
       state='cMut1'; % go to cMATRS (initialization for cMutation)
    elseif nc==0 && ni>0
       state='iMut1'; % go to iMATRS (initialization for iMutation)
    end
    if ni>0&& nc>0, MA.mixed.ae=mitune.ainit;end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATRS alternates calls to cMATRS, iMATRS, and miMATRS in this order %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save new point and its function value
xfinal=x; ffinal=f; info.f=f; 

while 1 % main loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% cMATRS: contniuous MATRS %%%%%%%%%%%%%%%%%%%
    %%%% cMATRS includes cMutation, selection, cRecom, and cTRS %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nc>0  % start of cMATRS
       if strcmp(state,'cMut1') % initialization for cMutation
            MA.cont.cit=MA.cont.cit+1;  
           if prt>=1
               disp('========================================')
               fprintf([num2str(MA.cont.cit),'th cMATRS\n'])
               disp('========================================')
           end
           if prt>=1 
             fprintf([num2str(MA.cont.cit),'th cMutation is done \n'])
           end
           if ~MA.cont.good, MA.cont.OffspringPop0=MA.cont.MutInfoList;
           else, MA.cont.OffspringPop0=[];
           end
           info.xbest0 =xbest; 
           MA.cont.goodm =0; MA.cont.MutInfoList=[]; MA.cont.cm=1;
           state='cMut2'; % go to start cMutation
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % cMutation: continuous mutation phase %
       % cMutation calls cLSS                 %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       incMut = cell(1,11); num=2:12; incMut(:) = {'cMut'};
       for ii=1:11, incMut{ii}=[incMut{ii},num2str(num(ii))];end
       while ismember(state,incMut) % cMutation is performing`
            if strcmp(state,'cMut2') % start of cMutation
               % compute continuous dustribution and mutation directions
               MA.cont.MutInfo.sd=MA.cont.D(:,MA.cont.cm); 
               MA.cont.MutInfo.p =MA.cont.M*MA.cont.MutInfo.sd;
               MA.cont.cxbest    =xbest(info.cI);  
               % check feasibility
               % find largest real step size 
               % find  mutation step size
               info.ls='c'; info.step='m';
               [MA,info]=getalpha(MA,info,ctune);
               % check whether or not cMutation ends
               if ~MA.cont.feasible
                    state='cMut7'; % go to try opposite direction
                    info.ftrial=[]; info.ftrial0=[];
               else % cMutation does not stop
                    info.ytrial=xbest;
                    % compute the first trial point and project it 
                    % into [low upp]
                    info.statep='c'; info.beta=MA.cont.alp0;
                    info.p=MA.cont.MutInfo.p;
                    [info]=projectStep(xbest,info);
                    state='cMut3'; % go to initialize extrapolation
                    x=info.ytrial;
                    MA.cont.SSext=info.beta;
                    return; % return MATRSstep to compute f at x
               end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % cLSS: continuous line search %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(state,'cMut3') % f at x was computed;
                                     % initialize extrapolation
                % save new point and its function value in the history
                info=UpdatePoints(info,mitune);
                MA.cont.extrapDone=0;
                MA.cont.SSext=[MA.cont.SSext info.beta];
                % check whether the first trial point is 
                % the first trial point of extrapolation
                if info.ftrial<fbest 
                    info.enf=info.nf-1;
                    MA.cont.extrapDone=1;
                    MA.cont.alpha=MA.cont.alp0; 
                    % calculate the second trial point
                    % along the mutation direction
                    info.ytrial =xbest;
                    % compute the step size
                    info.beta =...
                          min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                    % compute trial point and project it into [low upp]
                    info.statep='c';info.p=MA.cont.MutInfo.p;
                    [info]=projectStep(xbest,info);
                    x=info.ytrial;
                    state='cMut4'; % go to compute f at x
                    return; % return MATRSstep to compute f at x
                else % no decrease in f 
                    state='cMut7'; % go to try opposite direction
                end
            end
            if strcmp(state,'cMut4') % save new point and its f 
                                     % in the history
                info=UpdatePoints(info,mitune);
                MA.cont.SSext=[MA.cont.SSext info.beta];
                state='cMut5'; % go to do extrapolation 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % extrapolation along    %
            % the mutation direction %
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            while strcmp(state,'cMut5')||strcmp(state,'cMut6') 
                % extrapolation is doing
                if strcmp(state,'cMut5')
                    % increase step size
                    if (MA.cont.alpha<MA.cont.alpmax) && ...
                                                (info.ftrial < fbest)
                        % compute new feasible point
                        MA.cont.alpha=... 
                        min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                        % chekc whether alpmax reaches or not
                        if(MA.cont.alpha < MA.cont.alpmax)
                           info.ytrial=xbest; 
                           % compute new step size
                           info.beta  =...
                           min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                           % compute trial point and project it 
                           % into [low upp]
                           info.statep='c';info.p=MA.cont.MutInfo.p;
                           [info]=projectStep(xbest,info);
                           x=info.ytrial;
                           state='cMut6'; % go to compute f at x
                           return; % return MATRSstep to compute f at x
                        end
                    end
                    state='cMut12'; % go to update information
                    break
                end
               if strcmp(state,'cMut6') % f at x was computed; 
                                        % save x and f in the history
                   info=UpdatePoints(info,mitune);
                   state='cMut5'; % go back to continue extrapolation
                   MA.cont.SSext=[MA.cont.SSext info.beta];
               end
            end % end of extrapolation
            if strcmp(state,'cMut7') % opposite direction is tried
                MA.cont.MutInfo.p  =-MA.cont.MutInfo.p; 
                MA.cont.cxbest     =xbest(info.cI);       
                info.ls            ='c'; info.step='m';
                [MA,info]          =getalpha(MA,info,ctune);
                if ~MA.cont.feasible % check termination of cMutation
                    state='cMut12'; % go to update information
                else % cMutation is peroming
                     % compute new trial point
                     if~isempty(info.ftrial),info.ftrial0=info.ftrial;end
                     info.ytrial=xbest;
                     % compute trial point and project it into [low upp]
                     info.statep='c'; info.beta=MA.cont.alp0;
                     info.p=MA.cont.MutInfo.p;
                     [info]=projectStep(xbest,info);
                     x=info.ytrial;
                     state='cMut8'; % go to compute f at x
                     MA.cont.SSext=info.beta; 
                    return; % return MATRSstep to compute f at x
                end
            end
            % check whether extrapolation can be done or not
            if strcmp(state,'cMut8') 
                  % f at x was computed; save x and f in the history
                  info=UpdatePoints(info,mitune);
                  MA.cont.SSext=[MA.cont.SSext info.beta];
                  % check whether the first trial point is 
                  % the first trial point of extrapolation
                 if info.ftrial<fbest
                    MA.cont.extrapDone=1; 
                    info.enf=info.nf-1;
                    % initialize step size and best point
                    MA.cont.alpha=MA.cont.alp0;
                    % calculate a new feasible trial point
                    info.ytrial =xbest;
                    info.beta   =... 
                           min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                    % compute trial point and project it into [low upp]
                    info.statep='c'; info.p=MA.cont.MutInfo.p;
                    [info]=projectStep(xbest,info);
                    x=info.ytrial;
                    state='cMut9'; % go to compute f at x
                    return; % return MATRSstep to compute f at x
                 else % no decrease in f along opposite direction
                     MA.cont.MutInfo.sd=-MA.cont.MutInfo.sd;
                     state='cMut12'; % go to update information
                      if ~isempty(info.ftrial0)
                         if info.ftrial0<info.ftrial
                             % function value at mutation point
                             info.ftrial=info.ftrial0; 
                         end
                      end
                 end
            end
            if strcmp(state,'cMut9') % save x and f in the history
                info=UpdatePoints(info,mitune);
                state='cMut10'; % go to do extrapolation
                 MA.cont.SSext=[MA.cont.SSext info.beta];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extrapoltion along opposite %
            % mutation direction          %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            while strcmp(state,'cMut10')||strcmp(state,'cMut11') 
                % loop for extrapolation
                if strcmp(state,'cMut10')
                    % expansion step (increase stepsize)
                    if (MA.cont.alpha<MA.cont.alpmax) && ...
                                              (info.ftrial < fbest)
                        % step size calulation and best point updating
                        MA.cont.alpha=...
                          min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                        % next point to be tested
                        if(MA.cont.alpha < MA.cont.alpmax)
                           info.ytrial=xbest;  
                           % compute new step size
                           info.beta=min(MA.cont.alpmax,...
                                         ctune.cnu*MA.cont.alpha);
                           % compute trial point and project it
                           % into [low upp]
                           info.statep='c'; info.p=MA.cont.MutInfo.p;
                           [info]=projectStep(xbest,info);
                           x=info.ytrial;
                           state='cMut11'; % go to compute f at x
                           return; % return MATRSstep to compute f at x
                        end
                    end
                    % no decrease in f along both p and -p
                    state='cMut12'; % go to update information
                    break
                end
               if strcmp(state,'cMut11') % save x and f in the history
                   info=UpdatePoints(info,mitune);
                   state='cMut10'; % go back to do more extrapolation
                   MA.cont.SSext=[MA.cont.SSext info.beta];
               end
            end % end of ectrapolation along opposite mutation direction
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update information after     %
            % termination of extrapolation %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(state,'cMut12') 
                 if MA.cont.extrapDone 
                      % extrapolation was done
                      % extrapolation with amost two trial points
                      % one of which with the lowest inexact f
                      % is accepted as a new best point
                      X=info.XF(info.enf:info.nf-1,1:dim)';
                      F=info.XF(info.enf:info.nf-1,dim+1);
                      [ftrial,ib]=min(F);
                      ib=ib(1); ytrial=X(:,ib);
                      if ftrial < fbest, fbest=ftrial; xbest=ytrial; end
                      if prt>=1 
                        fprintf(['cMutation: fbest=',...
                                                 num2str(fbest),'\n'])
                      end
                      MA.cont.a(MA.cont.cm)=MA.cont.SSext(ib);
                      MA.cont.goodm=1;
                      if prt>=2
                         fprintf('best point was updated\n')
                      end
                 else % extrapolation cannot find a reduction of f
                      % update list of extrapolation step sizes
                      if MA.cont.feasible
                           MA.cont.a(MA.cont.cm)=info.beta;
                      else
                          MA.cont.a(MA.cont.cm)=MA.cont.sigma;
                      end
                      MA.cont.MutInfo.f  =info.ftrial;
                      MA.cont.MutInfoList=[MA.cont.MutInfoList ...
                                              MA.cont.MutInfo];
                 end
                 % check whether cMutation ends or not
                  MA.cont.cm=MA.cont.cm+1;
                  if MA.cont.cm>size(MA.cont.D,2)
                      state='cSel'; % go to change the set 
                                    % distribution directions
                      break;% end of cMutation
                  else 
                      state='cMut2'; % go back to continue cMutation
                  end
            end
       end % end of cMutation           
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % computation of the set of    %
       % distribution directions by   %
       % usequnce or random generator %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if strcmp(state,'cSel') % change the set D of distr. directions
          pinit      =xbest(info.cI)-info.xbest0(info.cI);
          normpinit  =norm(pinit,inf);
          changDc=(normpinit==0);
          if ~changDc % use usequence to generate D
               MA.cont.pinit=pinit;
               MA.cont.pinit=ctune.sc*MA.cont.pinit/normpinit;
               MA.cont.D    =MA.cont.kappa*usequence(nc,...
                            ctune.clambda,zeros(nc,1),MA.cont.pinit,0);
               if MA.cont.kappa<ctune.ckmax
                   MA.cont.kappa=MA.cont.kappa+1;
               else, MA.cont.kappa=1;
               end
          else % use the normal distribution to generate D
              MA.cont.D =randn(nc,ctune.clambda);
              MA.cont.kappa=1;    
          end
          %%%%%%%%%%%%%%%%%%%
          % selection phase %
          %%%%%%%%%%%%%%%%%%%
          MA.cont.lambda=length(MA.cont.MutInfoList);
          if ~MA.cont.good && MA.cont.lambda<6
             MA.cont.MutInfoList=...
                        [MA.cont.OffspringPop0,MA.cont.MutInfoList];
          end
          clambda=length(MA.cont.MutInfoList);
          MA.cont.good=MA.cont.goodm;
          MA.cont.goodr=0; MA.cont.okTR=~MA.cont.goodm;
          MA.cont.mu   =ceil(clambda/2);
          if clambda>=1 % initialization before selection
             cwi_raw      =log(clambda/2 + 0.5) - log((1:MA.cont.mu));
             MA.cont.wi   =cwi_raw/sum(cwi_raw);
             cmu_eff      =max(1+eps,1/sum(MA.cont.wi .^2));
             MA.cont.c_s  =min(1.999,(cmu_eff + 2)/(nc + cmu_eff + 5));
             MA.cont.c_1  =2/((nc+1.3)^2 + cmu_eff);
             MA.cont.c_mu =.... 
                    min( 1-MA.cont.c_1, 2*(cmu_eff - 2 + 1/cmu_eff)/...
                    ((nc + 2)^2 + cmu_eff));
             MA.cont.ds   =1 + MA.cont.c_s + ...
                                 2*max(0,sqrt((cmu_eff-1)/(nc+1))-1);
             MA.cont.sqrt_s=sqrt(MA.cont.c_s*(2-MA.cont.c_s)*cmu_eff);
             if prt>=1
                 fprintf([num2str(MA.cont.cit),...
                                   'th selection is done \n'])
             end
             MA.cont.permut=sortData(MA.cont.MutInfoList, ordering);
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % cRecom: continuous recombination phase %
          % cRecom calls cLSS                      %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if MA.cont.mu>=2
            % compute recombination distribution direction MA.cont.sumz 
            %         recombination mutation directions MA.cont.sumd 
            MA.cont.sumz=zeros(nc,1); MA.cont.sumd=zeros(nc,1); 
            for m=1:MA.cont.mu 
                MA.cont.sumz=MA.cont.sumz+MA.cont.wi(m)*...
                          MA.cont.MutInfoList(MA.cont.permut(m)).sd;
                MA.cont.sumd=MA.cont.sumd+MA.cont.wi(m)*...
                          MA.cont.MutInfoList(MA.cont.permut(m)).p;
             end
             MA.cont.okRec=(norm(MA.cont.sumd)~=0 && ~MA.cont.goodm);
             state='cRec1'; % go to do cRecom
           else
             state='cTR1'; % skip cRecom
          end
       end
       if strcmp(state,'cRec1') 
         if MA.cont.okRec % start of cRecom
              if prt>=2
                 fprintf([num2str(MA.cont.cit),...
                          'th cRecom is done\n'])
              end
             MA.cont.cxbest=xbest(info.cI);  nc=length(info.cI);
            okRecombonation=1;  MA.cont.goodr=0; cr=0; beta=1;
            % find alpmax and alp0
            while okRecombonation % scale recombination direction at most
                                  % nscale times so that some feasible
                                  % trial points can be found
               cr=cr+1; 
               if cr>1,beta=2*rand(nc,1)-1; beta=beta/norm(beta);end
               MA.cont.RecomInfo.p=beta.*MA.cont.sumd;
               okRecombonation=(norm(MA.cont.RecomInfo.p)~=0);
               if okRecombonation % find alpmax and alp0
                   info.ls='c'; info.step='r';
                   [MA,info]=getalpha(MA,info,ctune);
                   if MA.cont.alpmax>0
                       state='cRec2'; % go to start cLSS by cRecom
                       break;
                   end
               end
               if cr>=ctune.nscale
                   state='cTR1'; % skip cRecom
                   break;
               end
            end
         else
             state='cTR1'; % skip cRecom
         end
       end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      % start of cLSS by cRecom %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      if strcmp(state,'cRec2') % compute the trial recombination point
        % Build first point for starting linesearch
        info.ytrial =xbest;
        % project trial point into [low upp]
        info.statep='c'; info.beta=MA.cont.alp0;
        info.p=MA.cont.RecomInfo.p;
        [info]=projectStep(xbest,info);
        MA.cont.SSext=info.beta;
        x=info.ytrial;
        state='cRec3'; % go to compute f at x
        return; % return MATRSstep to compute f at x 
      end
      if strcmp(state,'cRec3') % save x and f in the history
        info=UpdatePoints(info,mitune);
        state='cRec4'; % go to check whether a reduction in f is found
        MA.cont.SSext=[MA.cont.SSext info.beta];
      end
      if strcmp(state,'cRec4')
             % check whether the first trial point can be 
             % the first trial point of extrapolation
             MA.cont.extrapDone=0;
            if info.ftrial<fbest
                info.enf=info.nf-1;
                MA.cont.extrapDone=1;
                % initialize MA.cont.alpha and best point
                MA.cont.alpha=MA.cont.alp0; 
                % calculate trial point
                info.ytrial=xbest;
                info.beta=min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                % project trial point into [low upp]
                info.statep='c'; info.p=MA.cont.RecomInfo.p;
                [info]=projectStep(xbest,info);
                x=info.ytrial; 
                state='cRec5'; % go to compute f at x
                return; % return MATRSstep to compute f at x
            else
                state='cRec12'; % go to do opposite direction
            end
      end
      if strcmp(state,'cRec5') % save x and f in the history
         state='cRec6'; info=UpdatePoints(info,mitune);
         MA.cont.SSext=[MA.cont.SSext info.beta];
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % start of extrapolation along     %
      % recombination mutation direction %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       while strcmp(state,'cRec6')||strcmp(state,'cRec7') 
           % extrapolation is doing
            if strcmp(state,'cRec6')
                % increase stepsize
                if (MA.cont.alpha<MA.cont.alpmax)&&(info.ftrial < fbest)
                    % step size calulation and best point updating
                    MA.cont.alpha=...
                        min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                    % next point to be tested
                    if(MA.cont.alpha < MA.cont.alpmax)
                       info.ytrial    =xbest;  
                       % compute a new step size
                       info.beta=...
                          min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                       % project trial point into [low upp]
                       info.statep='c'; info.p=MA.cont.RecomInfo.p;
                       [info]=projectStep(xbest,info);
                       x=info.ytrial;
                       state='cRec7'; % go to compute f at x
                       return; % return MATRSstep to compute f at x
                    else
                        state='cRec11'; % go to update information
                        break;
                    end
                else 
                    state='cRec11'; % go to update information
                    break;
                end
            end
           if strcmp(state,'cRec7') % save x and f in the history
               info=UpdatePoints(info,mitune);
               state='cRec6'; % go back to continue extrapolation
                MA.cont.SSext=[MA.cont.SSext info.beta];
           end
       end % end of extrapolation  
       if strcmp(state,'cRec12') % opposite direction is tried
            MA.cont.RecomInfo.p=-MA.cont.RecomInfo.p;
            % compute alpmax and alp0
            info.ls='c'; info.step='r';
            [MA,info]=getalpha(MA,info,ctune);
            if MA.cont.alpmax>0 
               % compute a new step size
               info.ytrial=xbest;
               % compute trial point and project it into [low upp]
               info.statep='c'; info.beta=MA.cont.alp0;
               info.p=MA.cont.RecomInfo.p;
               [info]=projectStep(xbest,info);
               MA.cont.SSext=info.beta;
               x=info.ytrial;
               state='cRec13'; % go to compute f at x
               return; % return MATRSstep to compute f at x
           else
               state='cRec11'; % go to update information
           end
       end
       if strcmp(state,'cRec13') % save x and f in the history
            info=UpdatePoints(info,mitune);
            MA.cont.SSext=[MA.cont.SSext info.beta];
            % check whether the first trial point can be 
            % the first trial point of extrapolation
            if info.ftrial<fbest
                MA.cont.extrapDone=1;
                info.enf=info.nf;
                % initialize step size 
                MA.cont.alpha=MA.cont.alp0; 
                info.ytrial    =xbest;
                % compute a new step size
                info.beta=min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                 % calculate a new trial recombination point
                 % compute trial point and project it into [low upp]
                 info.statep='c';info.p=MA.cont.RecomInfo.p;
                [info]=projectStep(xbest,info);
                state='cRec8';
                x=info.ytrial;
                return; % return MATRSstep to compute f at x
            else
                state='cRec11'; % go to update information
            end
       end
       if strcmp(state,'cRec8') % save x and f in the history
           info=UpdatePoints(info,mitune);
           state='cRec9'; % go to do extrapolation along -p
           MA.cont.SSext=[MA.cont.SSext info.beta];
       end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % start of extrapolation along opposite %
      % recombination mutation direction      %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       while strcmp(state,'cRec9')||strcmp(state,'cRec10') 
           % extrapolation is doing
            if strcmp(state,'cRec9')
                % increase stepsize
                if (MA.cont.alpha<MA.cont.alpmax) && (info.ftrial<fbest)
                    % MA.cont.alpha calulation and best point updating
                    MA.cont.alpha=...
                         min(MA.cont.alpmax, ctune.cnu*MA.cont.alpha);
                    % next point to be tested
                    if(MA.cont.alpha < MA.cont.alpmax)
                       info.ytrial =xbest; 
                       info.beta=...
                         min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                       % compute trial point and project it into [low,upp]
                       info.statep='c'; info.p=MA.cont.RecomInfo.p;
                       [info]=projectStep(xbest,info);
                       x=info.ytrial;
                       state='cRec10'; % go to compute f at x
                       return; 
                    else  
                        state='cRec11'; % go to upfdate information
                        break;
                    end
                else 
                    state='cRec11'; % go to upfdate information
                    break;
                end
            end
           if strcmp(state,'cRec10') % save x and f in the history
              info=UpdatePoints(info,mitune);
              state='cRec9'; % go back to do more extrapolation
              MA.cont.SSext=[MA.cont.SSext info.beta];
           end
       end  % end of extrapolation
       %%%%%%%%%%%%%%%%%%%%%%%
       % update information  %
       %%%%%%%%%%%%%%%%%%%%%%%
      if strcmp(state,'cRec11')  % update information                         
          if MA.cont.extrapDone
                % choose a trial point with the lowest inexact
                % function value among all evaluated points by
                % extrapolation
                X=info.XF(info.enf:info.nf-1,1:dim)';
                F=info.XF(info.enf:info.nf-1,dim+1);
                [ftrial,ib]=min(F);
                ib=ib(1); ytrial=X(:,ib);
                if ftrial < fbest, fbest=ftrial; xbest=ytrial; end
                if prt>=2 
                   fprintf(['cRecom: fbest=',num2str(fbest),'\n'])
                end
                MA.cont.ar=MA.cont.SSext(ib);
                MA.cont.goodr=1;
          else
              MA.cont.ar=info.beta;
          end % end of cLSS
          if MA.cont.lambda>=1 % update M and sigma
             info.var='c'; 
             [MA]=updateM(MA,info,ctune);
             % update sigma
             pow=norm(MA.cont.s)/MA.cont.echi - 1; 
             MA.cont.sigma=...
                    MA.cont.sigma*exp((MA.cont.c_s/MA.cont.ds)*pow);
             MA.cont.sigma=max(ctune.csigmamin,...
                             min(ctune.csigmamax,MA.cont.sigma));
          end
          MA.cont.good=MA.cont.good ||MA.cont.goodr;
          MA.cont.okTR=~MA.cont.okRec||~MA.cont.goodr;
          state='cTR1'; % go to perform cTRS
      end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % cTRS: continuous trust-region strategy %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if strcmp(state,'cTR1') % initialization for cTRS
           MA.cont.Del =min(ctune.Delmax,MA.cont.sigma);
           MA.cont.okTR=MA.cont.okTR && info.nf>nc+2; info.cTRSdone=1;
           if MA.cont.okTR % perform integer trust-region algorithm
               if prt>=2
                 fprintf([num2str(MA.cont.cit),'th cTRS is done\n'])
               end
               MA.cont.goodt=0;  
               MA.cont.M=adjustY(MA.cont.M,ctune);
               MA.cont.G=MA.cont.M'*MA.cont.M; 
               state='cTR2'; % go to perform cTRS
           else
              state='cMres'; % cTRS is skipped and go to rest M
           end
       end
       while ismember(state,{'cTR2' 'cTR3' 'cTR4'}) % loop for cTRS
            if strcmp(state,'cTR2')
                % choose the sample points
                XX=info.XF(info.nf-nc-2:info.nf-1,1:dim)';
                FF=info.XF(info.nf-nc-2:info.nf-1,dim+1)';
                MA.cont.cxbest=xbest(info.cI);
                % approximate gradient
                MA.cont.g =fitGrad(nc,XX(info.cI,:),...
                             FF,MA.cont.cxbest,fbest,nc+1,ctune);
                if norm(MA.cont.g)==0
                   I=randperm(info.nf-1,nc+1);
                   XX=info.XF(I,1:dim)';
                   FF=info.XF(I,dim+1)';
                    MA.cont.g =fitGrad(nc,XX(info.cI,:),...
                             FF,MA.cont.cxbest,fbest,nc+1,ctune);
                        
                    if norm(MA.cont.g)==0,
                        state='cMres'; % go to rest M
                        break;
                    end    
                end 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % solve the trust-region subproblem by minq8 %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               [A,D]   =eig(MA.cont.G);
               diagD    =diag(D);
               ind=diagD==0;
               if ~isempty(ind), [Dmax]=max(diagD);diagD(ind)=Dmax;end
               data.gam=fbest;
               data.c  =MA.cont.g;
               data.D  =diagD; 
               data.b  =zeros(nc,1);
               data.A  =A;
               xl =max(info.clow-MA.cont.cxbest,-MA.cont.Del);
               xu =min(info.cupp-MA.cont.cxbest,MA.cont.Del);
               % call minq8
               warning off
               [MA.cont.p,~,~,~,~]=minq8 ...
                    (data,xl,xu,xbest(info.cI),10000,1e-8,0);
               warning on
               if norm(MA.cont.p)==0, 
                   ptr=xbest(info.cI)-XX(info.cI,end);
                   if norm(ptr)~=0
                      MA.cont.p=MA.cont.Del*(ptr/norm(ptr));
                   else
                      state='cMres'; % skip cTRS and go to rest M
                       break;
                   end
               end
               state='cTR3'; % go to evaluate the trial point 
            end
            if strcmp(state,'cTR3')  % evaluate the trial point 
                info.ytrial    =xbest;
                info.statep='c'; info.beta=1;
                info.p=MA.cont.p;
                [info]=projectStep(xbest,info);
                x=info.ytrial;
                state='cTR4'; % go to compute f at x
                return; % return MATRSstep to compute f at x
            end
            if strcmp(state,'cTR4') % save x and f in the history
                info=UpdatePoints(info,mitune);
                gcsc=max(-MA.cont.g'*MA.cont.p,...
                       sqrt(eps)*abs(MA.cont.g)'*abs(MA.cont.p));
                if gcsc==0
                    if fbest>info.ftrial, succcTR=0;
                    else, succcTR=1;
                    end
                else
                    mueff=(fbest-info.ftrial)/gcsc;
                    succcTR =(mueff*abs(mueff-1) > ctune.czeta);
                end
                if info.cTRSdone
                  info.ncTR=info.ncTR+1; info.cTRSdone=0;
                end
                % decrease condition
                % Updating iterate and TR radius.
                if succcTR  % successful iteration
                   fbest=info.ftrial; xbest=info.ytrial; 
                    if prt>=2 
                       fprintf(['cTRS: fbest=',num2str(fbest),'\n'])
                    end
                   MA.cont.goodt=1;
                   % expand TR radius
                   MA.cont.Del=max(norm(MA.cont.p,inf),...
                                 ctune.ctheta)*MA.cont.Del;
                  MA.cont.succ_cTRS=MA.cont.succ_cTRS+1; 
                else % unsuccessful iteration
                    MA.cont.unsucc_cTRS=MA.cont.unsucc_cTRS+1;  
                    % reduce TR radius
                    MA.cont.Del=min(norm(MA.cont.p,inf),...
                                  MA.cont.Del)/ctune.ctheta;
                    if MA.cont.Del<=ctune.cDeltamin
                       MA.cont.good=MA.cont.good ||MA.cont.goodt;
                       state='cMres'; % go to rest M
                        break; 
                    end
                end
                state='cTR2'; % go back to form the trust-region
                              % subproblem, solve it
            end
       end % end of cTRS
       if strcmp(state,'cMres')  
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % change affine scaling matrix %
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if norm(MA.cont.M,inf)>=ctune.mmax
              MA.cont.s=zeros(nc, 1); MA.cont.M=eye(nc); 
          end
          if ni>0,state='iMut1'; % go to perform iMATRS
          else, state='cMut1'; % go to next iteration of cMATRS
          end
       end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% iMATRS: integer MATRS %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ni>0  
        if strcmp(state,'iMut1')  % initialization for iMutation
             MA.int.iit=MA.int.iit+1;  
             if prt>=1
               disp('========================================')
               fprintf([num2str(MA.int.iit),'th iMATRS\n'])
               disp('========================================')
            end
            if prt>=2
                 fprintf([num2str(MA.int.iit),'th iMutation is done \n'])
             end
             if ~MA.int.good, MA.int.OffspringPop0=MA.int.MutInfoList;
             else, MA.int.OffspringPop0=[];
             end
             info.xbest0i=xbest; 
             MA.int.goodm=0; MA.int.MutInfoList=[];
             state='iMut2'; % go to start iMutation
             MA.int.im=1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% iMutation: integer mutation phase %%%%%%%%
        %%%%%%%%%%%%%%% iMutation calls iLSS %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iniMut = cell(1,13); num=2:14; iniMut(:) = {'iMut'};
        for ii=1:13, iniMut{ii}=[iniMut{ii},num2str(num(ii))];end
        while ismember(state,iniMut) % loop for iMutation
            if strcmp(state,'iMut2') 
                % compute distribution and mutation directions
                MA.int.MutInfo.sd=MA.int.D(:,MA.int.im); 
                MA.int.MutInfo.p =...
                                round(MA.int.M*MA.int.MutInfo.sd);
                MA.int.ixbest=xbest(info.iI);    
                % check feasibility and find largest allowed step size
                info.ls='i'; info.step='m';
                [MA,info]=getalpha(MA,info,itune);
                if ~MA.int.feasible
                     state='iMut8'; % go to try opposite direction
                     info.ftrial=[];  info.ftrial0=[];
                else % compute the first trial point
                     info.ytrial =xbest;
                     % compute trial point and project it into [low upp]
                     info.statep='i'; info.beta=MA.int.alp0;
                     info.p=MA.int.MutInfo.p;
                     [info]=projectStep(xbest,info);
                     x=info.ytrial;
                     sxf =size(info.XF,1);
                     diff=... 
                     (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                     [mn,ind]=min(sum(diff,2));
                     if (mn<=10^-16) % x is not a new feasible point
                        fval=info.XF(ind(1),info.dim+1); 
                        info.ftrial=fval;
                        state='iMut4'; % go to reuse f and check whether
                                       % a reduction in f is found
                        if cputime-initTime>=secmax  % secmax reached
                           x=''; % this is needed to avoid ifinity loop
                                 % when no feasible point exists
                           return
                        end
                    else % x is a new point
                        MA.int.SSext=info.beta;
                        state='iMut3'; % go to compute f at x
                        return; % return MATRSstep to compute f at x
                    end
                end
            end
           if strcmp(state,'iMut3')  % save x and f in the history
               info=UpdatePoints(info,mitune);
               state='iMut4'; % go to possibly start an extrapolation
               MA.int.SSext=[MA.int.SSext info.beta];
           end
           if strcmp(state,'iMut4') 
                % check whether the first trial point can be 
                % the first trial point of extrapolation
                MA.int.extrapDone=0;
                if info.ftrial<fbest
                    MA.int.extrapDone=1;
                    info.enf=info.nf-1;
                    % initialize alpha 
                    MA.int.alpha=MA.int.alp0;
                    % calculate the second trial point
                    info.ytrial=xbest;
                    info.beta  =...
                          min(MA.int.alpmax,itune.inu*MA.int.alpha);
                    % compute trial point and project it into [low upp]
                    info.statep='i'; info.p=MA.int.MutInfo.p;
                    [info]=projectStep(xbest,info);
                    x=info.ytrial;
                    sxf =size(info.XF,1);
                    diff=... 
                    (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                    [mn,ind]=min(sum(diff,2));
                    if (mn<=10^-16) % x is not a new feasible point
                       fval=info.XF(ind(1),info.dim+1); 
                       info.ftrial=fval;
                       state='iMut6'; % go to reuse f and check whether
                                      % a reduction in f is found
                        if cputime-initTime>=secmax% secmax reached
                            x=''; % this is needed to avoid ifinity loop
                                 % when no feasible point exists
                            return
                        end
                    else % x is a new point
                        state='iMut5'; % go to compute f at x
                        return; % return MATRSstep to compute f at x
                    end
                else
                    state='iMut8'; % go to try opposite direction
                end
           end
           if strcmp(state,'iMut5') % check whether x is a new point
                info=UpdatePoints(info,mitune);
                state='iMut6'; % go to do extrapolation
                MA.int.SSext=[MA.int.SSext info.beta];
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % start of extrapolation along %
           % integer mutation direction   %
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
           while strcmp(state,'iMut6')||strcmp(state,'iMut7') 
               % loop for extrapolation
                if strcmp(state,'iMut6') 
                    % increase stepsize
                    if (MA.int.alpha<MA.int.alpmax)&&(info.ftrial<fbest)
                        % step size calulation and best point updating
                        MA.int.alpha=...
                           min(MA.int.alpmax,itune.inu*MA.int.alpha);
                        % next point to be tested
                        if(MA.int.alpha < MA.int.alpmax)
                           info.ytrial =xbest; 
                           % compute a new step size
                           info.beta   =...
                              min(MA.int.alpmax,itune.inu*MA.int.alpha);
                            % compute trial point and project it 
                            % into [low upp]
                           info.statep='i'; info.p=MA.int.MutInfo.p;
                           [info]=projectStep(xbest,info);
                            sxf =size(info.XF,1);
                            diff=... 
                            (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                            [mn,ind]=min(sum(diff,2));
                            if (mn<=10^-16) 
                              % x is not a new feasible point
                              fval=info.XF(ind(1),info.dim+1); 
                              info.ftrial=fval;  state='iMut6';
                               if cputime-initTime>=secmax
                                   x=''; 
                                   % secmax reached
                                   % this is needed to avoid ifinity 
                                   % loop when no feasible point exists
                                  return
                               end
                            else % x is a new point
                                state='iMut7'; % go to compute f at x
                                % return MATRSstep to compute f at x
                                return; 
                            end
                        else
                             state='iMut14'; % go to update information
                             break; 
                        end
                    else
                         state='iMut14'; % go to update information
                         break; 
                    end 
                end
               if strcmp(state,'iMut7') % check whether x is 
                                        % a new feasible point
                    info=UpdatePoints(info,mitune);
                    state='iMut6'; % go back to do more extrapolation
                    MA.int.SSext=[MA.int.SSext info.beta];
               end
           end
           if strcmp(state,'iMut8') % opposite direction is tried
                MA.int.MutInfo.p=-MA.int.MutInfo.p; 
                MA.int.ixbest     =xbest(info.iI);       
                info.ls           ='i'; info.step='m';
                [MA,info]         =getalpha(MA,info,itune);
                if ~MA.int.feasible
                    state='iMut14'; % go to update information
                else % compute the first trial point along opposite dir
                    if ~isempty(info.ftrial)
                       info.ftrial0  =info.ftrial;
                    end
                     info.ytrial   =xbest;
                     % compute trial point and project it into [low upp]
                     info.statep='i'; info.beta=MA.int.alp0;
                     info.p=MA.int.MutInfo.p;
                     [info]=projectStep(xbest,info);
                     x=info.ytrial;
                     sxf =size(info.XF,1);
                     diff=... 
                     (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                     [mn,ind]=min(sum(diff,2));
                     if (mn<=10^-16) % x is not a new feasible point
                        fval=info.XF(ind(1),info.dim+1);
                        info.ftrial=fval;
                        state='iMut10'; % go to reuse f from the history
                        if cputime-initTime>=secmax % secmax reached
                           x=''; % this is needed to avoid ifinity 
                                 % loop when no feasible point exists
                           return
                        end
                    else % x is a new point
                        MA.int.SSext=info.beta;
                        state='iMut9'; % go to compute f at x
                        return; % return MATRSstep to compute f at x
                    end
                end
           end
           if strcmp(state,'iMut9') % check whether x is
                                    % a new feasible point
                info=UpdatePoints(info,mitune);
                state='iMut10'; % go to possibly do extrapolation
                MA.int.SSext=[MA.int.SSext info.beta];
           end
           if strcmp(state,'iMut10')
                % check whether the first trial point can be 
                % the first trial point of extrapolation
                MA.int.extrapDone=0;
                if info.ftrial<fbest
                    MA.int.extrapDone=1;
                    info.enf=info.nf-1;
                    MA.int.alpha=MA.int.alp0; 
                    info.ytrial=xbest;
                    info.beta  =...
                            min(MA.int.alpmax,itune.inu*MA.int.alpha);
                    % compute trial point and project it into [low upp]
                    info.statep='i'; 
                    info.p=MA.int.MutInfo.p;
                    [info]=projectStep(xbest,info);
                     x=info.ytrial;
                     sxf =size(info.XF,1);
                     diff=... 
                     (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                     [mn,ind]=min(sum(diff,2));
                     if (mn<=10^-16) % x is not a new feasible point
                        fval=info.XF(ind(1),info.dim+1);
                        info.ftrial=fval;
                        state='iMut12'; % go to reuse f form the history
                        if cputime-initTime>=secmax % secmax reached
                            x=''; % this is needed to avoid ifinity 
                            % loop when no feasible point exists
                           return
                        end
                    else % x is a new point
                       state='iMut11'; % go to compute f at x
                        return; % return MATRSstep to compute f at x
                    end
                else % no decrease in f
                     MA.int.MutInfo.sd=-MA.int.MutInfo.sd;
                     state='iMut14'; % go to update information
                     if ~isempty(info.ftrial0)
                         if info.ftrial0<info.ftrial
                             info.ftrial=info.ftrial0; 
                         end
                     end
                end
            end
            if strcmp(state,'iMut11') % check whether x is 
                                      % a new feasible point
                info=UpdatePoints(info,mitune);
                state='iMut12'; % go to possibly fo extrapolation along -p
                MA.int.SSext=[MA.int.SSext info.beta];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % start of extrapolation along %
            % opposite direction           %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            while strcmp(state,'iMut12')||strcmp(state,'iMut13') 
                % loop for extrapolation
                if strcmp(state,'iMut12') % 
                    % expansion step (increase stepsize)
                    okstep=(MA.int.alpha<MA.int.alpmax && ...
                              info.ftrial < fbest);
                    if okstep
                        % MA.int.alpha calulation and best point updating
                        MA.int.alpha=...
                           min(MA.int.alpmax, itune.inu*MA.int.alpha);
                        % next point to be tested
                        if(MA.int.alpha < MA.int.alpmax)
                           info.ytrial=xbest;  
                           info.beta =...
                              min(MA.int.alpmax,itune.inu*MA.int.alpha);
                           % compute trial point and project it 
                           % into [low upp]
                           info.statep='i'; info.p=MA.int.MutInfo.p;
                           [info]=projectStep(xbest,info);
                           x=info.ytrial;
                           sxf =size(info.XF,1);
                           diff=... 
                           (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                           [mn,ind]=min(sum(diff,2));
                           if (mn<=10^-16)% x is not a new feasible point
                              fval=info.XF(ind(1),info.dim+1); 
                              info.ftrial=fval;
                              % go to reuse f from the history
                              state='iMut12';
                               if cputime-initTime>=secmax 
                                  % secmax reached
                                  x='';% this is needed to avoid ifinity 
                                  % loop when no feasible point exists
                                  return
                               end
                           else % x is a new point
                               state='iMut13'; % go to compute f at x
                               % return MATRSstep to compute f at x
                               return; 
                           end
                        else 
                           state='iMut14'; % go to update information
                           break;
                        end
                    else 
                       state='iMut14'; % go to update information
                       break;
                    end
                end
                if strcmp(state,'iMut13') % check whether x is 
                                          % a new feasible point
                   info=UpdatePoints(info,mitune);
                   state='iMut12'; % go back to do more extrapolation
                   MA.int.SSext=[MA.int.SSext info.beta];
                end
            end
             if strcmp(state,'iMut14')
                 % update information after ending extrapolation
                 if MA.int.extrapDone 
                      X=info.XF(info.enf:info.nf-1,1:dim)';
                      F=info.XF(info.enf:info.nf-1,dim+1);
                      [ftrial,ib]=min(F);
                      ib=ib(1); ytrial=X(:,ib);
                      if ftrial < fbest,fbest=ftrial;xbest=ytrial;end
                      if prt>=2 
                         fprintf(['iMutation: fbest=',...
                                  num2str(fbest),'\n'])
                      end
                      MA.int.a(MA.int.im)=MA.int.SSext(ib);
                      MA.int.goodm=1;
                 else
                      if MA.int.feasible
                         MA.int.a(MA.int.im)=info.beta; 
                      else
                         MA.int.a(MA.int.im)=MA.int.sigma;
                      end
                      MA.int.MutInfo.f   =info.ftrial;
                      MA.int.MutInfoList =...
                              [MA.int.MutInfoList MA.int.MutInfo];
                 end
                  MA.int.im=MA.int.im+1;
                  if MA.int.im>size(MA.int.D,2)
                      state='iSel'; 
                      break;% end of iMutation
                  else % go back to the next iteration of iMutation
                      state='iMut2'; 
                  end
             end
        end % end of iMutation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update of the set of distribution directions %
        % by usequnce or random generator              %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if strcmp(state,'iSel') 
            % generate a set of integer distribution directions
            MA=igeneratorD(xbest,ni,MA,itune,info);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % selection phase in iMATRS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            MA.int.lambda=length(MA.int.MutInfoList);
           if ~MA.int.good && MA.int.lambda<6
               MA.int.MutInfoList=...
                   [MA.int.OffspringPop0,MA.int.MutInfoList];
           end
           MA.int.lambda=length(MA.int.MutInfoList);
           MA.int.good=MA.int.goodm;
           MA.int.okTR=~MA.int.goodm; MA.int.goodr=0; 
           MA.int.mu  =ceil(MA.int.lambda/2);
           if MA.int.lambda>=1 
             % initialization before selection
             iwi_raw  =log(MA.int.lambda/2 + 0.5)-log((1:MA.int.mu));
             MA.int.wi=iwi_raw/sum(iwi_raw);
             imu_eff  =max(1+eps,1/sum(MA.int.wi .^2));
             MA.int.c_s =min(1.999,(imu_eff + 2)/(ni + imu_eff + 5));
             MA.int.c_1 =2/((ni+1.3)^2 + imu_eff);
             MA.int.c_mu=.... 
             min(1-MA.int.c_1,2*(imu_eff-2+1/imu_eff)/((ni+2)^2+imu_eff));
             MA.int.d_s     =...
                    1+MA.int.c_s+2*max(0,sqrt((imu_eff-1)/(ni+1))-1);
             MA.int.sqrt_s=sqrt(MA.int.c_s*(2-MA.int.c_s)*imu_eff);
              if prt>=1
                  fprintf([num2str(MA.int.iit),...
                          'th selection is done\n'])
              end
             MA.int.permut=sortData(MA.int.MutInfoList, ordering);
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % iRecom: integer recombination phase %
           %%%%%%%%%% iRecom calls iLSS %%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if MA.int.mu>=2
              % compute recombination distribution and mutation directions
              MA.int.sumz=zeros(ni,1);
              MA.int.sumd=zeros(ni,1);
              for m=1:MA.int.mu 
                 MA.int.sumz=MA.int.sumz+MA.int.wi(m)*...
                      MA.int.MutInfoList(MA.int.permut(m)).sd;
                 MA.int.sumd=MA.int.sumd+MA.int.wi(m)*...
                      MA.int.MutInfoList(MA.int.permut(m)).p;
              end
              MA.int.okRec=norm(MA.int.sumd)~=0 && ~MA.int.goodm; 
              state='iRecom1'; % go to perform iRecombination
           else
              state='iTR1'; % skip iRecombination
           end
       end
       if strcmp(state,'iRecom1') 
          if MA.int.okRec % find largest allowed step szie
              if prt>=1
                 fprintf([num2str(MA.int.iit),...
                     'th iRecombination is done\n'])
              end
             MA.int.ixbest=xbest(info.iI);  ni=length(info.iI);
             okRecombonation=1; MA.int.goodr=0; ir=0; beta=1;
            while okRecombonation % scale recombintion mutation direction
                                  % at most nscale times as long as
                                  % there is no feasible point
               ir=ir+1; 
               if ir>1,beta=2*rand(ni,1)-1;beta=beta/norm(beta);end
               MA.int.RecomInfo.p=ceil(beta.*MA.int.sumd);
               okRecombonation=(norm(MA.int.RecomInfo.p)~=0);
               if okRecombonation
                  info.ls='i'; info.step='r';
                  [MA,info]=getalpha(MA,info,itune);
                   if MA.int.alpmax>=1
                       state='iRecom2'; % go to perform iRecombination
                       break;
                   end
               end
               if ir>=itune.nscale
                   state='iTR1';  % skip iRecombination
                   break;
               else
                   state='iRecom2'; % go to perform iRecombination
               end
            end
         else
             state='iTR1';  % skip iRecombination
         end
      end
      if strcmp(state,'iRecom2')  % compute the trial point
         info.ytrial =xbest;
         % compute trial point and project it into [low upp]
         info.statep='i'; info.beta=MA.int.alp0;
         info.p=MA.int.RecomInfo.p;
         [info]=projectStep(xbest,info);
         x=info.ytrial;
         sxf =size(info.XF,1);
         diff=... 
         (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
         [mn,ind]=min(sum(diff,2));
         if (mn<=10^-16) % x is not a new feasible point
            fval=info.XF(ind(1),info.dim+1); info.ftrial=fval;
            state='iRecom4'; % go to reuse f from the history
            if cputime-initTime>=secmax % secmax reached
                x=''; % this is needed to avoid ifinity 
                % loop when no feasible point exists
               return
            end
        else % x is a new point
            MA.int.SSext=info.beta;
            state='iRecom3'; % go to compute f at x
            return; % return MATRSstep to compute f at x
        end
      end
      if strcmp(state,'iRecom3')  % check whether x is 
                                  % a new feasible point
        info=UpdatePoints(info,mitune);
        state='iRecom4'; % go to possibly do extrapolation
        MA.int.SSext=[MA.int.SSext,info.beta];
      end
      if strcmp(state,'iRecom4') 
            % check whether the first trial point can be 
            % the first trial point of extrapolation
            MA.int.extrapDone=0;
            if info.ftrial<fbest % decrease in f was found
                MA.int.extrapDone=1;
                info.enf=info.nf-1;
                MA.int.alpha=MA.int.alp0; 
                info.ytrial =xbest;
                % compute a new step size
                info.beta=min(MA.int.alpmax,itune.inu*MA.int.alpha);
                % compute trial point and project it into [low upp]
                info.statep='i'; info.p=MA.int.RecomInfo.p;
                [info]=projectStep(xbest,info);
                x=info.ytrial;
                sxf =size(info.XF,1);
                diff=... 
                (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                [mn,ind]=min(sum(diff,2));
                if (mn<=10^-16) % x is not a new feasible point
                    fval=info.XF(ind(1),info.dim+1); info.ftrial=fval;
                    state='iRecom6'; % go to reuse f from the history
                     if cputime-initTime>=secmax % secmax reached
                        x=''; % this is needed to avoid ifinity 
                        % loop when no feasible point exists
                       return
                    end
                else % x is a new point
                    state='iRecom5'; % go to compute f at x
                    return; % return MATRSstep to compute f at x
                end
            else
                state='iRecom8'; % go to try opposite direction
            end
      end
      if strcmp(state,'iRecom5')  % check whether x is 
                                  % a new feasible point
        info=UpdatePoints(info,mitune);
        state='iRecom6'; % go to do extrapolation
        MA.int.SSext=[MA.int.SSext,info.beta];
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % start of extrapolation along     %
      % recombination mutation direction %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      while strcmp(state,'iRecom6')||strcmp(state,'iRecom7') 
          % loop for extrapolation
            if strcmp(state,'iRecom6') 
                % increase stepsize
                if (MA.int.alpha<MA.int.alpmax) && (info.ftrial<fbest)
                    MA.int.alpha=...
                           min(MA.int.alpmax,itune.inu*MA.int.alpha);
                    % next point to be tested
                    if(MA.int.alpha < MA.int.alpmax)
                       info.ytrial    =xbest; 
                       info.beta=...
                             min(MA.int.alpmax,itune.inu*MA.int.alpha);
                       % compute trial point and project it into [low upp]
                       info.statep='i'; info.p=MA.int.RecomInfo.p;
                       [info]=projectStep(xbest,info);
                       x=info.ytrial;
                       sxf =size(info.XF,1);
                       diff=... 
                       (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                       [mn,ind]=min(sum(diff,2));
                       if (mn<=10^-16) % x is not a new feasible point
                            fval=info.XF(ind(1),info.dim+1);
                            info.ftrial=fval;
                            state='iRecom6'; % go to reuse f from 
                                             % the history
                            if cputime-initTime>=secmax % secmax reached
                               x=''; % this is needed to avoid ifinity 
                               % loop when no feasible point exists
                               return
                            end
                       else % x is a new point
                            state='iRecom7'; % go to compute f at x
                            return; % return MATRSstep to compute f at x
                       end
                    else
                        state='iRecom14'; % go to update information
                        break; 
                    end
                else
                    state='iRecom14'; % go to update information
                    break; 
                end
           end
           if strcmp(state,'iRecom7') % check whether x is 
                                      % a new feasible point
               info=UpdatePoints(info,mitune);
               state='iRecom6'; % go back to continue extrapolation
                MA.int.SSext=[MA.int.SSext,info.beta];
           end
      end  % end of extrapolation
      if strcmp(state,'iRecom8') % opposite direction is tried
           MA.int.RecomInfo.p =-MA.int.RecomInfo.p;
           info.ls='i'; info.step='r';
           [MA,info]=getalpha(MA,info,itune);
           if MA.int.alpmax>=1
               info.ytrial  =xbest;
               % compute trial point and project it into [low upp]
               info.statep='i'; info.beta=MA.int.alp0;
               info.p=MA.int.RecomInfo.p;
               [info]=projectStep(xbest,info);
               x=info.ytrial;
               sxf =size(info.XF,1);
               diff=... 
               (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
               [mn,ind]=min(sum(diff,2));
               if (mn<=10^-16) % x is not a new feasible point
                   fval=info.XF(ind(1),info.dim+1); info.ftrial=fval;
                   state='iRecom10'; % go to reuse f from the history
                   if cputime-initTime>=secmax % secmax reached
                        x=''; % this is needed to avoid ifinity 
                        % loop when no feasible point exists
                        return
                   end
               else % x is a new point
                   MA.int.SSext=info.beta;
                   state='iRecom9'; % go to compute f at x
                   return; % return MATRSstep to compute f at x
               end
           else
               state='iRecom14'; % go to update information
           end
      end
      if strcmp(state,'iRecom9') % check whether x is a new trial point
         info=UpdatePoints(info,mitune);
         state='iRecom10'; % go to possibly do extrapolation
         MA.int.SSext=[MA.int.SSext,info.beta];
      end
      if  strcmp(state,'iRecom10')
            % check whether the first trial point can be 
            % the first trial point of extrapolation
            if info.ftrial<fbest % decrese in f was found
                MA.int.extrapDone=1;
                info.enf=info.nf-1;
                % initialize MA.int.alpha and best point
                MA.int.alpha=MA.int.alp0; 
                % calculate trial point
                info.ytrial =xbest;
                info.beta =min(MA.int.alpmax,itune.inu*MA.int.alpha);
                % compute trial point and project it into [low upp]
                info.statep='i';  info.p=MA.int.RecomInfo.p;
                [info]=projectStep(xbest,info);
                x=info.ytrial;
                sxf =size(info.XF,1);
                diff=... 
                (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
               [mn,ind]=min(sum(diff,2));
                if (mn<=10^-16) % x is not a new feasible point
                    fval=info.XF(ind(1),info.dim+1);info.ftrial=fval;
                    state='iRecom12'; % go to reuse f from the history
                     if cputime-initTime>=secmax % secmax reached
                        x=''; % this is needed to avoid ifinity 
                        % loop when no feasible point exists
                        return
                     end
                else % x is a new point
                    state='iRecom11'; % go to compute f at x
                    return; % return MATRSstep to compute f at x
                end
            else
               state='iRecom14'; % go to update information
            end
      end
      if strcmp(state,'iRecom11') % check whether x is 
                                  % a new feasible point
        info=UpdatePoints(info,mitune);
        state='iRecom12';
        MA.int.SSext=[MA.int.SSext,info.beta];
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % start of extrapolation along opposite %
      % recombination mutation direction      %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      while strcmp(state,'iRecom12')||strcmp(state,'iRecom13') 
          % loop for extrapolation
            if strcmp(state,'iRecom12') 
                % expansion step (increase stepsize)
                if (MA.int.alpha<MA.int.alpmax) && (info.ftrial < fbest)
                    MA.int.alpha=...
                         min(MA.int.alpmax, itune.inu*MA.int.alpha);
                    % next point to be tested
                    if(MA.int.alpha < MA.int.alpmax)
                       info.ytrial =xbest;  
                       info.beta =...
                          min(MA.int.alpmax,itune.inu*MA.int.alpha);
                       % compute trial point and project it into [low,upp]
                       info.statep='i'; info.p=MA.int.RecomInfo.p;
                       [info]=projectStep(xbest,info);
                       x=info.ytrial;
                       sxf =size(info.XF,1);
                       diff=... 
                       (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                       [mn,ind]=min(sum(diff,2));
                       if (mn<=10^-16) % x is not a new feasible point
                          info.ftrial=info.XF(ind(1),info.dim+1);
                          state='iRecom12'; % go to reuse known f 
                          if cputime-initTime>=secmax % secmax reached
                             x=''; % this is needed to avoid ifinity 
                             % loop when no feasible point exists
                             return
                          end
                       else % x is a new point
                           state='iRecom13'; % go to compute f at x
                           return; % return MATRSstep to compute f at x
                       end
                    else
                       % go to update information
                       state='iRecom14'; 
                       break;
                    end
                else
                   % go to update information
                   state='iRecom14'; 
                   break;
                end
            end
           if strcmp(state,'iRecom13')  % check whether x is 
                                        % a new feasible point
               info=UpdatePoints(info,mitune);
               state='iRecom12';
               MA.int.SSext=[MA.int.SSext,info.beta];
           end
      end  % end of extrapolation
      %%%%%%%%%%%%%%%%%%%%%%
      % update information %
      %%%%%%%%%%%%%%%%%%%%%%
      if strcmp(state,'iRecom14')              
          if MA.int.extrapDone  % extrapolation was done
              % extrapolation was generated at least two trial points
              % one of which with with the lowest inexact
              % function value among all evaluated points by
              % extrapolation
              X=info.XF(info.enf:info.nf-1,1:dim)';
              F=info.XF(info.enf:info.nf-1,dim+1);
              [ftrial,ib]=min(F);
              ib=ib(1); ytrial=X(:,ib);
              if ftrial < fbest, fbest=ftrial; xbest=ytrial; end
              if prt>=2 
                 fprintf(['iRecom: fbest=',num2str(fbest),'\n'])
               end
              MA.int.goodr=1;
              MA.int.ar=MA.int.SSext(ib);
          else
              MA.int.ar=info.beta;
          end
          MA.int.good  =MA.int.good ||MA.int.goodr;
          MA.int.okTR  =~MA.int.okRec||~MA.int.goodr;
          state='iTR1'; 
          if MA.int.lambda>=1 % update M and sigma
             info.var='i'; % integer
             [MA]=updateM(MA,info,itune);
              % update sigma
              pow=norm(MA.int.s)/MA.int.echi - 1;
              MA.int.sigma=...
                       MA.int.sigma*exp((MA.int.c_s/MA.int.d_s)*pow);
              MA.int.sigma=...
                      round(min(itune.isigmamax,max(1,MA.int.sigma)));
          end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % iTRS: integer trust-region strategy %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if strcmp(state,'iTR1')  % initialization for iTRS
           MA.int.Deli=...
                     min(itune.Delmax,max(itune.Delmin,MA.int.sigma)); 
           MA.int.okTR=MA.int.okTR && info.nf>ni+2; info.iTRSdone=1;
           if MA.int.okTR % perform integer trust-region algorithm
                MA.int.istuck=0; 
                MA.int.M=adjustY(MA.int.M,itune);
                if rank(MA.int.M) < size(MA.int.M,2)
                   state='iMres'; % skip iTRS
                else
                    state='iTR2'; % go to the loop of iTRS
                    warning off
                    MA.int.invM=inv(MA.int.M);
                    warning on
                    MA.int.changePTR=0;
                    if prt>=1
                       disp('=====================')
                       fprintf('iTRS is tried \n')
                       disp('======================')
                    end
                 end
           else
               state='iMres'; % skip iTRS
           end
      end
      
     while ismember(state,{'iTR2' 'iTR3' 'iTR4' 'iTR5'}) % loop iTRS
            % pick sometimes randomly n+1 sample points 
            if strcmp(state,'iTR2') 
                if MA.int.changePTR 
                   I=randperm(info.nf-1,ni+1);
                   XXi=info.XF(I,1:dim)';
                   FFi=info.XF(I,dim+1)';
                else % n+1 last evaluated point
                   XXi=info.XF(info.nf-ni-2:info.nf-1,1:dim)';
                   FFi=info.XF(info.nf-ni-2:info.nf-1,dim+1)';
                end
               MA.int.ixbest=xbest(info.iI);
               MA.int.g=...
                     fitGrad(ni,XXi(info.iI,:),FFi,MA.int.ixbest,...
                                                     fbest,ni+1,itune);
              if norm(MA.int.g)==0
                   I=randperm(info.nf-1,ni+1);
                   XXi=info.XF(I,1:dim)';
                   FFi=info.XF(I,dim+1)';
                   MA.int.g= ...
                         fitGrad(ni,XXi(info.iI,:),FFi,...
                                       MA.int.ixbest,fbest,ni+1,itune);
                    if norm(MA.int.g)==0, state='iMres'; break ;end    
              end  
               ll =max(info.ilow-MA.int.ixbest,-MA.int.Deli);
               uu =min(info.iupp-MA.int.ixbest,MA.int.Deli);
               igls=0;
                for i=1 : ni
                   if ll(i) >=uu(i), igls=1; end
                end
                if ~igls
                   r=-MA.int.invM*MA.int.g;
                   U=qr([MA.int.M,r]);
                   R=triu(U(1:ni,1:ni));
                   y=U(1:ni,ni+1);
                   MA.int.p=obils_search(R,y,ll,uu);
                   if isempty(MA.int.p)||norm(MA.int.p)==0
                       ptr=xbest(info.iI)-XXi(info.iI,end);
                       if norm(prt)~=0
                          MA.int.p=MA.int.Deli*round(ptr/norm(ptr));
                       else
                           state='iMres'; % skip iTRS
                           break;
                       end
                       if isempty(MA.int.p)||norm(MA.int.p)==0
                           state='iMres'; % skip iTRS
                           break
                       end
                   end
                else 
                    state='iMres'; % skip iTRS
                    break
                end
                state='iTR3'; % go to compute a new feasible trial point
            end
            if strcmp(state,'iTR3') 
                info.ytrial    =xbest;
                info.statep='i'; info.beta=1;
                info.p=MA.int.p;
                [info]=projectStep(xbest,info);
                sxf =size(info.XF,1);
                diff=(info.XF(:,1:dim)-repmat(info.ytrial',sxf,1)).^2;
                [mn,~]    =min(sum(diff,2));
                if (mn>10^-16) % go to compute a new trial point
                   state='iTR4'; MA.int.changePTR=0;
                   x=info.ytrial; 
                   return; % return MATRSstep to compute f at x
                else
                    state='iTR5'; % the feasibe trial point already
                                % has been evaluated
                    if cputime-initTime>=secmax % secmax reached
                       x=''; % this is needed to avoid ifinity 
                       % loop when no feasible point exists
                       return
                    end
                end
            end
            if  strcmp(state,'iTR5') % update TR radius
                 MA.int.istuck=MA.int.istuck+1; MA.int.changePTR=1;
                 if sign(rand-0.5)
                    MA.int.Deli=...
                    abs(MA.int.Deli+sign(rand-0.5)*...
                                 randi([1,itune.Delmin],1));
                 else
                     MA.int.Deli=round(MA.int.Deli/...
                                      randi([1,itune.Delmin],1));
                 end
                 if MA.int.Deli<=0||MA.int.istuck>itune.stuckmax
                    state='iMres'; % skip iTRS
                    break; 
                 end
                 state='iTR2'; % go back to the next iteration of iTRS
            end
            if strcmp(state,'iTR4')   % save x and f in the history
               info=UpdatePoints(info,mitune);
               gisci =max(-MA.int.g'*MA.int.p,...
                             sqrt(eps)*abs(MA.int.g)'*abs(MA.int.p));
                if gisci==0
                    if fbest>info.ftrial, suciTR=0;
                    else, suciTR=1;
                    end
                else
                    mueff  =(fbest-info.ftrial)/gisci;
                    suciTR =(mueff*abs(mueff-1) > itune.izeta);
                end
                if info.iTRSdone
                   info.niTR=info.niTR+1; info.iTRSdone=0;
                end
                if suciTR % iteration is successful
                   fbest=info.ftrial; xbest=info.ytrial;  
                   if prt>=2 
                       fprintf(['iTRS: fbest=',num2str(fbest),'\n'])
                   end
                   MA.int.goodt=1;
                   MA.int.good=MA.int.good ||MA.int.goodt;
                   MA.int.succ_iTRS=MA.int.succ_iTRS+1; 
                   if MA.int.Deli<=itune.iDeltabar
                        MA.int.Deli=MA.int.Deli+1;
                    else
                        MA.int.Deli=...
                            floor(itune.itheta*max(norm(MA.int.p,inf),...
                                          MA.int.Deli));
                    end
                else % iteration is unsuccessful
                    MA.int.unsucc_iTRS=MA.int.unsucc_iTRS+1;
                    if MA.int.Deli<=itune.iDeltabar
                        MA.int.Deli=MA.int.Deli-1;
                    else
                        MA.int.Deli=...
                            floor(min(norm(MA.int.p,inf),...
                                          MA.int.Deli)/itune.itheta);
                    end
                end
                if MA.int.Deli<=0 % iTRS ends
                    state='iMres'; % end iTRS
                    break; 
                end
                % go back to form the trust-region subproblem and solve it
                state='iTR2'; 
            end
     end % end of iTRS
     if strcmp(state,'iMres')
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % change affine scaling matrix %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
       if norm(MA.int.M,inf)>=itune.mmax
         MA.int.s=zeros(ni, 1); MA.int.M=eye(ni); 
       end
       if nc>0 && ni==0, state='cMut1'; % go back to cMATRS
       elseif ni>0 && nc==0,state='iMut1';  % go back to iMATRS
       elseif ni>0 && nc>0, state='miLS1'; % go to miMATRS
       end
     end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% miMATRS: mixed-integer MATRS %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(state,'miLS1') % compute combination direction
        p0=xbest-info.xbest0; nF=length(info.F);
        if nF>=2
            for i=1:nF, dX(:,i)=info.X(:,i)-xbest;end
            sc=abs(max(dX')'); sc(sc==0)=1;
            if norm(sc)==0, sc=ones(dim,1); end
        else
            sc=ones(dim,1);
        end
        sc(info.iI)=max(1,1./sc(info.iI));
        sc(info.cI)=min(1,1./sc(info.cI));
        changDic=norm(p0)==0;
        if ~changDic 
           p0         =sc.*p0;
           MA.cont.p  =p0(info.cI); 
           MA.int.p   =ceil(p0(info.iI));
           MA.mixed.p =p0;
           MA.mixed.p(info.iI)=MA.int.p;
           if norm(MA.int.p)~=0
              % find largest allowed step sizes
              info.step='mi'; info.ls='i';
              MA.int.ixbest=xbest(info.iI);
              [MA,info]=getalpha(MA,info);
           end
           if norm(MA.cont.p)~=0
             % find largest allowed step sizes
             info.step='mi'; info.ls='c';
             MA.cont.cxbest=xbest(info.cI);
             [MA,info]=getalpha(MA,info);
           end
           if norm(MA.int.p)~=0||norm(MA.cont.p)~=0
             state='miLS2'; % go to mixed-integer phase
           else % skip mixed-integer phase and go back to cMATRS
             state='cMut1'; % skip miMATRS and go back to do cMutation
           end
        else
            MA.int.p=[]; MA.cont.p=[];
            for kk=1:length(info.F)-2
                p0       =xbest-info.X(:,end-kk);
               p0        =sc.*p0;
                MA.cont.p=p0(info.cI); 
                MA.int.p =ceil(p0(info.iI));
                if norm(MA.int.p)~=0 || norm(MA.cont.p)~=0,
                     MA.mixed.p =p0; MA.mixed.p(info.iI)=MA.int.p;
                    break;
                end
            end
            if (norm(MA.int.p)==0 && norm(MA.cont.p)==0) || ...
                    isempty(MA.int.p)&&isempty(MA.cont.p)
                % skip mixed-integer phase and go back to cMATRS
                state='cMut1'; % skip miMATRS and go back to do cMutation
            else
                 if norm(MA.int.p)~=0 || norm(MA.cont.p)~=0
                   if norm(MA.int.p)~=0
                      % find largest allowed step sizes
                      info.step='mi'; info.ls='i';
                      MA.int.ixbest=xbest(info.iI);
                      [MA,info]=getalpha(MA,info);
                   end
                   if norm(MA.cont.p)~=0
                     % find largest allowed step sizes
                     info.step='mi'; info.ls='c';
                     MA.cont.cxbest=xbest(info.cI);
                     [MA,info]=getalpha(MA,info);
                   end
                    state='miLS2'; % go to continue mixed-integer phase
                 else
                     % skip mixed-integer phase and go back to cMATRS
                     state='cMut1'; 
                 end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % miLSS: mixed-integer line search %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(state,'miLS2')
      info.ytrial =xbest;
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
           sxf =size(info.XF,1);
           diff=... 
           (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
           [mn,ind]=min(sum(diff,2));
           if (mn<=10^-16) 
                 % x is not a new feasible point
                fval=info.XF(ind(1),info.dim+1); info.ftrial=fval;
                state='miLS4'; % reuse know f from the history
           else % x is a new point
                state='miLS3'; % compute f at x
                return
           end
      else  % skip mixed-integer phase and go back to cMATRS
           state='cMut1';
       end
   end
   if strcmp(state,'miLS3') % check whether x is a new feasible point
       info.number_miMATRS=info.number_miMATRS+1;
        if prt>=1
           disp('===========================================')
           fprintf([num2str(info.number_miMATRS),'th mixed-integer\n'])
           disp('===========================================')
        end
        info=UpdatePoints(info,mitune);
        state='miLS4'; % go to do extrapolation
   end
   if strcmp(state,'miLS4')
        MA.mixed.extrapDone=0;
        if info.ftrial<fbest % extrapolation may be tried
            info.enf=info.nf-1;
            MA.mixed.extrapDone=1;
            if MA.int.feasible
                MA.cont.alpha=inf; MA.int.alpha=MA.int.alp0;
                % calculate trial point
                info.ytrial=xbest;
                info.beta =min(MA.int.alpmax,itune.inu*MA.int.alpha);
                % compute trial point and project it into [low upp]
                info.statep='i'; info.p=MA.int.p; 
                [info]=projectStep(xbest,info);
            end
            if MA.cont.feasible
                MA.cont.alpha=MA.cont.alp0; MA.int.alpha=inf;
                % calculate trial point
                info.ytrial=xbest;
                info.beta=min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                % compute trial point and project it into [low upp]
                info.statep='c'; info.p=MA.cont.p; 
                [info]=projectStep(xbest,info);
            end
            if MA.int.feasible||MA.cont.feasible
               x=info.ytrial;
               sxf =size(info.XF,1);
               diff=... 
               (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
               [mn,ind]=min(sum(diff,2));
               if (mn<=10^-16) 
                    % x is not a new feasible point
                    fval=info.XF(ind(1),info.dim+1); info.ftrial=fval;
                    state='miLS6'; % go to reuse known f from the history
               else % x is a new point
                    state='miLS5'; % go to compute f at x
                    return
               end
            else % no descrease in f
                state='miLS14'; % go to update information
            end
        else
            state='miLS8'; % go to try opposite direction
        end
   end
   if strcmp(state,'miLS5')  % check whether x is a new feasible point
        info=UpdatePoints(info,mitune);
        state='miLS6'; % go to do extrapolation
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%
   % extrapolation along p %
   %%%%%%%%%%%%%%%%%%%%%%%%%
   while strcmp(state,'miLS6')||strcmp(state,'miLS7') 
       % loop for extrapolation
        if strcmp(state,'miLS6')
              okLS=((MA.int.alpha<MA.int.alpmax || ...
                    MA.cont.alpha<MA.cont.alpmax) && ...
                    info.ftrial < fbest);
             if okLS
                info.ytrial =xbest;  
                MA.cont.alpha=...
                         min(MA.cont.alpmax, ctune.cnu*MA.cont.alpha); 
                if(MA.cont.alpha < MA.cont.alpmax)
                    info.beta=...
                         min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                    % compute trial point and project it into [low upp]
                    info.statep='c'; info.p=MA.cont.p; 
                    [info]=projectStep(xbest,info);
                end
                MA.int.alpha=...
                         min(MA.int.alpmax, itune.inu*MA.int.alpha);
                if (MA.int.alpha < MA.int.alpmax)  
                    info.beta =...
                          min(MA.int.alpmax,itune.inu*MA.int.alpha);
                    % project trial point into [low upp]
                    info.statep='i'; info.p=MA.int.p; 
                    [info]=projectStep(xbest,info);
                end
                if (MA.cont.alpha < MA.cont.alpmax)||...
                                         (MA.int.alpha < MA.int.alpmax) 
                   x=info.ytrial;
                   sxf =size(info.XF,1);
                   diff=... 
                   (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                   [mn,ind]=min(sum(diff,2));
                   if (mn<=10^-16)
                       % x is not a new feasible point
                       fval=info.XF(ind(1),info.dim+1); 
                       info.ftrial=fval;
                       state='miLS6'; % go to reuse f from the history
                   else % x is a new point
                       state='miLS7'; % go to compute f at x
                       return
                   end
                else
                    state='miLS14'; % go to update information
                    break;
                end
             else
                 state='miLS14'; % go to update information
                 break;
             end
        end
       if strcmp(state,'miLS7') % check whether x is a new feasible point
            info=UpdatePoints(info,mitune);
            state='miLS6'; % go to continue extrapolation
       end
   end % end of extrapolatio
   if strcmp(state,'miLS8') % opposite direction is tried
      MA.int.p=- MA.int.p; MA.cont.p=-MA.cont.p; 
      if norm(MA.int.p)~=0
          % find largest allowed step sizes
          info.step='mi'; info.ls='i';
          MA.int.ixbest=xbest(info.iI);
          [MA,info]=getalpha(MA,info);
       end
       if norm(MA.cont.p)~=0
         % find largest allowed step sizes
         info.step='mi'; info.ls='c';
         MA.cont.cxbest=xbest(info.cI);
         [MA,info]=getalpha(MA,info);
       end
       % compute the first triap point along opposite direction
      info.ytrial=xbest;
       if MA.int.feasible
             MA.int.alpha=0;
             % compute trial point and project it into [low upp]
             info.statep='i'; info.p=MA.int.p; info.beta=MA.int.alp0;
             [info]=projectStep(xbest,info);
       end
       if MA.cont.feasible
             MA.cont.alpha=0; 
             % compute trial point and project it into [low upp]
             info.statep='c'; info.p=MA.cont.p; info.beta=MA.cont.alp0;
             [info]=projectStep(xbest,info);
       end
       if MA.cont.feasible ||MA.int.feasible % trial point is feasible
           x=info.ytrial;
           sxf =size(info.XF,1);
           diff=... 
           (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
           [mn,ind]=min(sum(diff,2));
           if (mn<=10^-16)
                % x is not a new feasible point
                fval=info.XF(ind(1),info.dim+1); info.ftrial=fval;
                state='miLS10'; % go reuse f form the history
           else % x is a new point
                state='miLS9'; % go to compute f at x
                return
           end
       else
           state='cMut1'; % end miMATRS and go back to cMATRS
       end
   end
   if strcmp(state,'miLS9') % check whether x is a new point
        info=UpdatePoints(info,mitune);
        state='miLS10'; % go to do extrapolation
   end
   if strcmp(state,'miLS10')
       MA.mixed.extrapDone=0;
    if info.ftrial<fbest % the first trial point is the first 
                         % trial point of extrapolation
        MA.mixed.extrapDone=1;   
        info.enf=info.nf-1;
        % initialize alpha 
        info.ytrial=xbest;
       if MA.int.feasible
            MA.cont.alpha=inf; MA.int.alpha=MA.int.alp0;
            % compute step size
            info.beta     =min(MA.int.alpmax,itune.inu*MA.int.alpha);
            % compute trial point and project it into [low upp]
            info.statep='i'; info.p=MA.int.p; 
            [info]=projectStep(xbest,info);
       end
        if MA.cont.feasible
            MA.cont.alpha=MA.cont.alp0; MA.int.alpha=inf;
            % calculate step size
            info.beta=min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
            % project trial point into [low upp]
            info.statep='c'; info.p=MA.cont.p; 
            [info]=projectStep(xbest,info);
        end
        if MA.int.feasible||MA.cont.feasible
           x=info.ytrial;
           sxf =size(info.XF,1);
           diff=... 
           (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
           [mn,ind]=min(sum(diff,2));
           if (mn<=10^-16)
                % x is not a new feasible point
                fval=info.XF(ind(1),info.dim+1); info.ftrial=fval;
                state='miLS12'; % go to reuse f form the history
           else % x is a new point
                state='miLS11'; % go to compute f at x
                return
           end
        else
            state='miLS14'; % go to update information
        end
    else
        state='miLS14'; % go to update information
    end
   end
   if strcmp(state,'miLS11') % check whether x is a new feasible point
        info=UpdatePoints(info,mitune);
        state='miLS12'; % go to do extrapolation
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%
   % extrapolation along -p %
   %%%%%%%%%%%%%%%%%%%%%%%%%%
   while strcmp(state,'miLS12')||strcmp(state,'miLS13') 
       % loop for extrapolation
        if strcmp(state,'miLS12')
            % expansion step (increase stepsize)
            okLS=(MA.int.alpha<MA.int.alpmax || ...
                       MA.cont.alpha<MA.cont.alpmax) && ...
                                                info.ftrial < fbest;
            if okLS
                info.ytrial =xbest;  
                 MA.cont.alpha=...
                         min(MA.cont.alpmax, ctune.cnu*MA.cont.alpha); 
                if(MA.cont.alpha < MA.cont.alpmax)
                   
                    info.beta=...
                         min(MA.cont.alpmax,ctune.cnu*MA.cont.alpha);
                    % compute trial point and project it into [low upp]
                    info.statep='c'; info.p=MA.cont.p; 
                    [info]=projectStep(xbest,info);
                end
                MA.int.alpha=...
                         min(MA.int.alpmax, itune.inu*MA.int.alpha);
                if (MA.int.alpha < MA.int.alpmax)  
                   info.beta =...
                         min(MA.int.alpmax,itune.inu*MA.int.alpha);
                   % compute trial point and project it into [low upp]
                   info.statep='i'; info.p=MA.int.p; 
                   [info]=projectStep(xbest,info);
                end
                if (MA.cont.alpha < MA.cont.alpmax)...
                                      ||(MA.int.alpha < MA.int.alpmax) 
                       x=info.ytrial;
                       sxf =size(info.XF,1);
                       diff=... 
                       (info.XF(:,1:info.dim)-repmat(x',sxf,1)).^2;
                       [mn,ind]=min(sum(diff,2));
                       if (mn<=10^-16)
                          % x is not a new feasible point
                          fval=info.XF(ind(1),info.dim+1); 
                          info.ftrial=fval;
                          state='miLS12'; % go to reuse f from the history
                       else % x is a new point
                          state='miLS13'; % go to compute f at x
                          return
                       end 
                else
                   state='miLS14'; % go to update information
                   break; 
                end
            else
                state='miLS14';  % go to update information
                break;
            end
        end
        if strcmp(state,'miLS13') % check whether x is
                                  % a new feasible point
            info=UpdatePoints(info,mitune);
            state='miLS12'; % go back to do more extrapolation
        end
   end % end of extrapolation
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  update information after %
   %  the end of extrapolation %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if strcmp(state,'miLS14') 
     if MA.mixed.extrapDone % extrapolation 
           % choose a trial point with lowest 
           % f computed by extrapolation
           X=info.XF(info.enf:info.nf-1,1:dim)';
           F=info.XF(info.enf:info.nf-1,dim+1);
           [ftrial,ib]=min(F); ib=ib(1); ytrial=X(:,ib);
           if ftrial < fbest
               % save step size
               MA.mixed.ae=norm(xbest-ytrial)/norm(MA.mixed.p);
               % update the best point
               fbest=ftrial; xbest=ytrial;
           end
           if prt>=2 
              fprintf(['Mixed-integer: fbest=',num2str(fbest),'\n'])
           end
     end
     state='cMut1'; % go back to do cMutation in cMATRS
   end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%