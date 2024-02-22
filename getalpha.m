%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% getalpha.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% getalpha provides requirements for cLSS, iLSS, and miLSS
%
function [MA,info]=getalpha(MA,info,tune)
switch info.step
    case 'm' % mutation
        switch info.ls
            case 'c' % continuous search
                cI  = info.cI; clow  = info.clow; cupp  = info.cupp;
                nc = length(cI);
                % find MA.cont.alpmax, largest continuous step size
                ii=0; MA.cont.feasible = 0; MA.cont.alp0=0;
                MA.cont.alpmax=0;
                while 1
                   ii=ii+1;
                   ind = (MA.cont.MutInfo.p > 0 & MA.cont.cxbest<cupp);
                   aupp = min((cupp(ind)-MA.cont.cxbest(ind))...
                        ./MA.cont.MutInfo.p(ind));
                   if isempty(aupp), aupp=inf; end
                   ind  = (MA.cont.MutInfo.p < 0 & MA.cont.cxbest>clow);
                   alow  =min((clow(ind)-MA.cont.cxbest(ind))...
                           ./MA.cont.MutInfo.p(ind));
                   if isempty(alow),alow=inf; end
                   MA.cont.alpmax = min(alow,aupp);
                   if MA.cont.alpmax>0, MA.cont.feasible=1;break;end
                   if ii>=tune.ndd, break; end
                   MutInfo.sd  = randn(nc,1);
                   MA.cont.MutInfo.p = MA.cont.M*MutInfo.sd;
                end
                if MA.cont.feasible
                   amin  = sqrt(MA.cont.sigma*MA.cont.a(MA.cont.cm)); 
                   MA.cont.alp0 = min(MA.cont.alpmax,amin);
                end
            case 'i' % integer search
                iI    = info.iI; ilow  = info.ilow;iupp  = info.iupp;
                ni    = length(iI);
                % find MA.int.alpmax, largest integer step size
                ii=0; MA.int.feasible = 0; MA.int.alp0=0; MA.int.alpmax=0;
                while 1
                   ii=ii+1;
                   if norm(MA.int.MutInfo.p)~=0
                       ind  = (MA.int.MutInfo.p > 0&MA.int.ixbest<iupp);
                       aupp = min(( iupp(ind)-MA.int.ixbest(ind) )...
                              ./MA.int.MutInfo.p(ind));
                       if isempty(aupp), aupp=inf; end
                       ind  = (MA.int.MutInfo.p<0&MA.int.ixbest>ilow);
                       alow = min((ilow(ind)-MA.int.ixbest(ind))...
                               ./MA.int.MutInfo.p(ind));
                       if isempty(alow), alow=inf; end
                       MA.int.alpmax = max(1,floor(min(alow,aupp)));
                       if MA.int.alpmax>=1,MA.int.feasible = 1;break;end
                   end
                   if ii>=tune.ndd, break; end
                    psize  = round(2+log(ii));
                    MA.int.MutInfo.sd = ...
                               randi([-psize,psize],ni,1);
                    MA.int.MutInfo.p = ...
                               round(MA.int.M*MA.int.MutInfo.sd);
                end
                if MA.int.feasible
                    amin  = ...
                    max(1,round(sqrt(MA.int.sigma*MA.int.a(MA.int.im))));
                    MA.int.alp0 = min(MA.int.alpmax,amin); 
                end
        end
    case 'r' % recombination
          switch info.ls
              case 'c' % continuous search
                   % initialize vector alpha_max
                   MA.cont.feasible = 0;
                   ind = ...
                    (MA.cont.RecomInfo.p > 0&MA.cont.cxbest<info.cupp);
                   aupp = min((info.cupp(ind)-MA.cont.cxbest(ind))...
                          ./MA.cont.RecomInfo.p(ind));
                   if isempty(aupp), aupp=inf; end
                   ind = ...
                    (MA.cont.RecomInfo.p < 0&MA.cont.cxbest<info.clow);
                   alow = min((info.clow(ind)-MA.cont.cxbest(ind))...
                          ./MA.cont.RecomInfo.p(ind));
                   if isempty(alow), alow=inf; end
                   MA.cont.alpmax=min(alow,aupp);
                   if MA.cont.alpmax>0 
                      MA.cont.feasible = 1; 
                      amin = sqrt(MA.cont.sigma*MA.cont.ar); 
                      MA.cont.alp0 = min(MA.cont.alpmax,amin);
                   end
              case 'i' % integer search
                   % initialize vector alpha_max
                   MA.int.feasible = 0;
                   ind = ...
                     (MA.int.RecomInfo.p > 0 & MA.int.ixbest<info.iupp);
                   aupp = min((info.iupp(ind)-MA.int.ixbest(ind))...
                          ./MA.int.RecomInfo.p(ind));
                   if isempty(aupp), aupp=inf; end
                   ind = ...
                    (MA.int.RecomInfo.p < 0 & MA.int.ixbest<info.ilow);
                   alow = min((info.ilow(ind)-MA.int.ixbest(ind))...
                         ./MA.int.RecomInfo.p(ind));
                   if isempty(alow), alow=inf; end
                   MA.int.alpmax = max(1,floor(min(alow,aupp)));
                   if MA.int.alpmax>=1 
                      MA.int.feasible = 1;
                      amin  = ...
                      max(1,round(sqrt(MA.int.sigma*MA.int.ar)));
                      MA.int.alp0 = min(MA.int.alpmax,amin); 
                   end
          end
    case 'mi' % mixed-integer phase
       switch info.ls
            case 'i' % integer search 
                ilow  = info.ilow; iupp   = info.iupp;
                ixbest = MA.int.ixbest; MA.int.feasible = 0;
                MA.int.alp0=0;
                if norm(MA.int.p)~=0
                    ind = (MA.int.p > 0 & ixbest<iupp);
                    aupp= min((iupp(ind)-ixbest(ind))./ MA.int.p(ind));
                    if isempty(aupp), aupp=inf; end
                    ind = (MA.int.p < 0 & ixbest>ilow);
                    alow=min((ilow(ind)-ixbest(ind))./ MA.int.p(ind));
                    if isempty(alow),alow=inf; end
                    MA.int.alpmax = max(1,floor(min(alow,aupp)));
                    if MA.int.alpmax>=1
                      MA.int.feasible = 1;
                      amin  = ...
                      max(1,round(sqrt(MA.int.sigma*MA.mixed.ae)));
                      MA.int.alp0 = min(MA.int.alpmax,amin); 
                    end
                end
           case 'c' % continuous search
                clow = info.clow; cupp = info.cupp;
                MA.cont.feasible=0;
                cxbest = MA.cont.cxbest;MA.cont.alp0=0;
                if norm(MA.cont.p)~=0
                    ind = (MA.cont.p > 0 & cxbest<cupp);
                    aupp=min((cupp(ind)-cxbest(ind) )./MA.cont.p(ind));
                    if isempty(aupp), aupp=inf; end
                    ind = (MA.cont.p < 0 & cxbest>clow);
                    alow=min((clow(ind)-cxbest(ind))./MA.cont.p(ind));
                    if isempty(alow), alow=inf; end
                    MA.cont.alpmax=min(alow,aupp);
                    if MA.cont.alpmax>0
                        MA.cont.feasible = 1;
                        amin = sqrt(MA.cont.sigma*MA.mixed.ae); 
                        MA.cont.alp0 = min(MA.cont.alpmax,amin);
                    end
                end
       end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%