%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% requirMixedLSS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% requirMixedLSS provides requirements for line searches
% in the mixed-integer phase  
%
function [MA,finfo]= requirMixedLSS(xbest,MA,finfo)
cI     = finfo.cI;
iI     = finfo.iI;
clow   = finfo.clow;
cupp   = finfo.cupp;
ilow  = finfo.ilow;
iupp  = finfo.iupp;
% find MA.int.alpmax, largest integer step size
ixbest = xbest(iI); MA.int.feasible = 0; MA.int.alp0=0; MA.cont.alp0=0;
if norm(MA.int.p)~=0
    ind = (MA.int.p > 0 & ixbest<iupp);
    aupp= min(( iupp(ind) - ixbest(ind) )./ MA.int.p(ind));
    if isempty(aupp), aupp=inf; end
    ind = ( MA.int.p < 0 & ixbest>ilow);
    alow=min((ilow(ind) - ixbest(ind) )./ MA.int.p(ind));
    if isempty(alow), alow=inf; end
    MA.int.alpmax = max(1,floor(min(alow,aupp)));
    if isfinite(MA.int.alpmax)&& MA.int.alpmax>=1
        MA.int.feasible = 1;
    end
    if MA.int.feasible
        MA.int.alp0 = max(1,round(sqrt(MA.int.alpmax*MA.int.sigma)));
    end
end
% find MA.cont.alpmax, largest real step size
MA.cont.feasible=0; cxbest = xbest(cI);
if norm(MA.cont.p)~=0
ind = (MA.cont.p > 0 & cxbest<cupp);
aupp= min(( cupp(ind) - cxbest(ind) )./ MA.cont.p(ind));
if isempty(aupp), aupp=inf; end
ind = (MA.cont.p < 0 & cxbest>clow);
alow=min(( clow(ind) - cxbest(ind) )./ MA.cont.p(ind));
if isempty(alow), alow=inf; end
MA.cont.alpmax=min(alow,aupp);
if (MA.cont.alpmax>0 && ~isinf(MA.cont.alpmax))
    MA.cont.feasible = 1;
end
% find MA.cont.alp0, the real initial step size
if MA.cont.feasible
    MA.cont.alp0 = sqrt(MA.cont.alpmax*MA.cont.sigma);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

