%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% reqcontLSS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% reqcontLSS provides requirements for cLSS
%
function [MA,finfo]=reqcontLSS(MA,finfo)

cI     = finfo.cI;
clow  = finfo.clow;
cupp  = finfo.cupp;
nc = length(cI);
% find MA.cont.alpmax, largest continuous step size
ii=0; MA.cont.feasible = 0; MA.cont.alp0=0;MA.cont.alpmax=0;
while 1
   ii=ii+1;
   ind = (MA.cont.Offspring.p > 0 & MA.cont.cxbest<cupp);
   aupp = min((cupp(ind)-MA.cont.cxbest(ind))...
        ./MA.cont.Offspring.p(ind));
   if isempty(aupp), aupp=inf; end
   ind  = (MA.cont.Offspring.p < 0 & MA.cont.cxbest>clow);
   alow  =min((clow(ind)-MA.cont.cxbest(ind))...
           ./MA.cont.Offspring.p(ind));
   if isempty(alow),alow=inf; end
   MA.cont.alpmax = min(alow,aupp);
   if isfinite(MA.cont.alpmax)&&MA.cont.alpmax>0
       MA.cont.feasible=1;
       break;
   end
   if ii>=1, break; end
   Offspring.sd        = randn(nc,1);
   MA.cont.Offspring.p = MA.cont.M*Offspring.sd;
end
if MA.cont.feasible
    amin  = sqrt(MA.cont.sigma*MA.cont.a(MA.cont.cm)); 
    MA.cont.alp0 = min(MA.cont.alpmax,amin);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%