%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% reqintLSS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% reqintLSS provides requirements for iLSS
%  
function [MA,finfo]= reqintLSS(MA,finfo)
iI    = finfo.iI;
ilow  = finfo.ilow;
iupp  = finfo.iupp;
ni    = length(iI);
% find MA.int.alpmax, largest integer step size
ii=0; MA.int.feasible = 0; MA.int.alp0=0; MA.int.alpmax=0;
while 1
   ii=ii+1;
   if norm(MA.int.Offspring.p)~=0
       ind   = (MA.int.Offspring.p > 0 & MA.int.ixbest<iupp);
       aupp  = min(( iupp(ind)-MA.int.ixbest(ind) )...
              ./MA.int.Offspring.p(ind));
       if isempty(aupp), aupp=inf; end
       ind   = (MA.int.Offspring.p < 0 & MA.int.ixbest>ilow);
       alow  = min((ilow(ind)-MA.int.ixbest(ind))...
               ./MA.int.Offspring.p(ind));
       if isempty(alow), alow=inf; end
       MA.int.alpmax = max(1,floor(min(alow,aupp)));
        if isfinite(MA.int.alpmax)&& MA.int.alpmax>=1
            MA.int.feasible = 1;
           break; 
        end
   end
   if ii>=100, break; end
    psize  = round(2+log(ii));
    MA.int.Offspring.sd = randi([-psize,psize],ni,1);
    MA.int.Offspring.p  = round(MA.int.M*MA.int.Offspring.sd);
end
if MA.int.feasible
    amin  = max(1,round(sqrt(MA.int.sigma*MA.int.a(MA.int.im))));
    MA.int.alp0 = min(MA.int.alpmax,amin); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




