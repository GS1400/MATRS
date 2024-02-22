%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRSstep.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update the affine scaling matrix
%
function [MA]=updateM(MA,info,tune)
if strcmp(info.var,'c') % in continuous searches
    MA.cont.s = (1-MA.cont.c_s)*MA.cont.s + MA.cont.sqrt_s*MA.cont.sumz;
    MA.cont.M = (1 - 0.5*MA.cont.c_1 - 0.5*MA.cont.c_mu) * MA.cont.M + ...
              (0.5*MA.cont.c_1)*(MA.cont.M*MA.cont.s)*MA.cont.s';
    for m = 1:MA.cont.mu
     u1 = MA.cont.MutInfoList(MA.cont.permut(m)).p; 
     u2 = MA.cont.MutInfoList(MA.cont.permut(m)).sd; 
     MA.cont.M = MA.cont.M + ((0.5*MA.cont.c_mu*MA.cont.wi(m))*u1)*u2';
    end  
    MA.cont.M=adjustY(MA.cont.M,tune);
else % in integer searches
    MA.int.s = (1-MA.int.c_s)*MA.int.s + MA.int.sqrt_s*MA.int.sumz;
    MA.int.M = (1 - 0.5*MA.int.c_1 - 0.5*MA.int.c_mu) * MA.int.M + ...
              (0.5*MA.int.c_1)*(MA.int.M*MA.int.s)*MA.int.s';
    for m = 1:MA.int.mu
     u1 = MA.int.MutInfoList(MA.int.permut(m)).p; 
     u2 = MA.int.MutInfoList(MA.int.permut(m)).sd; 
     MA.int.M = MA.int.M + ((0.5*MA.int.c_mu*MA.int.wi(m))*u1)*u2';
    end  
    MA.int.M=adjustY(MA.int.M,tune);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

