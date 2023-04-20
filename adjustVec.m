%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adjustVec.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [step]=adjustVec(step,tune)
% remove NaN or inf from a vector

function [vec]=adjustVec(vec,tune)

for i=1:length(vec)
    Inan  = isnan(vec);
    if any(Inan), vec(Inan) = (1+rand)*tune.gammav;end
    Iinf= (vec>=tune.gammav);
    if any(Iinf), vec(Iinf) = (1+rand)*tune.gammav;end
    Iinf= (vec<=-tune.gammav);
    if any(Iinf), vec(Iinf) = -(1+rand)*tune.gammav;end  
end  

Iz = (vec==0);
if any(Iz)
   maxsc=max(vec);
   vec(Iz)=maxsc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%