
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adjustX.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point]=adjustX(point,tune);
% remove NaN or inf from X
%

function [X]=adjustY(X,tune)

for i=1:size(X,1)
    Xi        = X(i,:);
    Inan      = isnan(Xi);
    if any(Inan), Xi(Inan) = (1+rand)*tune.gammaX; X(i,:)=Xi; end
    Iinf= (Xi>=tune.gammaX);
    if any(Iinf), Xi(Iinf) = (1+rand)*tune.gammaX; X(i,:)=Xi; end
    Iinf= (Xi<=-tune.gammaX);
    if any(Iinf), Xi(Iinf) = -(1+rand)*tune.gammaX; X(i,:)=Xi; end  
    if X(i,i)==0, X(i,i)=1; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

