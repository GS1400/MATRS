%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% projectStep.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projectStep projects a trial point into [low upp]
%
function [info]=projectStep(xbest,info)
switch info.statep
    case 'c'
      info.ytrial(info.cI) = max(info.clow,min(info.cupp,...
                    xbest(info.cI)+info.beta*info.p));
    case 'i'
      info.ytrial(info.iI) = max(info.ilow,min(info.iupp,...
                    xbest(info.iI)+info.beta*info.p));
    otherwise
        error('info.statep must be either c or i')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

