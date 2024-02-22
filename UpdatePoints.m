%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% UpdatePoints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UpdatePoints save m last points and their function values and all 
% evaluated trial points and their function value
%
function info=UpdatePoints(info,tune)
% create and update the first list which uses to approximate gradient of
% the objective model function of trust-region subproblems and check
% whether or not the evaluated feasible point is a new feasible point
info.ftrial  = info.f; 
info.XF(info.nf,1:info.dim) = info.ytrial'; 
info.XF(info.nf,info.dim+1) = info.ftrial;
info.nf=info.nf+1;
% creat and update the second list which use to compute combination
% directions
if info.nc>0 && info.ni>0 
    if size(info.X,2)>=tune.m
       iw     = floor(tune.m/2)+randperm(floor(tune.m/2),1);
       xs0    = info.X; fs0=info.F;
       info.X = [xs0(:,1:iw-1) xs0(:,iw+1:end)];
       info.F = [fs0(:,1:iw-1) fs0(:,iw+1:end)];
    end
    xs0 = info.X; fs0=info.F;
    I1 = find(info.F<info.ftrial);
    I2 = find(info.F>info.ftrial);
    if ~isempty(I1) && ~isempty(I2)
       info.X = [xs0(:,I1) info.ytrial xs0(:,I2)];
       info.F = [fs0(I1) info.ftrial fs0(I2)];
    elseif ~isempty(I1)
       info.X = [xs0(:,I1) info.ytrial];
       info.F = [fs0(I1) info.ftrial];
    elseif ~isempty(I2)
       info.X = [info.ytrial xs0(:,I2)];
       info.F = [info.ftrial fs0(I2)];
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%