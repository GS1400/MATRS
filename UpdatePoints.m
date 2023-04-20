%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% UpdatePoints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UpdatePoints save 100 last points and their function values
%
function finfo=UpdatePoints(finfo)
if size(finfo.X,2)>=100
     iw  = 50+randperm(50,1);
     xs0 = finfo.X; fs0=finfo.F;
     finfo.X = [xs0(:,1:iw-1) xs0(:,iw+1:end)];
     finfo.F = [fs0(:,1:iw-1) fs0(:,iw+1:end)];
end
xs0 = finfo.X; fs0=finfo.F;
I1 = find(finfo.F<finfo.fworst);
I2 = find(finfo.F>finfo.fworst);
if ~isempty(I1) && ~isempty(I2)
   finfo.X  = [xs0(:,I1) finfo.yworst xs0(:,I2)];
   finfo.F  = [fs0(I1) finfo.fworst fs0(I2)];
elseif ~isempty(I1)
   finfo.X  = [xs0(:,I1) finfo.yworst];
   finfo.F  = [fs0(I1) finfo.fworst];
elseif ~isempty(I2)
   finfo.X  = [finfo.yworst xs0(:,I2)];
   finfo.F  = [finfo.fworst fs0(I2)];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


