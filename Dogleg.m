%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dogleg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dogleg finds an approximate solution of the trust region subproblem
%
function s=Dogleg(Delta,s,G,g)
sTGs     = s'*G*s;
pred     = -g'*s-.5*sTGs;  % Model reduction for Newton Pt
gTGg     = g'*G*g;
if pred >0
   if norm(s)>Delta, sc  = - (g'*g/gTGg)*g;
       if norm(sc)>=Delta, s=-Delta*(g/norm(g,inf));
       else
          a = sc; b = s; c = a'*(b-a); nba = norm(b-a)^2;
          if c<=0
             t = (-c+sqrt(c^2+nba*(Delta^2-norm(a)^2)))/nba;
          else
             Dela = Delta^2-norm(a)^2; t = Dela/(c+sqrt(c^2+nba*Dela));
         end
         t = real(t); s = a+t*(b-a);
       end
   end
else
    if gTGg>5*eps, sc  = - (g'*g/gTGg)*g;
       if norm(sc) <= Delta,s= sc;
       else, s = Delta*sc/norm(sc,inf);
       end
    else, s = -Delta*g/norm(g,inf);
    end
end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
