function [xpFnew, ypFnew] = adsorb(xpFnew, ypFnew, N, BoxL, BoxW, sig)
BoxLhalf = BoxL/2.0;

BoxW_R = BoxW - sig;
Boxhalf_R = BoxW_R/2.0;
for icell = 1: N
% if x- direction has periodicity
  if (xpFnew(icell)  > +BoxLhalf)   
      xpFnew(icell) = xpFnew(icell) - BoxL;                
  elseif (xpFnew(icell)  < -BoxLhalf)   
      xpFnew(icell) = xpFnew(icell) + BoxL;
  end
% % Fixed Boundary in x-direction 
%   if (xpFnew(icell)  > +Boxhalf_R)   
%       xpFnew(icell) =  BoxL_R - xpFnew(icell);                
%   elseif (xpFnew(icell)  < -Boxhalf_R)   
%       xpFnew(icell) = - xpFnew(icell) - BoxL_R;                
%   end           
%  
  if (ypFnew(icell)  > +BoxLhalf)   
      ypFnew(icell) =  BoxLhalf;                
  elseif (ypFnew(icell)  < -BoxLhalf)   
      ypFnew(icell) = -BoxLhalf;                
  end      

end