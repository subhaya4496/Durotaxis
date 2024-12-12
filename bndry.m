function [xpFnew, ypFnew] = bndry(xpFnew, ypFnew, N, BoxL)
BoxLhalf = BoxL/2.0;
for icell = 1: N
  if (ypFnew(icell)  > +BoxLhalf-1/2)   
      ypFnew(icell) =  BoxLhalf-1/2;                
  elseif (ypFnew(icell)  < -BoxLhalf+1/2)   
      ypFnew(icell) = -BoxLhalf+1/2;                
  end
  if (xpFnew(icell)  > +BoxLhalf)   
      xpFnew(icell) =  xpFnew(icell) - BoxL;                
  elseif (xpFnew(icell)  < -BoxLhalf)   
      xpFnew(icell) = BoxL + xpFnew(icell);                
  end
end