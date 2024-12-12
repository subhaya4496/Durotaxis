function [d1, d2] = wall_dist(d1,d2 ,N ,xpFold ,ypFold , BoxL, BoxW)
BoxLhalf = BoxL/2.0;
BoxWhalf = BoxW/2.0;
for cell = 1:N
   if xpFold(cell)> 0
       d1(cell) =  BoxLhalf -xpFold(cell);
   else
       d1(cell) = -BoxLhalf -xpFold(cell);
   end
   if ypFold(cell)> 0
       d2(cell) =  BoxWhalf -ypFold(cell);
   else
       d2(cell) = -BoxWhalf -ypFold(cell);
   end
end