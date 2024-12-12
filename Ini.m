function [xpFold, ypFold, theFold] = Ini1D(N, BoxL,sig) 
    BoxLhalf = BoxL/2.0;
    ypFold = rand(1,N)*10*sig - BoxLhalf;
    xpFold = zeros(1,N);
    theFold = rand(1,N)*2*pi;
 