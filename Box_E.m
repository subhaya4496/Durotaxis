% wall interaction
function [FyWall, TzWall] = Box_E(Vmax, BoxL, ypFold, theFold,  N, sig, nu, B)
%     FyWall = (Vmax*pi/BoxL)*(sin(2*pi*ypFold/BoxL));
    d2 = zeros(1,N);
    FyWall = zeros(1,N);
    TzWall = zeros(1,N);
    % Clamped boundary
    a_nu = -(1+nu)*(15+32*nu*(nu-1));
    b_nu = -(1+nu)*(34+32*(nu^2) - 72*nu);
    c_nu = -(1+nu)*(7-8*nu);
    abc_d = (1-nu)*(3-4*nu);
    for cell = 1:N
        if ypFold(cell)> 0
           d2(cell) =  BoxL/2 -ypFold(cell);
        else
           d2(cell) = -BoxL/2 -ypFold(cell);
        end
        FyWall(cell)= 3*Vmax/(abs(d2(cell)) + sig)^(4)*sign(d2(cell))*(1.0/abc_d)*(a_nu + b_nu*(sin(theFold(cell)))^2 + c_nu*(sin(theFold(cell)))^4);
        FyWall(cell)= FyWall(cell)*abc_d/(a_nu + b_nu/2 + 3*c_nu/8);
        TzWall(cell)= B*(1/(abs(d2(cell))+sig))^(3)*(1.0/abc_d)*(b_nu*sin(pi- 2*theFold(cell)) +c_nu*(0.5*sin(-4*theFold(cell))+sin(pi-2*theFold(cell))));
        TzWall(cell)= TzWall(cell)*abc_d/(b_nu + c_nu);
    end