%% wall interaction
% function [TzWall] = wall_E(d1, d2, theFold, con4, N, r_cutE, mu, sig, k_wall)
% % FxWall = zeros(1,N);
% % FyWall = zeros(1,N);
% TzWall = zeros(1,N);
% 
% % % Free boundary
% % a_nu = (5+2*mu*(6*mu-1));
% % b_nu = (22+4*mu*(2*mu-9));
% % c_nu = (13*(1-2*mu)+12*mu^2);
% % abc_d = (1-mu);
% 
% % Clamped boundary
% a_nu = -(15+32*mu*(mu-1));
% b_nu = -(34+32*mu^2-72*mu);
% c_nu = -(7-8*mu);
% abc_d = (1-mu)*(3-4*mu);
% 
% for cell = 1: N
% %     if (d1(cell)< r_cutE)
% %         FxWall(cell)=(con3/16.0)*(1/abs(d1(cell)))^5*(1.0/abc_d)*(a_nu + b_nu*(cos(theFold(cell)))^2 + c_nu*(cos(theFold(cell)))^4)*d1(cell);
% %         TzWall(cell)= TzWall(cell)- (con4/32.0)*(1/abs(d1(cell)))^3*(1.0/abc_d)*(b_nu*sin(2*theFold(cell)) +c_nu*(0.5*sin(4*theFold(cell))+sin(2*theFold(cell))));
% %     end
%     
% %     FyWall(cell)= -(con3/16.0)*(1/(abs(d2(cell))+sig))^4*(1.0/abc_d)*(a_nu + b_nu*0.5 + c_nu*0.375)*sign(d2(cell));
%     TzWall(cell)= TzWall(cell)+ (con4/32.0)*(1/(abs(d2(cell))+sig))^3*(1.0/abc_d)*(b_nu*sin(pi- 2*theFold(cell)) +c_nu*(0.5*sin(-4*theFold(cell))+sin(pi-2*theFold(cell))));
% 
% end

%% Double well potential
function [FyWall] = wall_E(ypFold, BoxL, Vmax)
    FyWall= (Vmax*pi/BoxL)*sin(2*pi*ypFold/BoxL);