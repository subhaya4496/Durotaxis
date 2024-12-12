function [N, Nhist, dt, nT, sig, Nsave, BoxL, D_T, D_R, nu] = Parameter_1D
%  Parameters
nu = 0.3;                    % Poisson's Ratio
N = 500;                   % Number of cells in simulation
D_T = 1;                     % Translational Diffusion Constant
D_R = 1;                     % Rotational Diffusion Constant
nT = 2.4e5;
Nhist = 10;              % interval for count of histogram 
dt = 1e-3;                   % Time-step duration
sig = 1.0;                   % Size of cell
% sizedif = 0.0;               % Difference in size for two cell types                                           
Nsave = Nhist;               % Number of steps after which a figure is saved
BoxL = 30.0*sig;           % Box length chosen in SIGMA units