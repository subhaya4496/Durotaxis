function [N, Nhist, mu, dt, nT, sig, sizedif, r_cut, r_cutE, r_list, delta, Nsave, BoxL, BoxW, D_T, D_R] = Parameter_file(Pe)

%  Parameters

N = 10000;                   % Number of cells in simulation
mu = 0.3;                   % Poisson's Ratio
% alpha = 25;                  % Non-dimensionalized factor 'dipole term'
% v0 = 3;                     % self propulsion velocity
D_T = 1;                    % Translational Diffusion Constant
D_R = 1;                  % Rotational Diffusion Constant
nT = 2e6;
Nhist = 1e3;                % interval for count of histogram 
dt = 1e-3;                  % Time-step duration


sig = 1.0;                  % Size of cell
sizedif = 0.0;              % Difference in size for two cell types
delta = 0.0;                % Buffer length                                              
Nsave = Nhist;                % Number of steps after which a figure is saved

factor = 1;
BoxL = 12.0*(sig+delta);     % Box length chosen in SIGMA units
BoxW = BoxL*factor;         % Box width

MAXHSD = sig + sizedif;     % This is equal to EHSDBB

%% Cut off lengths
r_cut =  1.0*MAXHSD;        % Maximum cut-off LJ and LS  
r_cutE = 7.0*MAXHSD;        % Cut off - Elasticity                                              
r_list = 7.5*MAXHSD;        % Cut off - List

%% Prefactors
% con3 = alpha*(1+mu)*(3/(16*pi));     % force prefactor term
% con4 = -B*(1+mu)/(8*pi);          % torque prefactor term

