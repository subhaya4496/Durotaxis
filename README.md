# Durotaxis
The codes are divided into 2 subgroups.
Subgroup 1 - Steady state behavior of motile cells
The codes for subgroup 1 are listed in folder 'Steady_State'
1. ES2b_T.m -> Executible file - Simulation with fixed/reflecting boundary conditions with and without wall interaction.
2. Parameter_file.m -> Parameter File - Called in the executable file.
3. wall_E.m -> Called in the executable file ES2b_T.m, calculates the boundary interactions (Free boundary/ Clamped boundary)
4. diluteB.m -> Defining Cell Position and Orientation - Called in the executable file. Random or lattice arrangement of cells. Also controls the Box size or BoxL.
5. adsorb.m -> Defining boundary conditions - Called in the executable file. Confinement in one direction and periodic in the orthogonal direction.
6. wall_dist.m -> Called in the executable file. Distance from the boundary is estimated here

Post-processing of images in ImageJ for movies.

Subgroup 2 - Escape time of motile cells assisted by diffusion
The codes for subgroup 2 are listed in folder 'Escape'
1. Main_sim.m - This is the main executable file for simulating multiple cell trajectories simultaneously
2. Parameter_file.m - List of parameters are declared here.
3. Ini.m - The initial location and the orientation of the cells are declared here.
4. BoxE.m - This is a function file where elastic interaction between the cells and the boundary is calculated.
5. bndry.m - This is a function file where boundary conditions are implemented.
6. Complete_1D.m - This is another executable file for simulating multiple cell trajectories in 1 dimension simultaneously.
7. DuroT.m - Different methods for calculating durotactic index are mentioned here. Escape times and rates are also calculated in this versatile file.
