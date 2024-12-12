# Durotaxis
The codes are divided into 2 subgroups.
Subgroup 1 - Steady state behavior of motile cells
The codes for subgroup 1 are listed in folder 'Steady_State'
1.

Subgroup 2 - Escape time of motile cells assisted by diffusion
The codes for subgroup 2 are listed in folder 'Escape'
1. Main_sim.m - This is the main executable file for simulating multiple cell trajectories simultaneously
2. Parameter_file.m - List of parameters are declared here.
3. Ini.m - The initial location and the orientation of the cells are declared here.
4. BoxE.m - This is a function file where elastic interaction between the cells and the boundary is calculated.
5. bndry.m - This is a function file where boundary conditions are implemented.
6. Complete_1D.m - This is another executable file for simulating multiple cell trajectories in 1 dimension simultaneously.
7. DuroT.m - Different methods for calculating durotactic index are mentioned here. Escape times and rates are also calculated in this versatile file.
