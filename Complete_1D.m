Tflip = 1;                 % Corresponding to B values in the flip time simulations
for Pe = 0
    dirOut = sprintf('Pe_Num=%g',Pe);
    % if exist(dirOut, 'dir')
    %     rmdir(dirOut, 's');
    % end
    % mkdir(dirOut);
    cd(dirOut)
    A = 1:2;
    for i = 1:length(A)
        Vmax = A(i);
        N = 500;                   % Number of cells in simulation
        D_T = 0.5;                     % Translational Diffusion Constant
        nT = 1e6;
        Nhist = 1000;              % interval for count of histogram 
        dt = 1e-3;                   % Time-step duration
        sig = 1.0;                   % Size of cell
        % sizedif = 0.0;               % Difference in size for two cell types                                           
        Nsave = Nhist;               % Number of steps after which a figure is saved
        BoxL = 30.0*sig;           % Box length chosen in SIGMA units
        Pflip = dt/Tflip;
        wall_dens = 0;                              
        t = 0.d0;                                                              
        savecount = 0;     
        histcount = 0;      
        nhist = nT/(Nhist);
        yhist = zeros(1,N*nhist);
        Yhist = zeros(1,N*nhist);
        thehist = zeros(1,N*nhist);
        BoxLhalf = BoxL/2.0;
        ypFold = rand(1,N)*0 - BoxLhalf;
        theFold = randi([0 1],1,N)*2*pi;
            % INITIALIZE WALL INTERACTIONS   
                Locate = 0;
                loc_inc = 1;
    
            % % The Wall Interaction
                    FyWall = zeros(1,N);
                    FTy = FyWall;
             % FIRST STEP - MOVE CELLS 
              ypFnew = ypFold ;
              ypU = ypFnew;
             % FIRST STEP - ENFORCE PERIODICITY
        BoxLhalf = BoxL/2.0;
        for icell = 1: N
          if (ypFnew(icell)  > +BoxLhalf)   
              ypFnew(icell) =  BoxLhalf;                
          elseif (ypFnew(icell)  < -BoxLhalf)   
              ypFnew(icell) = -BoxLhalf;                
          end
        end
             % UPDATE OLD VALUES FOR USE 
                ypFold = ypFnew;
            % START TIME STEPPER - use ADAPTIVE EXPLICIT EULER for time stepper
            tic
            for nTime = 1: nT
            % UPDATE LIST if needed 
             
             % DEFINE ARRAYS
             d2 = zeros(1,N);
             FyWall = zeros(1,N);
             for cell = 1:N
                 if ypFold(cell)> 0
                    d2(cell) =  BoxL/2 -ypFold(cell);
                 else
                    d2(cell) = -BoxL/2 -ypFold(cell);
                 end
                 FyWall(cell)= 3*Vmax*(abs(d2(cell)))/((d2(cell)^2) + sig^2)^(5/2)*sign(d2(cell));
             end
             FTy = FyWall;
%                  TTz = TzWall;
            % % DISPLACEMENT MADE 
                % Displacement due to Forces
                 FDY = dt.*(Pe*cos(theFold) + FTy);
%                  TDZ = dt.*(TTz);
                % Displacement due to diffusion
                 chh = rand(1,N);
                 flip = zeros(1,N);
                 flip(chh<=Pflip) = 1;
                 flip(chh>Pflip) = 0;
                 disY = sqrt(2.0*D_T*dt).*randn(1,N);
                 % disT = sqrt(2.0*D_R*dt).*randn(1,N);
                 ypU = ypU + disY + FDY;
                 ypFnew = ypU;
                 theFnew = mod(theFold + flip*pi, 2*pi);
    
            % ENFORCE WALL/FIXED BOUNDARY     
            for icell = 1: N
              if (ypFnew(icell)  > +BoxLhalf)   
                  ypFnew(icell) =  BoxLhalf;                
              elseif (ypFnew(icell)  < -BoxLhalf)   
                  ypFnew(icell) = -BoxLhalf;                
              end
            end
    
            % COUNT FOR HISTOGRAM
             if (mod(nTime,Nhist) == 0)
                 for icell=1:N
                     yhist(histcount*N+icell) = ypFnew(icell);
                     thehist(histcount*N+icell) = theFnew(icell);
                     Yhist(histcount*N+icell) = BoxL/2 - abs(ypFnew(icell));
    %                              if (ypFnew(icell)~=ypU(icell))
    %                                  wall_dens = wall_dens + 1;
    %                                  Locate(loc_inc) = histcount*N+icell;
    %                                  loc_inc = loc_inc+1;
    %                              end
                 end
                 histcount = histcount +1;
             end
    %           end
    
              % UPDATE CURRENT TIME
    
                t = t + dt;
             % UPDATE VARIABLES FOR NEXT RUN 
                    ypU = ypFnew;
                    ypFold = ypFnew;
                    theFold = theFnew;
    
            end
            toc
            save_file = sprintf('Vmax=%g.mat',Vmax);
%             matfile = fullfile(dirOut, save_file);
            save (save_file, 'Vmax', 'sig', 'Yhist', 'yhist','thehist', 'Pe', 'dt', 'N', 'nT', 'Nsave', 'Nhist', 'BoxL', 'histcount', 'nhist', '-v7.3');      %,
    %     end
    end
    cd ../
end