for Pe = [0.1]
    dirOut = sprintf('DI_Pe_Num=%g',Pe);
    
%     dirOut = 'D:\Durotaxis_Project\Data\B=0.1\FMI_Pe_Num=_uni_6_sigma';
%     dirOut = 'D:\Altered_complete\check1D_constflux';
    % if exist(dirOut, 'dir')
    %     rmdir(dirOut, 's');
    % end
    % mkdir(dirOut);
    % cd(dirOut)
    % parpool('AttachedFiles',{'Parameter_1D.m','Box_E.m', 'bndry.m', 'Ini1D.m', 'parsave.m'})
    A = [1.667,2.5,5];
    for i = 1:length(A)
        tic
        Vmax = A(i);
        B = Vmax;
        [N, Nhist, dt, nT, sig, Nsave, BoxL, D_T, D_R, nu] = Parameter_1D;
        wall_dens = 0;                              
            % ALLOCATE MEMORY 
            % initialize INITIAL CONFIGURATION 
                t = 0.d0;                                                              
                savecount = 0;     
                histcount = 0;      
                nhist = nT/(Nhist);
                yhist = zeros(1,N*nhist);
                xhist = zeros(1,N*nhist);
                Yhist = zeros(1,N*nhist);
                thehist = zeros(1,N*nhist);
                [xpFold, ypFold, theFold] = Ini1D(N, BoxL,sig) ;
            % INITIALIZE WALL INTERACTIONS   
                Locate = 0;
                loc_inc = 1;
    
            % % The Wall Interaction
                    FyWall = zeros(1,N);
                    FTy = FyWall;
             % FIRST STEP - MOVE CELLS 
              xpFnew = xpFold ;
              xpU = xpFnew;
              ypFnew = ypFold ;
              ypU = ypFnew;
             % FIRST STEP - ENFORCE PERIODICITY
             [xpFnew, ypFnew] = bndry(xpFnew, ypFnew, N, BoxL);
    
             % UPDATE OLD VALUES FOR USE 
                ypFold = ypFnew;
                FTx = zeros(1,N);
            % START TIME STEPPER - use ADAPTIVE EXPLICIT EULER for time stepper
            tic
            for nTime = 1: nT
            % UPDATE LIST if needed 
             
             % DEFINE ARRAYS     
               % Wall Interactions
                  [FyWall, TzWall] = Box_E(Vmax, BoxL, ypFold, theFold, N, sig, nu, B);
    %                 FyWall = -Vmax*ones(1,N);
                 %  SUM TO GET TOTAL FORCE and torque
    
                 FTy = FyWall;
                 TTz = TzWall;
            % % DISPLACEMENT MADE 
                % Displacement due to Forces
                 FDX = dt.*(Pe*cos(theFold) + FTx);
                 FDY = dt.*(Pe*sin(theFold) + FTy);
                 TDZ = dt.*(TTz);
                % Displacement due to diffusion
                 disX = sqrt(2.0*D_T*dt).*randn(1,N); 
                 disY = sqrt(2.0*D_T*dt).*randn(1,N);
                 disT = sqrt(2.0*D_R*dt).*randn(1,N);
                 xpU = xpU + disX + FDX;
                 ypU = ypU + disY + FDY;
                 xpFnew = xpFold + disX + FDX;
                 ypFnew = ypU;
                 theFnew = mod(theFold + disT + TDZ, 2*pi);
    
            % ENFORCE WALL/FIXED BOUNDARY     
            [xpFnew, ypFnew] = bndry(xpFnew, ypFnew, N, BoxL);
    
            % COUNT FOR HISTOGRAM
             if (mod(nTime,Nhist) == 0)
                 for icell=1:N
                     xhist(histcount*N+icell) = xpFnew(icell);
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
                    xpFold = xpFnew;
                    theFold = theFnew;
    
            end
            toc
            save_file = sprintf('Vmax=%g.mat',Vmax);
%             matfile = fullfile(dirOut, save_file);
            save (save_file, 'Vmax', 'sig', 'B', 'Yhist', 'xhist', 'yhist','thehist', 'Pe', 'dt', 'N', 'nT', 'Nsave', 'Nhist', 'BoxL', 'histcount', 'nhist', '-v7.3');      %,
    %     end
        toc
    end
    % cd ../
end