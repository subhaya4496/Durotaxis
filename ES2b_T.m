   % FIXED BOUNDARY PROBLEM

   % Program BD_N with Verlet neighbour lists - Aug30, 2021 revised
   %  Bidisperse spheres - Types A and B
   %  A has diameter eq 1 (diameter sigma = 1)  and velocity v0A
   %  VEL(N) AND RANCELL(N) - reflect this difference
   
   % Periodic boundary conditions in 2D,
   % Boxlength of {-L/2, L/2} in simulation units centered at origin 
   
   % (xpU, ypU, TpU) are unfolded coordinates
   % (xpF, ypF, TpF) are  folded coordinates
   
   % Time scaled with (3/D_r)
   
   % SET RANDOM NUMBER before using PARFOR
 
   % LOAD SAVED CONFIGURATION IF NEEDED

   % load ('./*.mat')

   % DEFINE PARAMETERS 
   
   %  cellrunID = 1;
 
   %  fileIDF = fopen('cellF1.dat','w');
   %  fileIDU = fopen('cellU1.dat','w');
   
for Pe = [0]
% opendir = sprintf ('Pe_Num=%g_DT=1',Pe);
%    mkdir (opendir)
%    cd (opendir)
for Vmax = [0.1,0.2,0.5,1, 2, 5]
%     B = 10000;
    [N, Nhist, mu, dt, nT, sig, sizedif, r_cut, r_cutE, r_list, delta, Nsave, BoxL, BoxW, D_T, D_R] = Parameter_file(Pe);
    vel = zeros(1,N);                                          
    wall_dens = 0;                              
    k_wall = 1/dt;

        % ALLOCATE MEMORY 
        % initialize INITIAL CONFIGURATION 

            t = 0.d0;                                                              
            savecount = 0;     
            histcount = 0;      
            nhist = nT/(Nhist);

            xhist = zeros(1,N*nhist);
            Xhist = zeros(1,N*nhist);
            yhist = zeros(1,N*nhist);
            Yhist = zeros(1,N*nhist);
            the_hist = zeros(1,N*nhist);

        %     xhist = zeros(1,N);     
%             yhist = zeros(1,N);

            [xpFold, ypFold,theFold] = diluteB(N, delta, sig, BoxL, BoxW);   % 
        %      xpFold = zeros(1,N);
        %      ypFold = zeros(1,N);
             %      Area = BoxL*BoxL;  

        % %  Density of cells

        %     50 % is CELL A and 50% is CELL B 
        %     Thus the formulae:
        %     AcellA = (N/2.0)*pi*(rsigA)*(rsigA);
        %     AcellB = (N/2.0)*pi*(rsigB)*(rsigB);
        %     densA = AcellA/Area;
        %     densB = AcellB/Area;
        %     density = densA + densB;

        %     xp1 = xpFold;
        %     yp1 = ypFold;
        %     thep1 = theFold;
        %     
        %     tp1 = 0.0;
        % INITIALIZE WALL INTERACTIONS   

            d1 = zeros(1,N);
            d2 = zeros(1,N);
            [d1, d2] = wall_dist(d1,d2 ,N ,xpFold ,ypFold , BoxL, BoxW);
            Locate = 0;
            loc_inc = 1;
        % CREATE LISTS FOR FIRST USE 
        %     [POINT, LIST, xpFsaved, ypFsaved, theFsaved, timesaved] = LISTUB2(r_list, N, xpFold, ypFold, theFold, BoxL, t);

        % INITIALIZE TOTAL FORCES AND TORQUES

        %  % First the Linear Soring Force

        % [LSFTx, LSFTy, LSTTz] = ForceLSB(POINT, LIST, r_cut,  N, xpFold, ypFold, BoxL, sig, dt);

                   LSFTx = zeros(1,N);
                   LSFTy = zeros(1,N);
                   LSTTz = zeros(1,N);

        %  % The Elastic forces

        %   [EFTx, EFTy, ETTz] = ForceEB(POINT, LIST, r_cutE,  N, xpFold, ypFold, theFold, BoxL, con3, con4, mu, radcell);

                   EFTx = zeros(1,N);
                   EFTy = zeros(1,N);
                   ETTz = zeros(1,N);

        % % The Wall Interaction

        %   [FxWall, FyWall, TzWall] = wall_E(d1, d2, theFold, con3, con4, N, r_cutE, mu);

                   FxWall = zeros(1,N);
                   FyWall = zeros(1,N);
                   TzWall = zeros(1,N);

                  % Sum to get total force

%                 FTx = LSFTx + EFTx + FxWall;
%                 FTy = LSFTy + EFTy + FyWall;
%                 TTz = LSTTz + ETTz + TzWall;

         % GENERATE trial displacement for FIRST TIME STEP

%          disX = (sqrt(2.0)*B_T*sqrt(dt)).*randn(1,N);
%          disY = (sqrt(2.0)*B_T*sqrt(dt)).*randn(1,N);
%          disT = (sqrt(2.0)*B_R*sqrt(dt)).*randn(1,N);

         % FIRST STEP - MOVE CELLS 

          xpFnew = xpFold ; 
          ypFnew = ypFold ; 
          theFnew = theFold ;


          xpU = xpFnew;         % Unfolded coordinates
          XH = xpU;
          ypU = ypFnew;         % Unfolded coordinates
          theU = theFnew;       % Unfolded coordinates

         % FIRST STEP - ENFORCE PERIODICITY
         [xpFnew, ypFnew] = adsorb(xpFnew, ypFnew, N, BoxL, BoxW, sig);

         % UPDATE OLD VALUES FOR USE 

            xpFold = xpFnew;
            ypFold = ypFnew;
            theFold = theFnew;

            %    fprintf(fileIDF, '%d %f %f %f\n', N, B_T, B_R, mu);
            %    fprintf(fileIDF, '%f %f %f %f\n',v0,alpha,delta,BoxL);

            %    fprintf(fileIDU, '%d %f %f %f\n', N, B_T, B_R, mu);
            %    fprintf(fileIDU, '%f %f %f %f\n',v0,alpha,delta,BoxL);

        % START TIME STEPPER - use ADAPTIVE EXPLICIT EULER for time stepper
        tic
        for nTime = 1: nT
        %     B_R = (B_R + (1 - B_R)*(exp(-mod(nTime,50000)/1000)));         %  Required for annealing
        %     B_T = (B_T + (1 - B_T)*(exp(-mod(nTime,50000)/1000)));         %  Required for annealing
            % UPDATE LIST if needed 
%             [d1, d2] = wall_dist(d1,d2 ,N ,xpFold ,ypFold , BoxL, BoxW);


        %     [POINT, LIST, xpFsaved, ypFsaved, theFsaved, timesaved] = LISTUB2(r_list, N, xpFold, ypFold, theFold, BoxL, t);
         % GENERATE RANDOM NUMBERS for noise 

             Xran = sqrt(2.0*D_T)*randn(1,N);
             Yran = sqrt(2.0*D_T)*randn(1,N);
             Tran = sqrt(2.0*D_R)*randn(1,N);

         % DEFINE ARRAYS
          % CALCULATE TOTAL FORCE USING LISTS BASED ON OLD positions

           % % First the Linear Spring force

        %       [LSFTx, LSFTy, LSTTz] = ForceLSB(POINT, LIST, r_cut, N, xpFold, ypFold, BoxL, sig, dt);

           % % Then  Elastic forces

%               [EFTx, EFTy, ETTz] = ForceEB(POINT, LIST, r_cutE,  N, xpFold, ypFold, theFold, BoxL, con3, con4, mu, radcell);

           % % Wall Interactions

%               [TzWall] = wall_E(d1, d2, theFold, con4, N, r_cutE, mu, sig, k_wall);
                [FyWall] = wall_E(ypFold, BoxL, Vmax);

             %  SUM TO GET TOTAL FORCE and torque

%                 FTx = LSFTx + EFTx + FxWall;
                FTy = LSFTy + EFTy + FyWall;
%                  TTz = TzWall;
        %         FTx = 0.0;
        %         FTy = 0.0;
        %         TTz = 0.0;
        % % DISPLACEMENT MADE 
            % Displacement due to Forces     
             FDX = dt.*(Pe*cos(theFold)); 
             FDY = dt.*(FTy + (Pe*sin(theFold))); 
%              TDT = dt * TTz;     
            % Displacement due to diffusion 
             disX =sqrt(dt).*Xran;
             disY =sqrt(dt).*Yran;   
             disT =sqrt(dt).*Tran;

             xpU = xpU + disX + FDX; 
             ypU = ypU + disY + FDY; 
             theU = theU + disT;
             XH = XH + disX + FDX;

            xpFnew = xpU; 
            ypFnew = ypU;
            theFnew = mod(theFold + disT, 2*pi);

        % ENFORCE WALL/FIXED BOUNDARY     
        [xpFnew, ypFnew] = adsorb(xpFnew, ypFnew, N, BoxL, BoxW, sig);

        % COUNT FOR HISTOGRAM

%           if (nTime>(nT/2))
             if (mod(nTime,Nhist) == 0)
                 for icell=1:N
                     xhist(histcount*N+icell) = xpFnew(icell);
                     Xhist(histcount*N+icell) = XH(icell);
                     yhist(histcount*N+icell) = ypFnew(icell);
                     Yhist(histcount*N+icell) = BoxL/2 - abs(ypFnew(icell));
                     the_hist(histcount*N+icell) = theFnew(icell);
                     if (ypFnew(icell)~=ypU(icell))
                         wall_dens = wall_dens + 1;
                         Locate(loc_inc) = histcount*N+icell;
                         loc_inc = loc_inc+1;
                     end
                 end
                 histcount = histcount +1;
             end
%           end

          % UPDATE CURRENT TIME

            t = t + dt;

          % save data, fig, mat and also binned velocity data to get averages
          % to do the binning move back to actual bins using BoxL  

%          if (mod(nTime,Nsave) == 0)
%            
%            savecount = savecount + 1;    
%            savedataB(t, xpFnew, ypFnew, theFnew, v0, D_R, D_T, LSFTx, LSFTy, LSTTz, EFTx, EFTy, ETTz, FTx, FTy, TTz, alpha, mu, N, savecount, BoxL, BoxW, radcell, vel);
%         
%           % reset values for next save and print
%          
%         %         xp1 = xp2;
%         %         yp1 = yp2;
%         %         tp1 = tp2;
%         
%          end

         % UPDATE VARIABLES FOR NEXT RUN 

                
                xpU = xpFnew;
                ypU = ypFnew;
                xpFold = xpFnew;
                ypFold = ypFnew;
                theFold = theFnew;     

        end
        % Count for Historam for Time Evolution
%             for icell=1:N
%         %              xhist(icell) = xpFnew(icell);
%                      yhist(icell) = ypFnew(icell);
%             end
        toc
         %    fclose(fileIDF);
         %    fclose(fileIDU);
        filename = sprintf('Vmax=%g.fig', Vmax);
        %% Create histogram
        f = figure('visible','off');
        xlim(axes,[-21 21]);
        h = histogram(yhist, 1000);
    %     hmax(sim) = max(h.Values);
    %     hbulk(sim) = mean(h.Values(1:150));
        set(gcf,'Position',[100,100,1000,900]);
        title('Probability distribution of the cells in steady state','fontweight','bold','fontsize',20)
        xlabel('Distance from the wall','fontweight','bold','fontsize',16) 
        ylabel('Probability','fontweight','bold','fontsize',16)
        set(f, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
        saveas(f, filename)
        % Saving Data

    %     end
        save_param = sprintf('Vmax=%g.mat', Vmax);        
        save (save_param, 'xhist', 'Xhist', 'Yhist', 'yhist', 'the_hist', 'Vmax', 'dt', 'N', 'nT', 'Nsave', 'Nhist', 'BoxL', 'histcount', 'nhist', 'wall_dens', 'Locate');
%     end
%     cd ..\
end

end
