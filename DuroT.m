%% Durotactic Index Calculation

DI =[];
% DI_2 = [];
T = [];
DIT = [];
for Vmax = 1:20
    tic
    filename = sprintf('Vmax=%g_B=0.1',Vmax);
    load(filename,'Yhist','N','BoxL', 'Pe')
    sigma = 1;
    Boxhalf = BoxL/2;
    NPos = 0;
    NNeg = 0;
    L_start = 6*sigma;
    Time = [];
    FMI_time = [];
    NP = 0;
    NN = 0;
%     L_right = 4*sigma;
    for step = 1:length(Yhist)-N
        if Yhist(step+N)-Yhist(step)<0
            NPos = NPos + 1;
        elseif Yhist(step+N)-Yhist(step)>0
            NNeg = NNeg + 1;
        end
        if mod(step, N*100)==0
            Time = [Time, step/(N*10)];
            FMI_time = [FMI_time,(NPos-NNeg)/(NPos+NNeg)];
            NP = NP + NPos;
            NN = NN + NNeg;
            NPos = 0;
            NNeg = 0;
        end
    end
    DuroIndex = (NP-NN)/(NP+NN);
    DI = [DI, DuroIndex];
%     end
    
    toc
    T = [T; Time];
    DIT = [DIT;FMI_time];
%     Nleft = sum(Yhist(1:50*N)<3*sigma);
%     Nright = sum(Yhist(1:50*N)>3*sigma);
%     D2 = (Nleft - Nright)/(Nleft + Nright);
% 
%     DI_2 = [DI_2, D2];
end
savefile = sprintf('Pe=%g.mat', Pe);
save(savefile,'DI', 'T', 'DIT')

%% Durotactic Index Calculation

DI =[];
% DI_2 = [];
T = [];
DIT = [];
DIIT = [];
NT = [];
for A = [0, 0.5, 1, 1.5, 1.667, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10]
    tic
    filename = sprintf('Vmax=%g',A);
    load(filename,'Yhist','N','BoxL', 'Pe', 'sig','dt')
    Boxhalf = BoxL/2;
    NPos = 0;
    NNeg = 0;
    NPoss = 0;
    NNegg = 0;
%     L_start = 6*sig;
    Time = [];
    DI_time = [];
    DII_time = [];
    N_time = [];
    NP = 0;
    NN = 0;
    L_right = 2*sig;
    L_left = 0*sig;
    nstep = 1;
    for step = 1:(length(Yhist)/(N*nstep))-1
        for particle = 1:N
            if (Yhist((step-1)*N*nstep + particle) < L_right) && (Yhist((step-1)*N*nstep + particle) > L_left) 
                if Yhist(step*N*nstep + particle) <=  Yhist((step-1)*N*nstep + particle)
                    NPos = NPos + 1;
                    NPoss = NPos*abs(Yhist(step*N*nstep + particle) - Yhist((step-1)*N*nstep + particle));
                elseif Yhist(step*N*nstep + particle) >  Yhist((step-1)*N*nstep + particle)
                    NNeg = NNeg + 1;
                    NNegg = NNeg*abs(Yhist(step*N*nstep + particle) - Yhist((step-1)*N*nstep + particle));
                end
%             elseif (Yhist((step-1)*N*nstep + particle) < L_left && (Yhist(step*N*nstep + particle) > L_left))
%                 NNeg = NNeg + 1;
            elseif (Yhist((step-1)*N*nstep + particle) > L_right && (Yhist(step*N*nstep + particle) < L_right))
                NPos = NPos + 1;
                NPoss = NPos*abs(Yhist(step*N*nstep + particle) - Yhist((step-1)*N*nstep + particle));
            end
        end
        DII = (NPoss-NNegg)/(NPoss+NNegg);
        Time = [Time, step*nstep*dt*10];
        DII_time = [DII_time, DII];
        DI_time = [DI_time,(NPos-NNeg)/(NPos+NNeg)];
        N_time = [N_time, NPos + NNeg];
        NP = NP + NPos;
        NN = NN + NNeg;
        NPoss = 0;
        NNegg = 0;
        NPos = 0;
        NNeg = 0;
    end
    DuroIndex = (NP-NN)/(NP+NN);
    DI = [DI, DuroIndex];
%     end
    
    toc
    T = [T; Time];
    NT = [NT; N_time];
    DIT = [DIT;DI_time];
    DIIT = [DIIT;DII_time];
%     Nleft = sum(Yhist(1:50*N)<3*sigma);
%     Nright = sum(Yhist(1:50*N)>3*sigma);
%     D2 = (Nleft - Nright)/(Nleft + Nright);
% 
%     DI_2 = [DI_2, D2];
end
T = 1:24;
Durr = [mean(DIIT(:,1:1000),2),mean(DIIT(:,1001:2000),2),mean(DIIT(:,2001:3000),2),mean(DIIT(:,3001:4000),2),mean(DIIT(:,4001:5000),2),mean(DIIT(:,5001:6000),2),mean(DIIT(:,6001:7000),2),mean(DIIT(:,7001:8000),2),mean(DIIT(:,8001:9000),2),mean(DIIT(:,9001:10000),2),mean(DIIT(:,10001:11000),2),mean(DIIT(:,11001:12000),2),mean(DIIT(:,12001:13000),2),mean(DIIT(:,13001:14000),2),mean(DIIT(:,14001:15000),2),mean(DIIT(:,15001:16000),2),mean(DIIT(:,16001:17000),2),mean(DIIT(:,17001:18000),2),mean(DIIT(:,18001:19000),2),mean(DIIT(:,19001:20000),2),mean(DIIT(:,20001:21000),2),mean(DIIT(:,21001:22000),2),mean(DIIT(:,22001:23000),2),mean(DIIT(:,23001:end),2)];
dddd = mean(Durr,2);
DI = mean(DIT(:,20:end),2);
savefile = sprintf('DI_ZZZ_Pe=%g_2.5_0-2.mat', Pe);
save(savefile,'DI', 'T', 'DIT', 'DIIT', 'NT','Durr','dddd')



%%
bins = 20;
A = 10;
h = histogram(ET(A,:),bins);
val = h.Values;
maxET = max(ET(A,:));
time = maxET*(1:bins)/bins - (1/2/bins);
t = 0:0.01:maxET;
v_spline = csaps(time,val,0.0000001,t);
figure
semilogy(time, val, 'LineStyle','none','Marker','*', 'Color', '[0.8500 0.3250 0.0980]', 'MarkerSize', 9 ,'DisplayName','T_{esc} distribution')
hold on
plot(t,v_spline, 'LineWidth', 2.5, 'Color', '[0.8500 0.3250 0.0980]','DisplayName', 'Fit')
plot(t, 2250*exp(-t/2350))
% gca.xlabel = 'T_{esc}';
% gca.ylabel = 'Number of escapes';
% hold off


%%
figure
semilogy(Pe, A1, 'LineStyle','None', 'Marker','*', 'MarkerSize', 9,'LineWidth',2, 'DisplayName','A = 1')
hold on
semilogy(Pe, A2, 'LineStyle','None', 'Marker','diamond', 'MarkerSize', 9,'LineWidth',2,'DisplayName','A = 2')
hold on
semilogy(Pe, A5, 'LineStyle','None', 'Marker','pentagram', 'MarkerSize', 9,'LineWidth',2, 'DisplayName','A = 5')
hold on
semilogy(Pe, A10, 'LineStyle','None', 'Marker','hexagram', 'MarkerSize', 9,'LineWidth',2, 'DisplayName','A = 10')
hold off
%% FMI calculation
T = [];
FMI =[];
FMIT = [];

A = [0, 0.5, 1, 1.5, 1.667, 2, 2.5, 3, 3.5, 4, 4.5, 5];
FMID = NaN(length(A), 500*10);
for k = 1:length(A)
    tic
    Vmax = A(k);
    filename = sprintf('Vmax=%g',Vmax);
    load(filename,'xhist','Yhist','N','BoxL', 'Pe', 'sig','dt')
    Boxhalf = BoxL/2;
    Time = [];
    FMI_time = [];
    Dx = zeros(1,N);
    Dy = zeros(1,N);
    L = zeros(1,N);
    nstep = 100;
    Tlim = 1;
    L_right = 10*sig;
    L_left = 0*sig;
    F_dist = [];
    for step = 1:100
        for particle = 1:N
            if (Yhist(Tlim*1000*N + particle -N) < L_right) % (Yhist((Tlim-1)*1000*N*nstep + particle) < L_right) && 
                dx = xhist(step*N*nstep + particle) - xhist((step-1)*N*nstep + particle);
                Dx(particle) = Dx(particle) + dx;
                dy = - Yhist(step*N*nstep + particle) + Yhist((step-1)*N*nstep + particle);
                if dy~=0
                    Dy(particle) = Dy(particle) + dy;
                    L(particle) = L(particle) + sqrt(dx^2 + dy^2);
%                 else
%                     Dy(particle) = Dy(particle) + dx;
%                     L(particle) = L(particle) + sqrt(dx^2 + dy^2);
                
                end
            end
        end
        if mod((step+1)*nstep, 1000)==0
            Tlim = Tlim+1;
            F = Dy./L; 
%             F(F==0)=1;
            FMI_time = [FMI_time, mean(F,'omitnan')];
            Time = [Time, step*nstep*dt*10];
            Dx = zeros(1,N);
            Dy = zeros(1,N);
            L = zeros(1,N);
            F_dist = [F_dist, F];

        end
    end
    T = [T;Time];
    FMIT = [FMIT; FMI_time];
    ave_FMI = mean(FMI_time,'omitnan');
    FMI = [FMI, ave_FMI];
    FMID(k,:) = F_dist;
    toc
end

% STDEV1 = std(FMID(1,:),'omitnan');
% STDEV4 = std(FMID(4,:),'omitnan');
% STDEV6 = std(FMID(6,:),'omitnan');
% STDEV11 = std(FMID(11,:),'omitnan');
% Mean1 = mean(FMID(1,:),'omitnan');
% Mean4 = mean(FMID(4,:),'omitnan');
% Mean6 = mean(FMID(6,:),'omitnan');
% Mean11 = mean(FMID(11,:),'omitnan');
savefile = sprintf('FMI_Pe=%g_0-10_10_maybe.mat', Pe);
save(savefile,'A', 'FMI', 'Time', 'FMI_time', 'FMIT', 'T', 'FMID')

% FMID_pos1 = FMID(FMID(1,:)>0);
% FMID_pos4 = FMID(FMID(4,:)>0);
% FMID_pos6 = FMID(FMID(6,:)>0);
% FMID_pos11 = FMID(FMID(11,:)>0);
% STDEV1 = std(FMID_pos1,'omitnan');
% STDEV4 = std(FMID_pos4,'omitnan');
% STDEV6 = std(FMID_pos6,'omitnan');
% STDEV11 = std(FMID_pos11,'omitnan');
% Mean1 = mean(FMID_pos1,'omitnan');
% Mean4 = mean(FMID_pos4,'omitnan');
% Mean6 = mean(FMID_pos6,'omitnan');
% Mean11 = mean(FMID_pos11,'omitnan');
% errorbar([A(1),A(4), A(6), A(11)], [Mean1, Mean4, Mean6, Mean11], [STDEV1/sqrt(length(FMID)), STDEV4/sqrt(length(FMID)), STDEV6/sqrt(length(FMID)), STDEV11/sqrt(length(FMID))])

%% Probability at the boundary
AI_all = [];
DI_all = [];
N_all = [];
PB_all = [];
for Pec = [0, 0.1, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10]% 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    direc = sprintf('DI_Pe_Num=%g',Pec);
    cd (direc)
    P_bound_A = [];
    N_bound_A = [];
    A = [0, 0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10];
    for V = 1:length(A)
        filename = sprintf('Vmax=%g',A(V));
        load(filename,'Yhist','N','BoxL', 'Pe', 'sig','dt')
        Y_steady = Yhist(length(Yhist)*2/5:end);
        Y_bound = Y_steady(Y_steady<=sig/2);
        P_bound = length(Y_bound)/length(Y_steady);
        P_bound_A = [P_bound_A,P_bound];
        if A(V) == 0
            Nb0 = length(Y_bound);
        end
        Nb = length(Y_bound);
        N_bound_A = [N_bound_A, Nb];
    end
    DI = (N_bound_A - Nb0)./(N_bound_A+Nb0);
    AI = (N_bound_A + Nb0)./(2*length(Y_steady));
    AI_all = [AI_all;AI];
    DI_all = [DI_all;DI];
    N_all = [N_all; N_bound_A];
    PB_all = [PB_all; P_bound_A];
    % plot(A,DI)
    % hold on
    % plot(A, (P_bound_A-P_bound_A(1))/(P_bound_A(end)-P_bound_A(1)))
    % hold on
    % [po,Durtactic_100] = find(DI>0.377);
    cd ../
end
%% Accumulation index
% X_duro = A;
X_duro = repmat(A,[8,1]);
Duro = [DI_all(1,:);DI_all(9,:);DI_all(10,:);DI_all(11,:);DI_all(12,:);DI_all(13,:);DI_all(14,:),DI_all(15,:)];
Duro(Duro<=DI_all(1,2)) = NaN;
X_duro(isnan(Duro)) = NaN;
i = 1;
% figure
for Pe = [0, 1, 2, 3, 4, 5, 10, 15]
    on = ones(1,length(A));
    xx = X_duro(i,:);
    plot(xx, Pe*on,'LineStyle','None','Marker', 'pentagram','MarkerSize',9,'MarkerEdgeColor','[0.4660 0.6740 0.1880]','MarkerFaceColor','[0.4660 0.6740 0.1880]')
    i = i+1;
    hold on
end

%%
X_duro = A;
X_duro(isnan(aduro)) = NaN;
i=1;
for Pe = 0.1
    on = ones(1,length(A));
    xx = X_duro(i,:);
    semilogy(xx, Pe*on,'LineStyle','None','Marker', '>','MarkerSize',9,'MarkerEdgeColor','[0.4940 0.1840 0.5560]','MarkerFaceColor','[0.4940 0.1840 0.5560]')
    i = i+1;
    hold on
end
%% FMI calculation averaged
T = [];
FMI =[];
FMI_time = [];
FMIT_std = [];
FMI_num = [];

A = [1.667, 2.5, 5];
for k = 1:length(A)
    tic
    Vmax = A(k);
    filename = sprintf('Vmax=%g',Vmax);
    load(filename,'xhist','Yhist','N','BoxL', 'Pe', 'sig','dt')
    Boxhalf = BoxL/2;
    Time = [];
    FMI_T_N = NaN(240,N);
    Dx = 0;
    Dy = 0;
    L = 0;
    nstep = 250;
    
    L_right = 6*sig;
    L_left = 0*sig;
    F_dist = [];
    
    for particle = 1:N
        Tlim = 1;
        step = 1;
        check = 1;
        while Yhist((Tlim-1)*4*nstep*N + particle) > 0 && step < length(Yhist)/(N*nstep) %
            flag = (check-1)*N*nstep + particle;
            if (Yhist(flag) < L_right)||(Yhist(Tlim*4*nstep*N + particle -N) < L_right)
                dx = xhist(step*N*nstep + particle) - xhist((step-1)*N*nstep + particle);
                Dx = Dx + dx;
                dy = - Yhist(step*N*nstep + particle) + Yhist((step-1)*N*nstep + particle);
                if dy~=0
                    Dy = Dy + dy;
                    L = L + sqrt(dx^2 + dy^2);
    %                 else
    %                     Dy(particle) = Dy(particle) + dx;
    %                     L(particle) = L(particle) + sqrt(dx^2 + dy^2);
                
                end
            end
            if mod(step, 4)==0
                F = Dy/L;
                FMI_T_N(Tlim,particle) = F;
                Dx = 0;
                Dy = 0;
                L = 0;
                Tlim = Tlim+1;
                check = step + 1;
            end
            step = step + 1;
        end
    end
    toc
    FMIT = mean(FMI_T_N,2,'omitnan');
    FMIT_S = std(FMI_T_N,0,2,'omitnan');
    FMI_N = sum(~isnan(FMI_T_N),2);
    FMIT_std = [FMIT_std, FMIT_S];
    FMI_time = [FMI_time, FMIT];
    FMI = [FMI,FMIT(1,:)];
    FMI_num = [FMI_num,FMI_N];
end
savefile = sprintf('FMI_Pe=%g_0-6_250_maybe.mat', Pe);
save(savefile,'A', 'FMI', 'FMIT', 'FMI_time','FMIT_std','FMI_num')
%% Escape time
ET = [];
MET = [];
S = [];
for V = 1:10
    et = [];
    file = sprintf('Vmax=%g.mat',V);
    load(file, 'Yhist', 'N', 'dt', 'Nhist')
    for i = 1:N
        for j = 1:length(Yhist)/N
            t = j*dt*Nhist;
            if Yhist(i + N*(j-1))>6
                break
            % elseif t == length(Yhist/N)*dt*Nhist
            %     break
            end
        end
        et = [et,t];
    end
    ET = [ET;et];
    MET = [MET,mean(et,'omitnan')];
    S = [S,std(et,'omitnan')];
end
%% Orientation vs Time
Pe = 1;
cntA = 1;
Theta_all = zeros(9, length(the_hist));
for A = [200, 500, 800, 1000, 1200, 1500, 1800, 2000, 2200]
    openfilename = sprintf('Clamped_Pe_Num_=%g_A=%g.mat', Pe, A);
    load(openfilename)
    Theta = cos(the_hist + pi/2);
    Theta_all(cntA,:) = Theta;
    T = 1:length(the_hist)/N;
    cntA = cntA + 1;
end

%% Flip Time
flips = zeros(1, 9);
FTC = zeros(N,9);
% A = [200,500,800,1000,1200,1500,1800,2000,2200];
for j = 1:9
    % load('Clamped_Pe_Num_=1_A=%g.mat',A(j))
    calc = zeros(N, length(the_hist)/N - 1);
    for cell = 1:N
        for i = 1:width(calc)
            calc(cell,i) = Theta_all(j, cell+(N*(i-1)))*Theta_all(j, cell+(N*i));
        end
        ctr = calc(calc(cell,:) < 0);
        FTcell = length(ctr);
        FTC(cell,j) = FTcell;
    end
    
end

% 
%     flips(j) = length(ctr);
% flip_freq = flips/length(the_hist);