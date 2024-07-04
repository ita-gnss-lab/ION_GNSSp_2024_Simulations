%% Comments
% Project Title: EKF_KF_VAR_MC_alt

% Author: Rodrigo de Lima Florindo

% Date: 22/06/2024

% Description: This code may be used to simulate the behaviour of
% multi-frequency KF-VAR and EKF-VAR carrier phase tracking loops during
% ionospheric scintillation through Monte Carlo Runs for different
% scintillation scenarios (Mild, Moderate and Severe).

% Its output are:
% The RMSE time series of all states of interests as 3 data sets for
% different ionopsheric scintillation scenarios and different Monte Carlo
% runs considering a fixed MFPSM realization and varying AWGN seeds.

clc;
clear all;
addpath([cd,'\MFPSM_Data']);
addpath([cd,'\Auxiliary_Codes']);
addpath([cd,'\Auxiliary_Codes\ARFIT']);

%% GENERAL SETUP
%Transmitted signal amplitudes
AmpL1 = 1;
AmpL2 = 1;
AmpL5 = 1;

%Transmitted signal frequencies
fcL1 = 154*10.23e6;
fcL2 = 120*10.23e6;
fcL5 = 115*10.23e6;

T_I = 0.01; % Integration period
time = 0:T_I:300; % Simulation Time
f_receiver_Hz = 20e6; % Receiver sampling frequency

% AWGN Parameters
CN0dBHzL1 = 42;
CN0dBHzL2 = 40;
CN0dBHzL5 = 44;
CN0L1 = 10^(CN0dBHzL1/10);
CN0L2 = 10^(CN0dBHzL2/10);
CN0L5 = 10^(CN0dBHzL5/10);

%Estimation of the variance of the noise from the output of the atan discriminators
sigma2Measure_valueL1 = (1/(2*CN0L1*T_I))*(1+(1/(2*CN0L1*T_I)));
sigma2Measure_valueL2 = (1/(2*CN0L2*T_I))*(1+(1/(2*CN0L2*T_I)));
sigma2Measure_valueL5 = (1/(2*CN0L5*T_I))*(1+(1/(2*CN0L5*T_I)));

% LOS paramaters
fd = 1000; % Frequency Shift
fdr = 0.94; % Frequency Drift
phiL1_0 = 0; % Initial phase shift for L1
phiL2_0 = 0; % Initial phase shift for L2
phiL5_0 = 0; % Initial phase shift for L5

M = 3; % Wiener process order

% Adaptative Module Configuration
% 1 - Non Adaptative
% 2 - Adaptative with noise
AdaptSwitch = 2;

%% MC Runs
MC = 1;
sigma2 = 2.6*10^-6;
TSTART_full = tic;
for s = 1:3
    TSTART_scen = tic;
    disp(['->Scenario: ',num2str(s)]);
    ScintScenario = s;
    if ScintScenario == 1
        load('S4_0p35_tau_2p0_200_L1_L2_L5_300s.mat');
    elseif ScintScenario == 2
        load('S4_0p57_tau_1p4_200_L1_L2_L5_300s.mat');
    elseif ScintScenario == 3
        load('S4_0p8_tau_0p8_200_L1_L2_L5_300s.mat');
    end

    Results = cell(1,7);

    RMSED_KF_Array = cell(MC,3);
    RMSEI_KF_Array = cell(MC,3);
    RMSET_KF_Array = cell(MC,3);

    RMSED_EKF_Array = cell(MC,3);
    RMSEI_EKF_Array = cell(MC,3);
    RMSET_EKF_Array = cell(MC,3);
    RMSERho_EKF_Array = cell(MC,3);
    for Topology = 1:3
        TSTART_topology = tic;
        disp(['--->Topology: ',num2str(Topology)]);
        TopologySelector = Topology;
        for i = 1:MC
            disp(['====>Monte Carlo Run: ',num2str(i)]);
            seedAWGN = i+10;
            scintSeed = 16;
            % KF-VAR - SPECIFIC MODEL CONFIGURATIONS AND SIMULATIONS
            run('KF_VAR_Specific_Config.m');

            % Received signal realization
            [yk,phi_D_true,rho_I_true,phi_I_true] = ReceivedSignal(TopologySelector,time,Scin_psi,AmpL1,AmpL2,AmpL5,T_I,f_receiver_Hz,fd,fdr,phiL1_0,phiL2_0,phiL5_0,seedAWGN,CN0dBHzL1,CN0dBHzL2,CN0dBHzL5,fcL1,fcL2,fcL5);
            
%             KF_tic = tic;
            [Px_0_0,x_hat_0_0,Qk,Fk,C] = KF_VAR_General_Config(L,M,NumSeries,P,T_I,sigma2,fd,fdr,deltaArray,Q_arfit,A,c,TopologySelector);
            [xkk1_KF, Pkk1_KF] = KF_VAR_PLL(Px_0_0, x_hat_0_0, yk, Fk, Qk, C, Rk, Hk, AdaptSwitch, T_I,CN0L1,CN0L2,CN0L5,TopologySelector);
            [~,~,~,ED_KF,EI_KF,ET_KF] = Evaluation_KFVAR(xkk1_KF,phi_D_true,phi_I_true,TopologySelector,L,M);
%             KF_toc = toc(KF_tic);
            
            RMSE_TS_ED_KF = zeros(length(time),Topology);
            RMSE_TS_EI_KF = zeros(length(time),Topology);
            RMSE_TS_ET_KF = zeros(length(time),Topology);
            SED_KF = ED_KF(:,:).^2;
            SEI_KF = EI_KF(:,:).^2;
            SET_KF = wrapToPi(ET_KF(:,:)).^2;
            for p = 1:length(time)
                RMSE_TS_ED_KF(p,:) = sqrt(mean(SED_KF(1:p,:),1));
                RMSE_TS_EI_KF(p,:) = sqrt(mean(SEI_KF(1:p,:),1));
                RMSE_TS_ET_KF(p,:) = sqrt(mean(SET_KF(1:p,:),1));
            end

            % EKF-VAR - SPECIFIC MODEL CONFIGURATIONS AND SIMULATIONS
            run('EKF_VAR_Specific_Config.m');
            
%             EKF_tic = tic;
            [Px_0_0,x_hat_0_0,Qk,Fk,C] = EKF_VAR_General_Config(L,M,NumSeries,P,T_I,sigma2,fd,fdr,deltaArray,Q_arfit,A,c,TopologySelector);
            [xkk1_EKF, Pkk1_EKF] = EKF_VAR_PLL(Px_0_0,x_hat_0_0,yk,Rk,Fk,Qk,M,L,NumSeries,P,AmpL1,AmpL2,AmpL5,Aux,C,TopologySelector,AdaptSwitch,CN0L1,CN0L2,CN0L5,T_I);
            [~,~,~,~,ED_EKF,EI_EKF,ET_EKF,ERho_EKF] = Evaluation_EKFVAR(xkk1_EKF,phi_D_true,phi_I_true,rho_I_true,TopologySelector,L,M);
%             EKF_toc = toc(EKF_tic);
            
%             disp(['ET - KF-VAR: ', num2str(KF_toc), '; ET - EKF-VAR: ', num2str(EKF_toc)]);
            
            RMSE_TS_ED_EKF = zeros(length(time),Topology);
            RMSE_TS_EI_EKF = zeros(length(time),Topology);
            RMSE_TS_ET_EKF = zeros(length(time),Topology);
            RMSE_TS_ERho_EKF = zeros(length(time),Topology);
            
            SED_EKF = ED_EKF(:,:).^2;
            SEI_EKF = EI_EKF(:,:).^2;
            SET_EKF = wrapToPi(ET_EKF(:,:)).^2;
            SERho_EKF = ERho_EKF(:,:).^2;
            for p = 1:length(time)
                RMSE_TS_ED_EKF(p,:) = sqrt(mean(SED_EKF(1:p,:),1));
                RMSE_TS_EI_EKF(p,:) = sqrt(mean(SEI_EKF(1:p,:),1));
                RMSE_TS_ET_EKF(p,:) = sqrt(mean(SET_EKF(1:p,:),1));
                RMSE_TS_ERho_EKF(p,:) = sqrt(mean(SERho_EKF(1:p,:),1));
            end

            % Cell assignment
            RMSED_KF_Array{i,Topology} = RMSE_TS_ED_KF;
            RMSEI_KF_Array{i,Topology} = RMSE_TS_EI_KF;
            RMSET_KF_Array{i,Topology} = RMSE_TS_ET_KF;
            
            RMSED_EKF_Array{i,Topology} = RMSE_TS_ED_EKF;
            RMSEI_EKF_Array{i,Topology} = RMSE_TS_EI_EKF;
            RMSET_EKF_Array{i,Topology} = RMSE_TS_ET_EKF;
            RMSERho_EKF_Array{i,Topology} = RMSE_TS_ERho_EKF;
        end
        TSTOP_topology = toc(TSTART_topology);
        disp(['Elapsed Time for this topology set: ', num2str(TSTOP_topology)]);
    end
    Results{1,1} = RMSED_KF_Array;
    Results{1,2} = RMSEI_KF_Array;
    Results{1,3} = RMSET_KF_Array;

    Results{1,4} = RMSED_EKF_Array;
    Results{1,5} = RMSEI_EKF_Array;
    Results{1,6} = RMSET_EKF_Array;
    Results{1,7} = RMSERho_EKF_Array;
    if s == 1
        save([cd,'\Results\Results_Mild_RMSE_TS_v2_1'], 'Results');
    elseif s == 2
        save([cd,'\Results\Results_Moderate_RMSE_TS_v2_1'], 'Results');
    elseif s == 3
        save([cd,'\Results\Results_Severe_RMSE_TS_v2_1'], 'Results');
    end
    TSTOP_scen = toc(TSTART_scen);
    disp(['Elapsed Time for this scenario: ', num2str(TSTOP_scen)]);
end
TSTOP_full = toc(TSTART_full);
disp(['Elapsed Time for all simulations: ', num2str(TSTOP_full)]);
