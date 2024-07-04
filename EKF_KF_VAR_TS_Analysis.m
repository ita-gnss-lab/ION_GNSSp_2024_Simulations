%% Comments
% Project Title: EKF_KF_VAR_TS_Analysis

% Author: Rodrigo de Lima Florindo

% Version: 1.0

% Date: 22/06/2024

% Description: This code may be used to simulate the behaviour of
% multi-frequency KF-VAR and EKF-VAR carrier phase tracking loops during
% ionospheric scintillation for a single seeds of AWGN and MFPSM.

% Its outputs are:
% 1 - The values of Total RMSE of the main states
% 2 - A plot with the residuals estimates of LOS dynamics,
% 3 - A plot that compares the estimates of scintillation phase of both
% KF-VAR and EKF-VAR with the true one.

%% Initial setup
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
sigma2Measure_valueL1 = 1*(1/(2*CN0L1*T_I))*(1+(1/(2*CN0L1*T_I)));
sigma2Measure_valueL2 = 1*(1/(2*CN0L2*T_I))*(1+(1/(2*CN0L2*T_I)));
sigma2Measure_valueL5 = 1*(1/(2*CN0L5*T_I))*(1+(1/(2*CN0L5*T_I)));

% LOS paramaters
fd = 1000; % Frequency Shift
fdr = 0.94; % Frequency Drift
phiL1_0 = 0; % Initial phase shift for L1
phiL2_0 = 0; % Initial phase shift for L2
phiL5_0 = 0; % Initial phase shift for L5

M = 3; % Wiener process order

% Adaptative Module Configuration
% 1 - Non Adaptive
% 2 - Adaptive with noise
AdaptSwitch = 2;

%Scintillation Scenarios
% 1 - MFPSM Weak
% 2 - MFPSM Moderate
% 3 - MFPSM Severe
ScintScenario = 3;
    if ScintScenario == 1
        load('S4_0p35_tau_2p0_200_L1_L2_L5_300s.mat');
    elseif ScintScenario == 2
        load('S4_0p57_tau_1p4_200_L1_L2_L5_300s.mat');
    elseif ScintScenario == 3
        load('S4_0p8_tau_0p8_200_L1_L2_L5_300s.mat');
    end
    
%% MC Runs

% This represents the value of the variance related to the satellite
% acceleration that is used to model the LOS process covariance matrices.
sigma2 = 2.6*10^-6;
% Select:
% *** TopologySelector=1 to simulate the topologies #1.1 and #2.1;
% *** TopologySelector=2 to simulate the topologies #1.2 and #2.2;
% *** TopologySelector=3 to simulate the topologies #1.3 and #2.3;
TopologySelector = 3;
% Change SeedAWGN and ScintSeed to observe the behaviour of the models to
% other seeds.
seedAWGN = 4;
scintSeed = 4; % 16

% KF-VAR - SPECIFIC MODEL CONFIGURATIONS AND SIMULATIONS
run('KF_VAR_Specific_Config.m');
% Simulating the received signal with LOS dynamics, AWGN and scintillation
[yk,phi_D_true,rho_I_true,phi_I_true] = ReceivedSignal(TopologySelector,time,Scin_psi,AmpL1,AmpL2,AmpL5,T_I,f_receiver_Hz,fd,fdr,phiL1_0,phiL2_0,phiL5_0,seedAWGN,CN0dBHzL1,CN0dBHzL2,CN0dBHzL5,fcL1,fcL2,fcL5);

% Simulating the KF-VAR Model
[Px_0_0,x_hat_0_0,Qk,Fk,C] = KF_VAR_General_Config(L,M,NumSeries,P,T_I,sigma2,fd,fdr,deltaArray,Q_arfit,A,c,TopologySelector);
[xkk1_KF, Pkk1_KF] = KF_VAR_PLL(Px_0_0, x_hat_0_0, yk, Fk, Qk, C, Rk, Hk, AdaptSwitch, T_I,CN0L1,CN0L2,CN0L5,TopologySelector);
[~,~,~,ED_KF,EI_KF,ET_KF] = Evaluation_KFVAR(xkk1_KF,phi_D_true,phi_I_true,TopologySelector,L,M);

%EKF-VAR - SPECIFIC MODEL CONFIGURATIONS AND SIMULATIONS
run('EKF_VAR_Specific_Config.m');

% Simulating the EKF-VAR Model
[Px_0_0,x_hat_0_0,Qk,Fk,C] = EKF_VAR_General_Config(L,M,NumSeries,P,T_I,sigma2,fd,fdr,deltaArray,Q_arfit,A,c,TopologySelector);
[xkk1_EKF, Pkk1_EKF] = EKF_VAR_PLL(Px_0_0,x_hat_0_0,yk,Rk,Fk,Qk,M,L,NumSeries,P,AmpL1,AmpL2,AmpL5,Aux,C,TopologySelector,AdaptSwitch,CN0L1,CN0L2,CN0L5,T_I);
[~,~,~,~,ED_EKF,EI_EKF,ET_EKF,ERho_EKF] = Evaluation_EKFVAR(xkk1_EKF,phi_D_true,phi_I_true,rho_I_true,TopologySelector,L,M);

% % EKF-VAR-SEPARATE - SPECIFIC MODEL CONFIGURATIONS AND SIMULATIONS
% run('EKF_VAR_Separate_Specific_Config.m');

% Simulating the EKF-VAR-SEPARATE Model
% [Px_0_0,x_hat_0_0,Qk,Fk,C] = EKF_VAR_Separate_General_Config(L,M,NumSeries,P_rho,P_phi,T_I,sigma2,fd,fdr,deltaArray,Q_arfit_rho,Q_arfit_phi,A_rho,A_phi,c_rho,c_phi,TopologySelector);
% [xkk1_EKF_Sep, Pkk1_EKF_Sep] = EKF_VAR_Separate_PLL(Px_0_0,x_hat_0_0,yk,Rk,Fk,Qk,M,L,NumSeries,P_rho,P_phi,AmpL1,AmpL2,AmpL5,Aux,C,TopologySelector,AdaptSwitch,CN0L1,CN0L2,CN0L5,T_I);
% [~,~,~,~,ED_EKF_Sep,EI_EKF_Sep,ET_EKF_Sep,ERho_EKF_Sep] = Evaluation_Separate_EKFVAR(xkk1_EKF_Sep,phi_D_true,phi_I_true,rho_I_true,TopologySelector,L,M,P_rho,P_phi);

%% Plots
% Select:
% ---Topology #1.1 and #2.1:
% *** Lp=1 to observe L1 states;
% ---Topology #1.2 and #2.2:
% *** Lp=1 to observe L1 states; Lp=2 to observe L5 states;
% ---Topology #1.3 and #2.3:
% *** Lp=1 to observe L1 states; Lp=2 to observe L2 states; and Lp=3 to observe L5 states
Lp = 1;

% Display the Total RMSE of each parameter of our interest
disp(['LOS phase RMSE EKF: ', num2str(sqrt(mean(ED_EKF(:,Lp).^2,1)))]);
disp(['LOS phase RMSE KF: ', num2str(sqrt(mean(ED_KF(:,Lp).^2,1)))]);
disp(['Scint phase RMSE EKF: ', num2str(sqrt(mean(EI_EKF(:,Lp).^2,1)))]);
disp(['Scint phase RMSE KF: ', num2str(sqrt(mean(EI_KF(:,Lp).^2,1)))]);
disp(['Total phase RMSE EKF: ', num2str(sqrt(mean(wrapToPi(ET_EKF(:,Lp)).^2,1)))]);
disp(['Total phase RMSE KF: ', num2str(sqrt(mean(wrapToPi(ET_KF(:,Lp)).^2,1)))]);
disp(['Scintillation Amplitude RMSE EKF: ', num2str(sqrt(mean(ERho_EKF(:,Lp).^2,1)))]);
% First plot with 4 subplots
figure('Position',[50,50,800,650]);
t1 = tiledlayout(2,2, 'TileSpacing','compact');
nexttile;
plot(time,[ED_EKF(:,Lp),ED_KF(:,Lp)], 'LineWidth',2);
ylabel('$$\hat{\phi}_{D,L_1} - \phi_{D,L_1}$$ [rad]', 'Interpreter', 'latex');
xlabel('Time [s]');
set(gca, 'FontSize', 12);
nexttile;
plot(time,[EI_EKF(:,Lp),EI_KF(:,Lp)], 'LineWidth',2);
ylabel('$$\hat{\phi}_{\mathcal{I},L_1} - \phi_{\mathcal{I},L_1}$$ [rad]', 'Interpreter', 'latex');
xlabel('Time [s]');
set(gca, 'FontSize', 12);
nexttile;
plot(time,[wrapToPi(ET_EKF(:,Lp)),wrapToPi(ET_KF(:,Lp))], 'LineWidth',2);
ylabel('$$\hat{\phi}_{TL_1L1} - \phi_{T,L_1}$$ [rad]', 'Interpreter', 'latex');
xlabel('Time [s]');
legend({'EKF-VAR','KF-VAR'}, 'Location', 'best');
set(gca, 'FontSize', 12);
nexttile;
plot(time,[ERho_EKF(:,Lp)], 'LineWidth',2);
ylabel('$$\hat{\rho}_{L_1} - \rho_{L_1}$$', 'Interpreter', 'latex');
xlabel('Time [s]');
set(gca, 'FontSize', 12);

title(t1,'Performance of #1.3 and #2.3 - L1 - Severe scint (seed 1) - EKF-VAR and KF-VAR');

% Second Plot
figure;
plot(time,[xkk1_EKF(:,(L+M-1)+3+Lp),xkk1_KF(:,(L+M-1)+Lp),phi_I_true(:,Lp)],'LineWidth',2);
xlabel('Time [s]');
ylabel('$$\hat{\phi}_{I,L_1}$$ [rad]', 'Interpreter', 'latex');
legend({'$$\hat{\phi}_{I,L_1,EKF-VAR}$$','$$\hat{\phi}_{I,L_1,KF-VAR}$$','$$\phi_{I,L_1,True}$$'}, 'Location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
title('$$\hat{\phi}_{I,L_1}$$ - EKF-VAR and KF-VAR', 'Interpreter', 'latex');
    
% Third Plot
figure;
plot(time,[xkk1_EKF(:,(L+M-1)+Lp),rho_I_true(:,Lp)],'LineWidth',2);
xlabel('Time [s]');
ylabel('$$\hat{\phi}_{I,L_1}$$ [rad]', 'Interpreter', 'latex');
legend({'$$\hat{\rho}_{I,L_1,EKF-VAR}$$','$$\rho_{I,L_1,True}$$'}, 'Location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
title('$$\hat{\phi}_{I,L_1}$$ - EKF-VAR and KF-VAR', 'Interpreter', 'latex');

% % Display the Total RMSE of each parameter of our interest
% disp(['LOS phase RMSE EKF: ', num2str(sqrt(mean(ED_EKF(:,Lp).^2,1)))]);
% disp(['LOS phase RMSE EKF-SEP: ', num2str(sqrt(mean(ED_EKF_Sep(:,Lp).^2,1)))]);
% disp(['LOS phase RMSE KF: ', num2str(sqrt(mean(ED_KF(:,Lp).^2,1)))]);
% disp(['Scint phase RMSE EKF: ', num2str(sqrt(mean(EI_EKF(:,Lp).^2,1)))]);
% disp(['Scint phase RMSE EKF-SEP: ', num2str(sqrt(mean(EI_EKF_Sep(:,Lp).^2,1)))]);
% disp(['Scint phase RMSE KF: ', num2str(sqrt(mean(EI_KF(:,Lp).^2,1)))]);
% disp(['Total phase RMSE EKF: ', num2str(sqrt(mean(wrapToPi(ET_EKF(:,Lp)).^2,1)))]);
% disp(['Total phase RMSE EKF-SEP: ', num2str(sqrt(mean(wrapToPi(ET_EKF_Sep(:,Lp)).^2,1)))]);
% disp(['Total phase RMSE KF: ', num2str(sqrt(mean(wrapToPi(ET_KF(:,Lp)).^2,1)))]);
% disp(['Scintillation Amplitude RMSE EKF: ', num2str(sqrt(mean(ERho_EKF(:,Lp).^2,1)))]);
% disp(['Scintillation Amplitude RMSE EKF-SEP: ', num2str(sqrt(mean(ERho_EKF_Sep(:,Lp).^2,1)))]);
% % First plot with 4 subplots
% figure('Position',[50,50,800,650]);
% t1 = tiledlayout(2,2, 'TileSpacing','compact');
% nexttile;
% plot(time,[ED_EKF(:,Lp),ED_EKF_Sep(:,Lp),ED_KF(:,Lp)], 'LineWidth',2);
% ylabel('$$\hat{\phi}_{D,L_1} - \phi_{D,L_1}$$ [rad]', 'Interpreter', 'latex');
% xlabel('Time [s]');
% set(gca, 'FontSize', 12);
% nexttile;
% plot(time,[EI_EKF(:,Lp),EI_EKF_Sep(:,Lp),EI_KF(:,Lp)], 'LineWidth',2);
% ylabel('$$\hat{\phi}_{\mathcal{I},L_1} - \phi_{\mathcal{I},L_1}$$ [rad]', 'Interpreter', 'latex');
% xlabel('Time [s]');
% set(gca, 'FontSize', 12);
% nexttile;
% plot(time,[wrapToPi(ET_EKF(:,Lp)),wrapToPi(ET_EKF_Sep(:,Lp)),wrapToPi(ET_KF(:,Lp))], 'LineWidth',2);
% ylabel('$$\hat{\phi}_{TL_1L1} - \phi_{T,L_1}$$ [rad]', 'Interpreter', 'latex');
% xlabel('Time [s]');
% legend({'EKF-VAR','EKF-VAR-Sep','KF-VAR'}, 'Location', 'best');
% set(gca, 'FontSize', 12);
% nexttile;
% plot(time,[ERho_EKF(:,Lp),ERho_EKF_Sep(:,Lp)], 'LineWidth',2);
% ylabel('$$\hat{\rho}_{L_1} - \rho_{L_1}$$', 'Interpreter', 'latex');
% xlabel('Time [s]');
% set(gca, 'FontSize', 12);
% 
% title(t1,'Performance of #1.3 and #2.3 - L1 - Severe scint (seed 1) - EKF-VAR and KF-VAR');
% 
% % Second Plot
% figure;
% plot(time,[xkk1_EKF(:,(L+M-1)+3+Lp),xkk1_EKF_Sep(:,(L+M-1)+L*P_rho+Lp),xkk1_KF(:,(L+M-1)+Lp),phi_I_true(:,Lp)],'LineWidth',2);
% xlabel('Time [s]');
% ylabel('$$\hat{\phi}_{I,L_1}$$ [rad]', 'Interpreter', 'latex');
% legend({'$$\hat{\phi}_{I,L_1,EKF-VAR}$$','$$\hat{\phi}_{I,L_1,EKF-VAR-SEP}$$','$$\hat{\phi}_{I,L_1,KF-VAR}$$','$$\phi_{I,L_1,True}$$'}, 'Location', 'best', 'Interpreter', 'latex');
% set(gca, 'FontSize', 12);
% title('$$\hat{\phi}_{I,L_1}$$ - EKF-VAR and KF-VAR', 'Interpreter', 'latex');
%     
% % Third Plot
% figure;
% plot(time,[xkk1_EKF(:,(L+M-1)+Lp),xkk1_EKF_Sep(:,(L+M-1)+Lp),rho_I_true(:,Lp)],'LineWidth',2);
% xlabel('Time [s]');
% ylabel('$$\hat{\phi}_{I,L_1}$$ [rad]', 'Interpreter', 'latex');
% legend({'$$\hat{\rho}_{I,L_1,EKF-VAR}$$','$$\hat{\rho}_{I,L_1,EKF-VAR-SEP}$$','$$\rho_{I,L_1,True}$$'}, 'Location', 'best', 'Interpreter', 'latex');
% set(gca, 'FontSize', 12);
% title('$$\hat{\phi}_{I,L_1}$$ - EKF-VAR and KF-VAR', 'Interpreter', 'latex');