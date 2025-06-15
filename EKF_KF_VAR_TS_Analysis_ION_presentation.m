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
ScintScenario = 1;
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
% Change SeedAWGN and ScintSeed to observe the behaviour of the models to
% other seeds.
seedAWGN = 1;
scintSeed = 16; % 16
phi_D_topos_arcs = zeros(30001,6);
phi_I_topos_arcs = zeros(30001,6);
rho_topos_arcs = zeros(30001,3);
for topo = 1:3
    TopologySelector = topo;
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
    
    phi_D_topos_arcs(:,1 + 2*(topo-1)) = xkk1_KF(:,1);
    phi_D_topos_arcs(:,2 + 2*(topo-1)) = xkk1_EKF(:,1);
    phi_I_topos_arcs(:,1 + 2*(topo-1)) = xkk1_KF(:,(L+M));
    phi_I_topos_arcs(:,2 + 2*(topo-1)) = xkk1_EKF(:,(L+M)+topo);
    rho_topos_arcs(:,topo) = xkk1_EKF(:,(L+M));
end

%% Plots
% First plot with 4 subplots
figure('Position',[50,50,800,650]);
t1 = tiledlayout(3,1, 'TileSpacing','compact');
nexttile;
% title(t1,'Performance of all models for L1 state estimates - Severe - S_4 = 0.8; \tau_0 = 0.8', 'FontSize', 17);
h1 = plot(time,phi_D_topos_arcs-repmat(phi_D_true(:,1),1,6),'LineWidth',3);
set(h1(1),'Color','#00F5CC','LineStyle','-');
set(h1(2),'Color','#F54500','LineStyle','-.');
set(h1(3),'Color','#2DB59E','LineStyle','-');
set(h1(4),'Color','#C05930','LineStyle','-.');
set(h1(5),'Color','#499185','LineStyle','-');
set(h1(6),'Color','#9F6751','LineStyle','-.');
black = yline(0);
set(black,'Color','b','LineStyle','--','LineWidth',3);
ylabel('LOS dyn. error [rad]', 'interpreter', 'latex');
% legend({'KF single-frequency','EKF single-frequency','KF dual-frequency','EKF dual-frequency','KF triple-frequency','EKF triple-frequency'},'Location','best');
set(gca, 'FontSize', 14);
set(gca, 'fontname', 'Times New Roman');

nexttile;
h2 = plot(time,[phi_I_topos_arcs,phi_I_true(:,1)],'LineWidth',3);
set(h2(1),'Color','#00F5CC','LineStyle','-');
set(h2(2),'Color','#F54500','LineStyle','-.');
set(h2(3),'Color','#2DB59E','LineStyle','-');
set(h2(4),'Color','#C05930','LineStyle','-.');
set(h2(5),'Color','#499185','LineStyle','-');
set(h2(6),'Color','#9F6751','LineStyle','-.');
set(h2(7),'Color','b','LineStyle','--','LineWidth',3);
ylabel('Ion. Phase [rad]', 'interpreter', 'latex')
%legend({'KF single-frequency','EKF single-frequency','KF dual-frequency','EKF dual-frequency','KF triple-frequency','EKF triple-frequency','Zero'},'Location','best');
set(gca, 'FontSize', 14);
set(gca, 'fontname', 'Times New Roman');

nexttile;
h3 = plot(time,10*log10([rho_topos_arcs,rho_I_true(:,1)].^2),'LineWidth',3);
set(h3(1),'Color','#F54500','LineStyle','-.');
set(h3(2),'Color','#C05930','LineStyle','-.');
set(h3(3),'Color','#9F6751','LineStyle','-.');
set(h3(4),'Color','b','LineStyle','--','LineWidth',3);
ylabel('Magnitude [dB]', 'interpreter', 'latex');
xlabel('Time [sec]', 'interpreter', 'latex');
%legend({'EKF single-frequency','EKF dual-frequency','EKF triple-frequency','Zero'},'Location','best');
set(gca, 'FontSize', 14);
set(gca, 'fontname', 'Times New Roman');

% Define figure size in inches
width = 16*0.67;   % Width in inches
height = 9*0.67;  % Height in inches

% Set the figure properties
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, width, height]);

% Set the paper size to match the figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width, height]);
set(gcf, 'PaperPosition', [0, 0, width, height]);
set(gcf, 'PaperPositionMode', 'manual');

if ScintScenario == 1
    print([cd,'\Results_ION\Results_Mild'], '-dpng', '-r300');
elseif ScintScenario == 2
    print([cd,'\Results_ION\Results_Moderate'], '-dpng', '-r300');
elseif ScintScenario == 3
    print([cd,'\Results_ION\Results_Severe'], '-dpng', '-r300');
end
