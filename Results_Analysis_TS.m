%% Comments
% Project Title: Results_Analysis_TS

% Author: Rodrigo de Lima Florindo

% Version: 1.0

% Date: 22/06/2024

% Description: This code may be used to observe the performances of
% multi-frequency EKF-VAR and KF-VAR for mild, moderate and severe
% ionospheric scintillation scenarios through the data sets generated using
% the code with name "EKF_KF_VAR_MC_alt", using ARMSE (Averaged RMSE) as
% explained in the paper.

% Its outputs are:
% A plot with 12 subplots exploiting the ARMSE of each state of interest
%% Initial Setup

clear all;
clc;
addpath([cd,'\Results']);
% Select the scintillation scenario that you want to analyze
% 1 - MFPSM Weak
% 2 - MFPSM Moderate
% 3 - MFPSM Severe
ScintScenario = 1;
if ScintScenario == 1
    load('Results_Mild_RMSE_TS_v2_1.mat');
elseif ScintScenario == 2
    load('Results_Moderate_RMSE_TS_v2_1.mat');
elseif ScintScenario == 3
    load('Results_Severe_RMSE_TS_v2_1.mat');
end

%% Generating the Subplots
ResultsAux = cell(1,7);
for i = 1:7
    ResultsAux{1,i} = cell(1,3);
end

% Please, set the value of MC runs corresponding to the one used to
% generate the data sets.
MCruns = 1;
for i = 1:7
    for m = 1:3
        sumTot = zeros(30001,1);
        for k = 1:MCruns
            sumTot = sumTot + Results{1,i}{k,m};
        end
        ResultsAux{1,i}{1,m} = sumTot/MCruns;
    end
end

time = 0:0.01:300;

%% Plots

% Phi-D
figure;
t1 = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(t1);
hold on;
h1 = semilogy(time,ResultsAux{1,1}{1,1}(:,1),'LineWidth',2);
h2 = semilogy(time,ResultsAux{1,1}{1,2}(:,1),'LineWidth',2);
h3 = semilogy(time,ResultsAux{1,1}{1,3}(:,1),'LineWidth',2);
h4 = semilogy(time,ResultsAux{1,4}{1,1}(:,1),'LineWidth',2);
h5 = semilogy(time,ResultsAux{1,4}{1,2}(:,1),'LineWidth',2);
h6 = semilogy(time,ResultsAux{1,4}{1,3}(:,1),'LineWidth',2);
hold off;
set(h1, 'Color', '#00F5CC');
set(h2, 'Color', '#2DB59E');
set(h3, 'Color', '#499185');
set(h4, 'Color', '#F54500');
set(h5, 'Color', '#C05930');
set(h6, 'Color', '#9F6751');
title('L1');
ylabel('ARMSE$$(\hat{\phi}_{D,l})$$', 'interpreter', 'latex');
% legend({'#1.1','#1.2','#1.3','#2.1','#2.2','#2.3'},'Location','best')
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
set(gca, 'fontname', 'Times New Roman');


nexttile(t1);
hold on;
h3 = semilogy(time,ResultsAux{1,1}{1,3}(:,2),'LineWidth',2);
h6 = semilogy(time,ResultsAux{1,4}{1,3}(:,2),'LineWidth',2);
hold off;
set(h3, 'Color', '#499185');
set(h6, 'Color', '#9F6751');
title('L2');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
set(gca, 'fontname', 'Times New Roman');


nexttile(t1);
hold on;
h2 = semilogy(time,ResultsAux{1,1}{1,2}(:,2),'LineWidth',2);
h3 = semilogy(time,ResultsAux{1,1}{1,3}(:,3),'LineWidth',2);
h5 = semilogy(time,ResultsAux{1,4}{1,2}(:,2),'LineWidth',2);
h6 = semilogy(time,ResultsAux{1,4}{1,3}(:,3),'LineWidth',2);
hold off;
set(h2, 'Color', '#2DB59E');
set(h3, 'Color', '#499185');
set(h5, 'Color', '#C05930');
set(h6, 'Color', '#9F6751');
title('L5');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
set(gca, 'fontname', 'Times New Roman');


%Phi-I
nexttile(t1);
hold on;
h1 = semilogy(time,ResultsAux{1,2}{1,1}(:,1),'LineWidth',2);
h2 = semilogy(time,ResultsAux{1,2}{1,2}(:,1),'LineWidth',2);
h3 = semilogy(time,ResultsAux{1,2}{1,3}(:,1),'LineWidth',2);
h4 = semilogy(time,ResultsAux{1,5}{1,1}(:,1),'LineWidth',2);
h5 = semilogy(time,ResultsAux{1,5}{1,2}(:,1),'LineWidth',2);
h6 = semilogy(time,ResultsAux{1,5}{1,3}(:,1),'LineWidth',2);
hold off;
set(h1, 'Color', '#00F5CC');
set(h2, 'Color', '#2DB59E');
set(h3, 'Color', '#499185');
set(h4, 'Color', '#F54500');
set(h5, 'Color', '#C05930');
set(h6, 'Color', '#9F6751');
ylabel('ARMSE$$(\hat{\phi}_{\mathcal{I},l})$$', 'interpreter', 'latex');
% legend({'#1.1','#1.2','#1.3','#2.1','#2.2','#2.3'},'Location','best')
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
set(gca, 'fontname', 'Times New Roman');


nexttile(t1);
hold on;
h3 = semilogy(time,ResultsAux{1,2}{1,3}(:,2),'LineWidth',2);
h6 = semilogy(time,ResultsAux{1,5}{1,3}(:,2),'LineWidth',2);
hold off;
set(h3, 'Color', '#499185');
set(h6, 'Color', '#9F6751');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
set(gca, 'fontname', 'Times New Roman');


nexttile(t1);
hold on;
h2 = semilogy(time,ResultsAux{1,2}{1,2}(:,2),'LineWidth',2);
h3 = semilogy(time,ResultsAux{1,2}{1,3}(:,3),'LineWidth',2);
h5 = semilogy(time,ResultsAux{1,5}{1,2}(:,2),'LineWidth',2);
h6 = semilogy(time,ResultsAux{1,5}{1,3}(:,3),'LineWidth',2);
hold off;
set(h2, 'Color', '#2DB59E');
set(h3, 'Color', '#499185');
set(h5, 'Color', '#C05930');
set(h6, 'Color', '#9F6751');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
set(gca, 'fontname', 'Times New Roman');


% Legend formatting
hold on;
h1 = semilogy(nan, nan, 'Color','#00F5CC', 'DisplayName', '#1.1','LineWidth',2);
h2 = semilogy(nan, nan, 'Color','#2DB59E', 'DisplayName', '#1.2','LineWidth',2);
h3 = semilogy(nan, nan, 'Color','#499185', 'DisplayName', '#1.3','LineWidth',2);
h4 = semilogy(nan, nan, 'Color','#F54500', 'DisplayName', '#2.1','LineWidth',2);
h5 = semilogy(nan, nan, 'Color','#C05930', 'DisplayName', '#2.2','LineWidth',2);
h6 = semilogy(nan, nan, 'Color','#9F6751', 'DisplayName', '#2.3','LineWidth',2);
hold off;
set(h1, 'Color', '#00F5CC');
set(h2, 'Color', '#2DB59E');
set(h3, 'Color', '#499185');
set(h4, 'Color', '#F54500');
set(h5, 'Color', '#C05930');
set(h6, 'Color', '#9F6751');
% Create a new axis for the legend
legend([h1, h2, h3, h4, h5, h6], 'Location', 'eastoutside');
set(gca, 'fontname', 'Times New Roman');


%Phi-T
nexttile(t1);
hold on;
h1 = semilogy(time,ResultsAux{1,3}{1,1}(:,1),'LineWidth',2);
h2 = semilogy(time,ResultsAux{1,3}{1,2}(:,1),'LineWidth',2);
h3 = semilogy(time,ResultsAux{1,3}{1,3}(:,1),'LineWidth',2);
h4 = semilogy(time,ResultsAux{1,6}{1,1}(:,1),'LineWidth',2);
h5 = semilogy(time,ResultsAux{1,6}{1,2}(:,1),'LineWidth',2);
h6 = semilogy(time,ResultsAux{1,6}{1,3}(:,1),'LineWidth',2);
hold off;
set(h1, 'Color', '#00F5CC');
set(h2, 'Color', '#2DB59E');
set(h3, 'Color', '#499185');
set(h4, 'Color', '#F54500');
set(h5, 'Color', '#C05930');
set(h6, 'Color', '#9F6751');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
ylabel('ARMSE$$(\hat{\phi}_{T,l})$$', 'interpreter', 'latex');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
set(gca, 'fontname', 'Times New Roman');
set(gca, 'fontname', 'Times New Roman');


nexttile(t1);
hold on;
h3 = semilogy(time,ResultsAux{1,3}{1,3}(:,2),'LineWidth',2);
h6 = semilogy(time,ResultsAux{1,6}{1,3}(:,2),'LineWidth',2);
hold off;
set(h3, 'Color', '#499185');
set(h6, 'Color', '#9F6751');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
set(gca, 'fontname', 'Times New Roman');


nexttile(t1);
hold on;
h2 = semilogy(time,ResultsAux{1,3}{1,2}(:,2),'LineWidth',2);
h3 = semilogy(time,ResultsAux{1,3}{1,3}(:,3),'LineWidth',2);
h5 = semilogy(time,ResultsAux{1,6}{1,2}(:,2),'LineWidth',2);
h6 = semilogy(time,ResultsAux{1,6}{1,3}(:,3),'LineWidth',2);
hold off;
set(h2, 'Color', '#2DB59E');
set(h3, 'Color', '#499185');
set(h5, 'Color', '#C05930');
set(h6, 'Color', '#9F6751');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
set(gca, 'fontname', 'Times New Roman');


%Rho
nexttile(t1);
hold on;
h4 = semilogy(time,ResultsAux{1,7}{1,1}(:,1),'LineWidth',2);
h5 = semilogy(time,ResultsAux{1,7}{1,2}(:,1),'LineWidth',2);
h6 = semilogy(time,ResultsAux{1,7}{1,3}(:,1),'LineWidth',2);
hold off;
set(h4, 'Color', '#F54500');
set(h5, 'Color', '#C05930');
set(h6, 'Color', '#9F6751');
ylabel('ARMSE$$(\hat{\rho}_{l})$$', 'interpreter', 'latex');
% legend({'#1.1','#1.2','#1.3','#2.1','#2.2','#2.3'},'Location','best')
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
set(gca, 'fontname', 'Times New Roman');


nexttile(t1);
hold on;
h6 = semilogy(time,ResultsAux{1,7}{1,3}(:,2),'LineWidth',2);
hold off;
set(h6, 'Color', '#9F6751');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
xlabel('Time [s]','interpreter', 'latex');
set(gca, 'fontname', 'Times New Roman');


nexttile(t1);
hold on;
h5 = semilogy(time,ResultsAux{1,7}{1,2}(:,2),'LineWidth',2);
h6 = semilogy(time,ResultsAux{1,7}{1,3}(:,3),'LineWidth',2);
hold off;
set(h5, 'Color', '#C05930');
set(h6, 'Color', '#9F6751');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca, 'TickLength', [.03,.03]);
set(gca, 'FontSize', 13);
set(gca, 'fontname', 'Times New Roman');


% Define figure size in inches
width = 10;   % Width in inches
height = 7.5;  % Height in inches

% Set the figure properties
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, width, height]);

% Set the paper size to match the figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width, height]);
set(gcf, 'PaperPosition', [0, 0, width, height]);
set(gcf, 'PaperPositionMode', 'manual');

if ScintScenario == 1
    print([cd,'\Results\Results_Mild'], '-dpdf', '-r300');
elseif ScintScenario == 2
    print([cd,'\Results\Results_Moderate'], '-dpdf', '-r300');
elseif ScintScenario == 3
    print([cd,'\Results\Results_Severe'], '-dpdf', '-r300');
end

