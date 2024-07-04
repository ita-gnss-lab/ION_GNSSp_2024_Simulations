% Project Title: TimeSeriesExamplePlots

% Author: Rodrigo de Lima Florindo

% Version: 1.0

% Date: 22/06/2024

% Description: This code plots the time series of amplitude and phases of a
% realization of the MFPSM

%% Initial setup
clear all;
clc;

% Select the Ionospheric scintillation sceanario that you want to analyze
% 1 - Mild
% 2 - Moderate
% 3 - Severe
ScintScenario = 3;
if ScintScenario == 1
    load('S4_0p35_tau_2p0_200_L1_L2_L5_300s.mat');
elseif ScintScenario == 2
    load('S4_0p57_tau_1p4_200_L1_L2_L5_300s.mat');
elseif ScintScenario == 3
    load('S4_0p8_tau_0p8_200_L1_L2_L5_300s.mat');
end

%% Plots
% Select the MFPSM seed that you want to analyze
ScintSeed = 16;

plots = Y_obs_full(:,:,ScintSeed);
scen = 3;

t = 0:0.01:300;

figure;

subplot(2,1,1);
hold on;
ampL5 = plot(t,10*log10(plots(:,3).^2));
ampL2 = plot(t,10*log10(plots(:,2).^2));
ampL1 = plot(t,10*log10(plots(:,1).^2));
set(gca, 'FontName', 'Times', 'FontSize', 14);
set(ampL1, 'Color', 'r', 'LineWidth', 2);
set(ampL2, 'Color', 'g', 'LineWidth', 2);
set(ampL5, 'Color', 'b', 'LineWidth', 2);
if ScintScenario == 1
    ylabel('10log$$_{10}(\rho_l^2[k])$$ [dB]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 14);
end
grid on;

hold off;


subplot(2,1,2);
hold on
phiL5 = plot(t,plots(:,6));
phiL2 = plot(t,plots(:,5));
phiL1 = plot(t,plots(:,4));
hold off

subplot(2,1,2);
set(gca, 'FontName', 'Times', 'FontSize', 14);
set(phiL1, 'Color', 'r', 'LineWidth', 2);
set(phiL2, 'Color', 'g', 'LineWidth', 2);
set(phiL5, 'Color', 'b', 'LineWidth', 2);
if ScintScenario == 1
    ylabel('$$\phi_{\mathcal{I},l}[k]$$ [rad]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 14);
end
if ScintScenario == 2
    xlabel('Time [s]')
    legend({'L5','L2','L1'}, 'Location','best')
end

grid on;

% Define figure size in inches
width = 5;   % Width in inches
height = 6;  % Height in inches
% 
% Set the figure properties
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, width, height]);

% Set the paper size to match the figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width, height]);
set(gcf, 'PaperPosition', [0, 0, width, height]);
set(gcf, 'PaperPositionMode', 'manual');

axesHandles = findall(gcf, 'Type', 'axes');

for i = 1:length(axesHandles)
    set(axesHandles(i), 'FontName', 'Times', 'FontSize', 13);
    textHandles = findall(axesHandles(i), 'Type', 'text');
    set(textHandles, 'FontName', 'Times', 'FontSize', 13);
end

if ScintScenario == 1
    print('Example_Mild', '-dpdf', '-r300');
elseif ScintScenario == 2
    print('Example_Moderate', '-dpdf', '-r300');
elseif ScintScenario == 3
    print('Example_Severe', '-dpdf', '-r300');
end