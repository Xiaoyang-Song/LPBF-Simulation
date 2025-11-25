clear; close all; clc;

% ---------------------------
% Material and simulation parameters
paramsStruct.k = 0.0067;          
paramsStruct.rho = 0.00433;       
paramsStruct.specificHeat = 0.526;
paramsStruct.ic = 300;            % Initial temperature (K)
paramsStruct.thick = 1.0;         
paramsStruct.width = 12;          
paramsStruct.height = 3;          
paramsStruct.hmax = 0.3;          
paramsStruct.squareSideFraction = 0.4;  
paramsStruct.scan_pattern = linspace(0, paramsStruct.width, 48);  
paramsStruct.style = "simultaneous";    

% Laser scanning parameters (fixed speed, variable power)
paramsStruct.params.SS = 600;     % Scanning speed (mm/s) fixed
paramsStruct.params.eeta = 0.3;   % Absorption efficiency
paramsStruct.params.r_b = 0.06;   % Beam radius (mm)
paramsStruct.params.H = 0.1;      % Hatch spacing (mm)

% Thermal simulation settings
paramsStruct.heatTime = 0.05;
paramsStruct.coolTime = 0.10;
paramsStruct.nTimeStepsHeat = 50;
paramsStruct.nTimeStepsCool = 50;
% paramsStruct.doPlot = true;  
paramsStruct.doPlot = false;  
paramsStruct.tempRange = [2500, 3000];  % Target temperature range

% ---------------------------
% Laser power values to test
laserPowerList = linspace(100, 600, 10);  % (W)
nSteps = length(laserPowerList);

% ---------------------------
% Fixed scanning speed
scanSpeed = paramsStruct.params.SS; % fixed speed

% Compute energy density for each laser power
% Common definition: ED = (LP) / (v * h * t)
% where LP = laser power (W), v = scanning speed (mm/s),
% h = hatch spacing (mm), t = layer thickness (mm)
energyDensityList = laserPowerList ./ ...
    (scanSpeed * paramsStruct.params.H * paramsStruct.thick);  % J/mm³

% ---------------------------
% Initialize storage
meanDeviationList = zeros(1, nSteps);
maxTempList = zeros(1, nSteps);

% ---------------------------
% Run one-layer simulations
for i = 1:nSteps
    LP = laserPowerList(i);
    paramsStruct.params.LP = LP;
    
    fprintf('Running case %d/%d: Laser Power = %.1f W, Energy Density = %.3f J/mm³\n', ...
        i, nSteps, LP, energyDensityList(i));
    
    % Single-layer simulation
    [uFinal, tAll, uAll, resultAll, model, meanDeviation] = simulateHeatingCooling(paramsStruct);
    
    % Store metrics
    meanDeviationList(i) = meanDeviation;
    maxTempList(i) = max(uFinal);
end

% ---------------------------
% Summary of results
fprintf('\nSimulation complete.\n');
for i = 1:nSteps
    fprintf('LP = %.1f W, Energy Density = %.3f J/mm³ -> Mean Deviation = %.2f K, Max Temp = %.2f K\n', ...
        laserPowerList(i), energyDensityList(i), meanDeviationList(i), maxTempList(i));
end

% ---------------------------
% Plot Reward vs. Energy Density
figure;
plot(energyDensityList, -meanDeviationList, 'o-', 'LineWidth', 1.8, 'MarkerSize', 8);
xlabel('Energy Density (J/mm³)');
ylabel('Reward (-Mean Deviation from Target Temp)');
title('Reward vs. Energy Density');
grid on;
