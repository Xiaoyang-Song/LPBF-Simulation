clear; close all; clc;

% ---------------------------
% Material and simulation parameters
paramsStruct.k = 0.0067;          
paramsStruct.rho = 0.00433;       
paramsStruct.specificHeat = 0.526;
paramsStruct.ic = 300;            % Initial temperature (K)
paramsStruct.thick = 2.0;         
paramsStruct.width = 12;          
paramsStruct.height = 3;          
paramsStruct.hmax = 0.3;          
paramsStruct.squareSideFraction = 0.4;  
paramsStruct.scan_pattern = linspace(0, paramsStruct.width, 48);  
paramsStruct.style = "simultaneous";    

% Laser parameters (fixed power, variable speed)
paramsStruct.params.LP = 100;     % Laser power (W) fixed
paramsStruct.params.eeta = 0.3;   % Absorption efficiency
paramsStruct.params.r_b = 0.06;   % Beam radius (mm)
paramsStruct.params.H = 0.1;      % Hatch spacing (mm)

% Thermal simulation settings
paramsStruct.heatTime = 0.05;
paramsStruct.coolTime = 0.10;
paramsStruct.nTimeStepsHeat = 50;
paramsStruct.nTimeStepsCool = 50;
paramsStruct.doPlot = false;  
paramsStruct.tempRange = [2500, 3000];  % Target temperature range

% ---------------------------
% Scanning speed values to test (mm/s)
scanSpeedList = linspace(200, 1200, 20);
nSteps = length(scanSpeedList);

% ---------------------------
% Compute energy density for each scanning speed
% ED = LP / (v * h * t)
laserPower = paramsStruct.params.LP;
energyDensityList = laserPower ./ ...
    (scanSpeedList * paramsStruct.params.H * paramsStruct.thick);  % J/mm³

% ---------------------------
% Initialize storage
meanDeviationList = zeros(1, nSteps);
maxTempList = zeros(1, nSteps);

% ---------------------------
% Run one-layer simulations
for i = 1:nSteps
    v = scanSpeedList(i);
    paramsStruct.params.SS = v;
    
    fprintf('Running case %d/%d: Scan Speed = %.1f mm/s, Energy Density = %.3f J/mm³\n', ...
        i, nSteps, v, energyDensityList(i));
    
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
    fprintf('Speed = %.1f mm/s, Energy Density = %.3f J/mm³ -> Mean Deviation = %.2f K, Max Temp = %.2f K\n', ...
        scanSpeedList(i), energyDensityList(i), meanDeviationList(i), maxTempList(i));
end

% ---------------------------
% Plot Reward vs. Energy Density
figure;
plot(energyDensityList, -meanDeviationList, 'o-', 'LineWidth', 1.8, 'MarkerSize', 8);
xlabel('Energy Density (J/mm³)');
ylabel('Reward (-Mean Deviation from Target Temp)');
title('Reward vs. Energy Density (Varying Scanning Speed)');
grid on;
