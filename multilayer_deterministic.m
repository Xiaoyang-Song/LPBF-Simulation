clear; close all;

paramsStruct.k = 0.0067;          
paramsStruct.rho = 0.00433;       
paramsStruct.specificHeat = 0.526;
paramsStruct.ic = 300; % Initial condition: set to 300K as ambient temperature        
paramsStruct.thick = 2.0;         
paramsStruct.width = 12;          
paramsStruct.height = 3;           
paramsStruct.hmax = 0.3;          
paramsStruct.squareSideFraction = 0.4;  
paramsStruct.scan_pattern = linspace(0,paramsStruct.width,48);  
paramsStruct.style = "simultaneous";    
paramsStruct.params.SS = 600;
paramsStruct.params.LP = 100;
paramsStruct.params.eeta = 0.3;
paramsStruct.params.r_b = 0.06;
paramsStruct.params.H = 0.1;
paramsStruct.heatTime = 0.05;
paramsStruct.coolTime = 0.10;
paramsStruct.nTimeStepsHeat = 50;
paramsStruct.nTimeStepsCool = 50;
paramsStruct.doPlot = true;  
paramsStruct.tempRange = [400, 450];  % target temperature range

% Single run test

% [uFinal, tAll, uAll, resultAll, model, meanDeviation] = simulateHeatingCooling(paramsStruct);
% 
% fprintf('Maximum temperature at the end of cooling: %.2f K\n', max(uFinal));
% fprintf('Mean deviation from range [%d, %d] K: %.2f K\n', ...
%     paramsStruct.tempRange(1), paramsStruct.tempRange(2), meanDeviation);


% ---------------------------
% Define laser powers and speeds for each step
laserPowerList = [400, 500, 600, 700];   
scanSpeedList  = [100, 100, 100, 100];    
nSteps = length(laserPowerList);

% ---------------------------
% Initialize storage
meanDeviationList = zeros(1,nSteps); % Reward list
maxTempList = zeros(1,nSteps); % Not needed actually
uFinalList = cell(1,nSteps);
tAllList = cell(1,nSteps);
uAllList = cell(1,nSteps);
resultsAllList = cell(1,nSteps);

% ---------------------------
% Run simulations step by step
for i = 1:nSteps
    fprintf('Step %d/%d: SS = %d, LP = %d\n', i, nSteps, laserPowerList(i), scanSpeedList(i));
    
    % Update laser parameters
    paramsStruct.params.SS = laserPowerList(i);
    paramsStruct.params.LP = scanSpeedList(i);
    
    % Pass previous final temperature as initial condition
    if i == 1
        % First step: use ambient as indicated above
        [uFinal, tAll, uAll, resultAll, model, meanDeviation] = simulateHeatingCooling(paramsStruct);
    else
        % Subsequent steps: use previous uFinal as IC
        paramsStruct.ic = resultAll(2); % Update initial condition for next step
        [uFinal, tAll, uAll, resultAll, model, meanDeviation] = simulateHeatingCooling(paramsStruct);
    end
    
    % Store results
    meanDeviationList(i) = meanDeviation;
    maxTempList(i) = max(uFinal);
    uFinalList{i} = uFinal;
    tAllList{i} = tAll;
    uAllList{i} = uAll;
    resultsAllList{i} = resultAll;
end

% ---------------------------
% Summary
fprintf('\nSimulation complete.\n');
for i = 1:nSteps
    fprintf('Step %d: SS=%d, LP=%d -> Mean Deviation=%.2f K, Max Temp=%.2f K\n', ...
        i, laserPowerList(i), scanSpeedList(i), meanDeviationList(i), maxTempList(i));
end