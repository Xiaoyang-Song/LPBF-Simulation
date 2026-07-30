function multilayer_random_v2(trajID, root_dir)
% MULTILAYER_RANDOM_V2  Same as multilayer_random.m, except:
%   - The reward (meanDeviation) is computed at the end of heating
%     (see simulateHeatingCooling_v2.m), not after cooling.
%   - coolTime is randomly sampled once per trajectory (from
%     coolTime_range below) instead of being a fixed constant, and the
%     sampled value is recorded per-layer and in the trajectory summary.

    % Heat_MultipleScans.m (and other shared helpers) live in the parent
    % LPBF-Simulation folder; add it to the path so it can still be found.
    addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

    paramsStruct.k = 0.0067;
    paramsStruct.rho = 0.00433;
    paramsStruct.specificHeat = 0.526;
    paramsStruct.ic = 300; % Initial condition: set to 300K as ambient temperature
    % paramsStruct.thick = 2.0;
    paramsStruct.thick = 4.0; % Thickness of the layer
    paramsStruct.width = 12;
    paramsStruct.height = 3;
    paramsStruct.hmax = 0.4;
    % Customize the shape here
    paramsStruct.squareSideFraction = 0.4;

    paramsStruct.scan_pattern = linspace(0,paramsStruct.width,48);
    paramsStruct.style = "simultaneous";
    paramsStruct.params.SS = 600;
    paramsStruct.params.LP = 100;
    paramsStruct.params.eeta = 0.3;
    paramsStruct.params.r_b = 0.06;
    paramsStruct.params.H = 0.1;
    % paramsStruct.heatTime = 0.05;
    paramsStruct.heatTime = 0.025;
    paramsStruct.nTimeStepsHeat = 50;
    paramsStruct.nTimeStepsCool = 50;
    paramsStruct.doPlot = false;
    paramsStruct.tempRange = [2000, 2800];  % target temperature range

    % ---------------------------
    % Randomly draw coolTime once for the whole trajectory (discrete, step 0.01)
    coolTime_values = 0.05:0.01:0.15;
    coolTime_traj = coolTime_values(randi(numel(coolTime_values)));
    paramsStruct.coolTime = coolTime_traj;

    % ---------------------------
    % Root folder for all trajectories
    rootFolder = root_dir;
    if ~exist(rootFolder,'dir')
        mkdir(rootFolder);
    end

    % ---------------------------
    % Define trajectory ID (increment automatically)
    existing = dir(fullfile(rootFolder,'trajectory_*'));
    trajectoryID = length(existing) + 1;
    trajectoryFolder = fullfile(rootFolder,sprintf('trajectory_%03d',trajectoryID));
    mkdir(trajectoryFolder);
    fprintf('Saving this trajectory in folder: %s\n', trajectoryFolder);
    fprintf('Cooling time for this trajectory: %.4f s\n', coolTime_traj);

    % ---------------------------
    % Start diary to save console output
    diaryFile = fullfile(trajectoryFolder,'console_output.txt');
    diary(diaryFile);
    diary on

    % ---------------------------
    % Define action ranges
    SS_range = [400, 1500];    % SS min/max
    LP_values = 100:10:400;  % discrete LP choices from 100 to 400 with step of 10

    nSteps = 12;  % number of steps

    % Define layer evolution
    initialFraction = 0.4;
    finalFraction = 0.5;
    fractions = linspace(initialFraction, finalFraction, nSteps);

    % ---------------------------
    % Initialize storage
    actions = zeros(nSteps,2);           % [SS, LP]
    meanDeviationList = zeros(1,nSteps); % Reward list (now end-of-heating)
    maxTempList = zeros(1,nSteps); % Not needed actually
    uFinalList = cell(1,nSteps);
    uHeatFinalList = cell(1,nSteps); % Temperature field at end of heating (reward input)
    tAllList = cell(1,nSteps);
    uAllList = cell(1,nSteps);
    resultsAllList = cell(1,nSteps);
    coolTimeList = zeros(1,nSteps); % Cooling time used for each layer

    % ---------------------------
    % Run simulations step by step
    for i = 1:nSteps
        % Random agent generates actions
        SS_action = SS_range(1) + rand()*(SS_range(2)-SS_range(1));
        LP_action = LP_values(randi(numel(LP_values)));

        actions(i,:) = [SS_action, LP_action];
        paramsStruct.params.SS = SS_action;
        paramsStruct.params.LP = LP_action;

        % Shape evolution
        paramsStruct.squareSideFraction = fractions(i);

        coolTimeList(i) = coolTime_traj; % same coolTime for every layer in this trajectory

        fprintf('Layer %d/%d -> SS: %.1f, LP: %.1f, coolTime: %.4f\n', i, nSteps, SS_action, LP_action, coolTime_traj);

        % Pass previous final temperature as initial condition
        if i == 1
            % First step: use ambient as indicated above
            [uFinal, tAll, uAll, resultAll, model, meanDeviation, uHeatFinal] = simulateHeatingCooling_v2(paramsStruct);
        else
            % Subsequent steps: use previous uFinal as IC
            paramsStruct.ic = resultAll(2); % Update initial condition for next step
            [uFinal, tAll, uAll, resultAll, model, meanDeviation, uHeatFinal] = simulateHeatingCooling_v2(paramsStruct);
        end

        fprintf('Layer %d/%d -> Mean Deviation (reward, end-of-heating): %.4f\n', i, nSteps, meanDeviation);

        % Store results
        meanDeviationList(i) = meanDeviation;
        maxTempList(i) = max(uFinal);
        uFinalList{i} = uFinal;
        uHeatFinalList{i} = uHeatFinal;
        tAllList{i} = tAll;
        uAllList{i} = uAll;
        resultsAllList{i} = resultAll;

        % Save per-step data
        coolTime_step = coolTime_traj;
        save(fullfile(trajectoryFolder,sprintf('layer_%d_data.mat',i)), ...
            'uFinal','uHeatFinal','tAll','uAll','resultAll','SS_action','LP_action','meanDeviation','coolTime_step');

        % Save final temperature figure: end-of-heating (left) vs end-of-cooling (right)
        % Figure is twice the default width (and default height) so each
        % subplot keeps the same size/scaling as the original single-panel plot.
        defaultFigPos = get(0,'defaultfigureposition');
        fig = figure('Visible','off','Position',[defaultFigPos(1) defaultFigPos(2) 2*defaultFigPos(3) defaultFigPos(4)]);

        subplot(1,2,1);
        pdeplot(model,'XYData',uHeatFinal,'Mesh','on','ColorMap','jet');
        colorbar; caxis([300 5000]);
        title('End of Heating');

        subplot(1,2,2);
        pdeplot(model,'XYData',uFinal,'Mesh','on','ColorMap','jet');
        colorbar; caxis([300 5000]);
        title(sprintf('End of Cooling (coolTime = %.4f s)', coolTime_traj));

        sgtitle(sprintf('Step %d: Heating vs Cooling', i));
        saveas(fig, fullfile(trajectoryFolder,sprintf('layer_%d_finalTemp.png',i)));
        close(fig);

    end

    % ---------------------------
    % Save trajectory summary
    save(fullfile(trajectoryFolder,'trajectory_summary.mat'), ...
        'actions','meanDeviationList','maxTempList','uFinalList','uHeatFinalList','tAllList','uAllList','resultsAllList', ...
        'coolTimeList','coolTime_traj');

    % ---------------------------
    % Summary
    fprintf('\nRandom Agent Simulation Complete.\n');
    for i = 1:nSteps
        uShape = size(uFinalList{i});  % get the shape of the temperature field
        fprintf('Step %d -> SS=%.1f, LP=%.1f | Mean Deviation (end-of-heating)=%.2f, Max Temp=%.2f | uFinal size: [%d x %d]\n', ...
            i, actions(i,1), actions(i,2), meanDeviationList(i), maxTempList(i), uShape(1), uShape(2));
    end
    fprintf('Cooling time used for this trajectory: %.4f s\n', coolTime_traj);

    % ---------------------------
    % End console log
    diary off

end
