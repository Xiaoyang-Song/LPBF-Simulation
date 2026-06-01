function multilayer_random_cooling(trajID, root_dir, coolTime)
% multilayer_random_cooling  Same as multilayer_random but coolTime is
% supplied by the caller rather than hardcoded, and root_dir is actually
% used as the output root folder.

    paramsStruct.k = 0.0067;
    paramsStruct.rho = 0.00433;
    paramsStruct.specificHeat = 0.526;
    paramsStruct.ic = 300;
    paramsStruct.thick = 2.0;
    paramsStruct.width = 12;
    paramsStruct.height = 3;
    paramsStruct.hmax = 0.4;
    paramsStruct.squareSideFraction = 0.4;

    paramsStruct.scan_pattern = linspace(0,paramsStruct.width,48);
    paramsStruct.style = "simultaneous";
    paramsStruct.params.SS = 600;
    paramsStruct.params.LP = 100;
    paramsStruct.params.eeta = 0.3;
    paramsStruct.params.r_b = 0.06;
    paramsStruct.params.H = 0.1;
    paramsStruct.heatTime = 0.05;
    paramsStruct.coolTime = coolTime;   % <-- parameterised
    paramsStruct.nTimeStepsHeat = 50;
    paramsStruct.nTimeStepsCool = 50;
    paramsStruct.doPlot = false;
    paramsStruct.tempRange = [2000, 2800];

    % ---------------------------
    rootFolder = root_dir;
    if ~exist(rootFolder,'dir')
        mkdir(rootFolder);
    end

    % ---------------------------
    existing = dir(fullfile(rootFolder,'trajectory_*'));
    trajectoryID = length(existing) + 1;
    trajectoryFolder = fullfile(rootFolder,sprintf('trajectory_%03d',trajectoryID));
    mkdir(trajectoryFolder);
    fprintf('Saving this trajectory in folder: %s\n', trajectoryFolder);

    % ---------------------------
    diaryFile = fullfile(trajectoryFolder,'console_output.txt');
    diary(diaryFile);
    diary on

    % ---------------------------
    SS_range = [400, 1500];
    LP_values = 150:10:400;

    nSteps = 12;

    initialFraction = 0.4;
    finalFraction = 0.5;
    fractions = linspace(initialFraction, finalFraction, nSteps);

    % ---------------------------
    actions = zeros(nSteps,2);
    meanDeviationList = zeros(1,nSteps);
    maxTempList = zeros(1,nSteps);
    uFinalList = cell(1,nSteps);
    tAllList = cell(1,nSteps);
    uAllList = cell(1,nSteps);
    resultsAllList = cell(1,nSteps);

    % ---------------------------
    for i = 1:nSteps
        SS_action = SS_range(1) + rand()*(SS_range(2)-SS_range(1));
        LP_action = LP_values(randi(numel(LP_values)));

        actions(i,:) = [SS_action, LP_action];
        paramsStruct.params.SS = SS_action;
        paramsStruct.params.LP = LP_action;
        paramsStruct.squareSideFraction = fractions(i);

        fprintf('Layer %d/%d -> SS: %.1f, LP: %.1f\n', i, nSteps, SS_action, LP_action);

        if i == 1
            [uFinal, tAll, uAll, resultAll, model, meanDeviation] = simulateHeatingCooling(paramsStruct);
        else
            paramsStruct.ic = resultAll(2);
            [uFinal, tAll, uAll, resultAll, model, meanDeviation] = simulateHeatingCooling(paramsStruct);
        end

        meanDeviationList(i) = meanDeviation;
        maxTempList(i) = max(uFinal);
        uFinalList{i} = uFinal;
        tAllList{i} = tAll;
        uAllList{i} = uAll;
        resultsAllList{i} = resultAll;

        save(fullfile(trajectoryFolder,sprintf('layer_%d_data.mat',i)), ...
            'uFinal','tAll','uAll','resultAll','SS_action','LP_action','meanDeviation');

        fig = figure('Visible','off');
        pdeplot(model,'XYData',uFinal,'Mesh','on','ColorMap','jet');
        colorbar; caxis([300 5000]);
        title(sprintf('Step %d: Cooling Final Temperature',i));
        saveas(fig, fullfile(trajectoryFolder,sprintf('layer_%d_finalTemp.png',i)));
        close(fig);
    end

    % ---------------------------
    save(fullfile(trajectoryFolder,'trajectory_summary.mat'), ...
        'actions','meanDeviationList','maxTempList','uFinalList','tAllList','uAllList','resultsAllList');

    % ---------------------------
    fprintf('\nRandom Agent Simulation Complete.\n');
    for i = 1:nSteps
        uShape = size(uFinalList{i});
        fprintf('Step %d -> SS=%.1f, LP=%.1f | Mean Deviation=%.2f, Max Temp=%.2f | uFinal size: [%d x %d]\n', ...
            i, actions(i,1), actions(i,2), meanDeviationList(i), maxTempList(i), uShape(1), uShape(2));
    end

    diary off

end
