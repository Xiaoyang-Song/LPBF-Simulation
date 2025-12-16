function multilayer_random(trajID, root_dir)

    paramsStruct.k = 0.0067;          
    paramsStruct.rho = 0.00433;       
    paramsStruct.specificHeat = 0.526;
    paramsStruct.ic = 300; % Initial condition: set to 300K as ambient temperature        
    paramsStruct.thick = 2.0;         
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
    paramsStruct.heatTime = 0.05;
    paramsStruct.coolTime = 0.10;
    paramsStruct.nTimeStepsHeat = 50;
    paramsStruct.nTimeStepsCool = 50;
    paramsStruct.doPlot = false;  
    paramsStruct.tempRange = [2500, 3000];  % target temperature range
    
    % Single run test
    
    % [uFinal, tAll, uAll, resultAll, model, meanDeviation] = simulateHeatingCooling(paramsStruct);
    % 
    % fprintf('Maximum temperature at the end of cooling: %.2f K\n', max(uFinal));
    % fprintf('Mean deviation from range [%d, %d] K: %.2f K\n', ...
    %     paramsStruct.tempRange(1), paramsStruct.tempRange(2), meanDeviation);
    
   
    % ---------------------------
    % Root folder for all trajectories
    rootFolder = 'RL_Dataset';
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
    
    % ---------------------------
    % Start diary to save console output
    diaryFile = fullfile(trajectoryFolder,'console_output.txt');
    diary(diaryFile);
    diary on
    
    % ---------------------------
    % Define action ranges
    SS_range = [400, 1500];    % SS min/max
    % LP_range = [100, 600];     % LP min/max
    % LP_values = [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];  % discrete LP choices
    LP_values = 100:20:600;

    nSteps = 8;  % number of steps

    % Define layer evolution
    initialFraction = 0.4;
    finalFraction = 0.5;
    fractions = linspace(initialFraction, finalFraction, nSteps);

    % ---------------------------
    % Initialize storage
    actions = zeros(nSteps,2);           % [SS, LP]
    meanDeviationList = zeros(1,nSteps); % Reward list
    maxTempList = zeros(1,nSteps); % Not needed actually
    uFinalList = cell(1,nSteps);
    tAllList = cell(1,nSteps);
    uAllList = cell(1,nSteps);
    resultsAllList = cell(1,nSteps);
    
    % ---------------------------
    % Run simulations step by step
    for i = 1:nSteps
        % Random agent generates actions
        SS_action = SS_range(1) + rand()*(SS_range(2)-SS_range(1));
        % LP_action = LP_range(1) + rand()*(LP_range(2)-LP_range(1));
        LP_action = LP_values(randi(numel(LP_values)));

        
        actions(i,:) = [SS_action, LP_action];
        paramsStruct.params.SS = SS_action;
        paramsStruct.params.LP = LP_action;

        % Shape evolution
        paramsStruct.squareSideFraction = fractions(i);
        
        fprintf('Layer %d/%d -> SS: %.1f, LP: %.1f\n', i, nSteps, SS_action, LP_action);
    
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
    
        % Save per-step data
        save(fullfile(trajectoryFolder,sprintf('layer_%d_data.mat',i)), ...
            'uFinal','tAll','uAll','resultAll','SS_action','LP_action','meanDeviation');
        % Save final temperature figure
        fig = figure('Visible','off');
        pdeplot(model,'XYData',uFinal,'Mesh','on','ColorMap','jet');
        colorbar; caxis([300 5000]);
        title(sprintf('Step %d: Cooling Final Temperature',i));
        saveas(fig, fullfile(trajectoryFolder,sprintf('layer_%d_finalTemp.png',i)));
        close(fig);
    
    end
    
    % ---------------------------
    % Save trajectory summary
    save(fullfile(trajectoryFolder,'trajectory_summary.mat'), ...
        'actions','meanDeviationList','maxTempList','uFinalList','tAllList','uAllList','resultsAllList');
    
    % ---------------------------
    % Summary
    fprintf('\nRandom Agent Simulation Complete.\n');
    for i = 1:nSteps
        uShape = size(uFinalList{i});  % get the shape of the temperature field
        fprintf('Step %d -> SS=%.1f, LP=%.1f | Mean Deviation=%.2f, Max Temp=%.2f | uFinal size: [%d x %d]\n', ...
            i, actions(i,1), actions(i,2), meanDeviationList(i), maxTempList(i), uShape(1), uShape(2));
    end
    
    % ---------------------------
    % End console log
    diary off

end