%% Script to run multilayer_random.m multiple times

clear; clc;

% ---------------------------
% Number of trajectories to generate
nTrajectories = 2000;

% Root folder to save trajectories
rootFolder = 'RL_Dataset';
if ~exist(rootFolder,'dir')
    mkdir(rootFolder);
end

% ---------------------------
% Loop over trajectories
for trajID = 1:nTrajectories
    fprintf('\nRunning trajectory %d/%d\n', trajID, nTrajectories);


    tic; 
    % Call the external function
    % Assume multilayer_random.m is written to accept 'trajID' and 'rootFolder'
    % and it handles saving its own files
    multilayer_random(trajID, rootFolder);

    elapsedTime = toc;  % end timer
    
    trajTimes(trajID) = elapsedTime;
    fprintf('Trajectory %d finished in %.2f seconds.\n', trajID, elapsedTime);
end

fprintf('\nAll %d trajectories completed.\n', nTrajectories);