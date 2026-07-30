%% Script to run multilayer_random_v2.m multiple times
% v2: reward is computed at the end of heating (not cooling), and
% coolTime is randomly perturbed per trajectory and recorded.

clear; clc;

% Heat_MultipleScans.m (and other shared helpers) live in the parent
% LPBF-Simulation folder; add it to the path so it can still be found.
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

% ---------------------------
% Number of trajectories to generate
nTrajectories = 5000;

% Root folder to save trajectories
rootFolder = 'RL_Dataset_v2';
if ~exist(rootFolder,'dir')
    mkdir(rootFolder);
end

% ---------------------------
% Loop over trajectories
for trajID = 1:nTrajectories
    fprintf('\nRunning trajectory %d/%d\n', trajID, nTrajectories);


    tic;
    % Call the external function
    % multilayer_random_v2.m accepts 'trajID' and 'rootFolder'
    % and handles saving its own files
    multilayer_random_v2(trajID, rootFolder);

    elapsedTime = toc;  % end timer

    trajTimes(trajID) = elapsedTime;
    fprintf('Trajectory %d finished in %.2f seconds.\n', trajID, elapsedTime);
end

fprintf('\nAll %d trajectories completed.\n', nTrajectories);
