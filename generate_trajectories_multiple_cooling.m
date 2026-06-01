%% Generate 300 trajectories for each of 11 cooling times (0.05 to 0.15 step 0.01)
% Each cooling time gets its own folder: RL_Dataset_ct_005, ..., RL_Dataset_ct_015

clear; clc;

nSamplesPerCool = 300;
coolTimes = 0.05:0.01:0.15;   % [0.05, 0.06, ..., 0.15]

for ci = 1:numel(coolTimes)
    ct = coolTimes(ci);

    % Folder name: RL_Dataset_ct_005 for 0.05, ct_015 for 0.15, etc.
    rootFolder = sprintf('RL_Dataset_ct_%03d', round(ct * 1000));

    if ~exist(rootFolder,'dir')
        mkdir(rootFolder);
    end

    fprintf('\n=== coolTime = %.2f  ->  folder: %s ===\n', ct, rootFolder);

    trajTimes = zeros(1, nSamplesPerCool);

    for trajID = 1:nSamplesPerCool
        fprintf('  [ct=%.2f] Trajectory %d/%d\n', ct, trajID, nSamplesPerCool);
        tic;
        multilayer_random_cooling(trajID, rootFolder, ct);
        trajTimes(trajID) = toc;
        fprintf('  [ct=%.2f] Trajectory %d finished in %.2f s\n', ct, trajID, trajTimes(trajID));
    end

    fprintf('coolTime=%.2f complete. Mean time per trajectory: %.2f s\n', ct, mean(trajTimes));
end

fprintf('\nAll cooling-time sweeps complete.\n');
