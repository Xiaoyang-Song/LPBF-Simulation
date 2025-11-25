rootFolder = 'RL_Dataset';
trajectoryFolders = dir(fullfile(rootFolder,'trajectory_*'));
nTrajectories = length(trajectoryFolders);

% Initialize storage for all trajectories
allActions = {};
allMeanDevs = {};
allMaxTemps = {};
allStates = {};
allNextStates = {};

fprintf('Found %d trajectories in %s\n', nTrajectories, rootFolder);

for t = 1:nTrajectories
    trajFolder = fullfile(rootFolder, trajectoryFolders(t).name);
    fprintf('\nLoading trajectory %s\n', trajectoryFolders(t).name);
    
    % Load summary
    summaryFile = fullfile(trajFolder,'trajectory_summary.mat');
    if exist(summaryFile,'file')
        data = load(summaryFile);
        
        % Store into cell arrays
        allActions{t} = data.actions;
        allMeanDevs{t} = data.meanDeviationList;
        allMaxTemps{t} = data.maxTempList;
        allStates{t} = data.uFinalList;  % last temperature per step (can also save tAll/uAll if needed)
        
        % Build next_states: shift by one step
        nextStates = data.uFinalList(2:end);
        % pad last step with itself (optional)
        nextStates{end+1} = data.uFinalList{end};
        allNextStates{t} = nextStates;
        
        fprintf('Loaded %d layers\n', size(data.actions,1));
    else
        warning('Summary file not found in %s', trajFolder);
    end
    
    % Optionally, read console output
    consoleFile = fullfile(trajFolder,'console_output.txt');
    if exist(consoleFile,'file')
        fid = fopen(consoleFile,'r');
        consoleLog = fread(fid,'*char')';
        fclose(fid);
        fprintf('Console log preview:\n%s\n', consoleLog(1:min(500,end)));  % first 500 chars
    end
end

%% ---------------------------
% Example usage:
% allActions{1} -> actions matrix for trajectory 1
% allStates{1}{2} -> temperature field after step 2 in trajectory 1
% allNextStates{1}{2} -> next state after step 2
% allMeanDevs{1} -> rewards per step for trajectory 1

disp(size(allStates{1}{1}))

disp(size(allActions{1}))
disp(allActions{1})
disp(size(allMeanDevs{1}))
disp(allMeanDevs{1})