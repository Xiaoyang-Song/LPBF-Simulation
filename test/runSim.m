
    cd('../../LPBF-Simulation/');
    paramsStruct = load('test/params.mat').paramsStruct;
    [uFinal, tAll, uAll, resultAll, model, meanDeviation] = simulateHeatingCooling(paramsStruct);
    save('test/results.mat','uFinal','tAll','uAll','meanDeviation');
    i=7;
    fig = figure('Visible','off');
    pdeplot(model,'XYData',uFinal,'Mesh','on','ColorMap','jet');
    colorbar; caxis([300 5000]);
    title(sprintf('Step %d: Cooling Final Temperature',i));
    saveas(fig, fullfile("../Offline-RL-Controller-in-AM/checkpoints/results/",sprintf('layer_%d_finalTemp.png',i)));
    close(fig);
    exit
    