
                    cd('../../LPBF-Simulation/');
                    paramsStruct = load('test/params.mat').paramsStruct;
                    if 8 > 0
                        prevData = load('test/results.mat', 'resultCool');
                        paramsStruct.ic = prevData.resultCool;
                    end
                    [uFinal, tAll, uAll, resultAll, model, meanDeviation] = simulateHeatingCooling(paramsStruct);
                    resultCool = resultAll(2);
                    save('test/results.mat','uFinal','tAll','uAll','meanDeviation','resultCool');
                    i=8;
                    fig = figure('Visible','off');
                    pdeplot(model,'XYData',uFinal,'Mesh','on','ColorMap','jet');
                    colorbar; caxis([300 5000]);
                    title(sprintf('Step %d: Cooling Final Temperature',i));
                    saveas(fig, fullfile("../Offline-RL-Controller-in-AM/checkpoints/results/",sprintf('layer_%d_finalTemp.png',i)));
                    close(fig);
                    exit
    