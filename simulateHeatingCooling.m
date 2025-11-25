function [uFinal, tAll, uAll, resultAll, model, meanDeviation] = simulateHeatingCooling(paramsStruct)
% SIMULATEHEATINGCOOLING Simulates heating and cooling of a 2D domain
%
% Inputs (inside paramsStruct):
%   .k, .rho, .specificHeat, .ta, .thick
%   .width, .height, .hmax, .squareSideFraction
%   .scan_pattern, .style, .params
%   .heatTime, .coolTime, .nTimeStepsHeat, .nTimeStepsCool
%   .doPlot      - true/false to enable plotting
%   .tempRange   - [T_l, T_h] temperature range for mean deviation calculation
%
% Outputs:
%   uFinal        - temperature distribution at the last step of cooling
%   tAll          - time vector for all steps
%   uAll          - nodal solutions at all time points
%   model         - PDE model object
%   meanDeviation - mean deviation of temperature inside the square from the range

%% ---------------------------
% Material parameters
c = paramsStruct.thick * paramsStruct.k;
d = paramsStruct.thick * paramsStruct.rho * paramsStruct.specificHeat;
a = 0.01;

%% ---------------------------
% Geometry
width = paramsStruct.width;
height = paramsStruct.height;
gdm = [3 4 0 width width 0 0 0 height height]';
dl = decsg(gdm,'S1',('S1')');

%% ---------------------------
% Mesh
hmax = paramsStruct.hmax;

%% ---------------------------
% Square scan region
cx = width/2; cy = height/2;
squareSide = min(width,height) * paramsStruct.squareSideFraction;
halfSide = squareSide/2;
sqX = [cx-halfSide, cx+halfSide, cx+halfSide, cx-halfSide];
sqY = [cy-halfSide, cy-halfSide, cy+halfSide, cy+halfSide];

scan_pattern = paramsStruct.scan_pattern;
style = paramsStruct.style;

%% ---------------------------
% Time vectors
tlistHeat = linspace(0, paramsStruct.heatTime, paramsStruct.nTimeStepsHeat);
tlistCool = linspace(0, paramsStruct.coolTime, paramsStruct.nTimeStepsCool);

%% ---------------------------
% Create PDE model
model = createpde(1);
geometryFromEdges(model, dl);
generateMesh(model,'Hmax',hmax);

%% ---------------------------
% Heat function for heating
fHeat = @(location,state) heatSource(location,state,paramsStruct.params,style,scan_pattern,sqX,sqY);
specifyCoefficients(model,'m',0,'d',d,'c',c,'a',a,'f',fHeat);
applyBoundaryCondition(model,'neumann','Edge',[1,2,3,4],'q',0,'g',0);
setInitialConditions(model, paramsStruct.ic);

%% ---------------------------
% Solve heating
resultHeat = solvepde(model, tlistHeat);
uHeat = resultHeat.NodalSolution;

%% ---------------------------
% Solve cooling
fCool = @(location,state) zeros(size(location.x));
specifyCoefficients(model,'m',0,'d',d,'c',c,'a',a,'f',fCool); 
setInitialConditions(model, resultHeat);
resultCool = solvepde(model, tlistCool);
uCool = resultCool.NodalSolution;

%% ---------------------------
% Combine results
tAll = [tlistHeat, tlistCool + paramsStruct.heatTime];
uAll = [uHeat, uCool];
resultAll = [resultHeat, resultCool];
uFinal = uCool(:,end);

%% ---------------------------
% Compute mean deviation inside the square
if isfield(paramsStruct,'tempRange')
    T_l = paramsStruct.tempRange(1);
    T_h = paramsStruct.tempRange(2);

    % Nodes coordinates
    nodes = model.Mesh.Nodes';
    inSquare = inpolygon(nodes(:,1), nodes(:,2), sqX, sqY);

    % Deviation from the range
    uInSq = uFinal(inSquare);
    deviation = zeros(size(uInSq));
    deviation(uInSq < T_l) = T_l - uInSq(uInSq < T_l);
    deviation(uInSq > T_h) = uInSq(uInSq > T_h) - T_h;
    meanDeviation = mean(deviation) / (T_h - T_l);
else
    meanDeviation = NaN;
end

%% ---------------------------
% Optional plotting
if isfield(paramsStruct,'doPlot') && paramsStruct.doPlot
    figure('Name','Temperature Evolution');
    for i = 1:length(tlistHeat)
        pdeplot(model,'XYData',uHeat(:,i),'Mesh','on','ColorMap','jet');
        colorbar; caxis([300 5000]);
        hold on; plot([sqX sqX(1)],[sqY sqY(1)],'m-','LineWidth',1.5); hold off;
        title(sprintf('Heating Phase: t = %.4f s', tlistHeat(i)));
        drawnow;
    end
    for i = 1:length(tlistCool)
        pdeplot(model,'XYData',uCool(:,i),'Mesh','on','ColorMap','jet');
        colorbar; caxis([300 5000]);
        hold on; plot([sqX sqX(1)],[sqY sqY(1)],'m-','LineWidth',1.5); hold off;
        title(sprintf('Cooling Phase: t = %.4f s', tlistCool(i)+paramsStruct.heatTime));
        drawnow;
    end

    % Global metrics
    mean_u = mean(uAll,1);
    max_u  = max(uAll,[],1);
    min_u  = min(uAll,[],1);
    figure('Name','Temperature Metrics');
    plot(tAll, mean_u,'b-','LineWidth',2); hold on;
    plot(tAll, max_u,'r--','LineWidth',2);
    plot(tAll, min_u,'g-.','LineWidth',2);
    xlabel('Time (s)'); ylabel('Temperature (K)');
    legend('Mean','Max','Min','Location','best');
    title('Temperature evolution in the domain');
    grid on;
end

%% ---------------------------
% Subfunction: Heat source
    function f = heatSource(location,state,params,style,scan_pattern,sqX,sqY)
        f = zeros(size(location.x));
        inSq = inpolygon(location.x, location.y, sqX, sqY);
        fvals = Heat_MultipleScans(location,state,params,style,scan_pattern);
        f(inSq) = fvals(inSq);
    end

end
