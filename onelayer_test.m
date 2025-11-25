clear; close all;

%% ---------------------------
%% Material properties
k = 0.0067;           % W/(mm-K)
rho = 0.00433;        % g/mm^3
specificHeat = 0.526; % J/(g-K)
ta = 300;             % ambient temperature
thick = 2.0;          % thickness

c = thick * k;
d = thick * rho * specificHeat;
a = 0.01;

%% ---------------------------
%% Geometry
width = 12; height = 3;
gdm = [3 4 0 width width 0 0 0 height height]';
dl = decsg(gdm,'S1',('S1')');

%% ---------------------------
%% Mesh
hmax = 0.2;

%% ---------------------------
%% Square scan region
cx = width/2; cy = height/2;
squareSide = min(width,height)*0.4;
halfSide = squareSide/2;
sqX = [cx-halfSide, cx+halfSide, cx+halfSide, cx-halfSide];
sqY = [cy-halfSide, cy-halfSide, cy+halfSide, cy+halfSide];

scan_pattern = linspace(0,width,48);
style = "simultaneous";

params.SS   = 600;   
params.LP   = 100;   
params.eeta = 0.3;  
params.r_b  = 0.06;  
params.H    = 0.1;  

%% ---------------------------
%% Time parameters
heatTime = 0.05;   % heating phase
coolTime = 0.10;   % cooling phase
tlistHeat = linspace(0, heatTime, 50);
tlistCool = linspace(0, coolTime, 50);

%% ---------------------------
%% Create PDE model
model = createpde(1);
geometryFromEdges(model, dl);
generateMesh(model,'Hmax',hmax);

%% ---------------------------
%% Heat function for heating
fHeat = @(location,state) heatSource(location,state,params,style,scan_pattern,sqX,sqY);
specifyCoefficients(model,'m',0,'d',d,'c',c,'a',a,'f',fHeat);
applyBoundaryCondition(model,'neumann','Edge',[1,2,3,4],'q',0,'g',0);
setInitialConditions(model, ta);

%% ---------------------------
%% Solve heating
resultHeat = solvepde(model, tlistHeat);
uHeat = resultHeat.NodalSolution;

%% ---------------------------
%% Solve cooling
fCool = @(location,state) zeros(size(location.x));
specifyCoefficients(model,'m',0,'d',d,'c',c,'a',a,'f',fCool);      % set heat to zero
setInitialConditions(model, resultHeat);  % last heating state as IC
resultCool = solvepde(model, tlistCool);
uCool = resultCool.NodalSolution;

%% ---------------------------
%% Visualization
figure;
for i = 1:length(tlistHeat)
    pdeplot(model,'XYData',uHeat(:,i),'Mesh','on','ColorMap','jet');
    colorbar; caxis([300 5000]);
    hold on; plot([sqX sqX(1)],[sqY sqY(1)],'k-','LineWidth',1.5); hold off;
    title(sprintf('Heating: t = %.4f s', tlistHeat(i)));
    drawnow;
end

for i = 1:length(tlistCool)
    pdeplot(model,'XYData',uCool(:,i),'Mesh','on','ColorMap','jet');
    colorbar; caxis([300 5000]);
    hold on; plot([sqX sqX(1)],[sqY sqY(1)],'k-','LineWidth',1.5); hold off;
    title(sprintf('Cooling: t = %.4f s', tlistCool(i)+heatTime));
    drawnow;
end

%% ---------------------------
%% Global metrics
tAll = [tlistHeat, tlistCool+heatTime];
uAll = [uHeat, uCool];
mean_u = mean(uAll,1);
max_u  = max(uAll,[],1);
min_u  = min(uAll,[],1);

figure;
plot(tAll, mean_u,'b-','LineWidth',2); hold on;
plot(tAll, max_u,'r--','LineWidth',2);
plot(tAll, min_u,'g-.','LineWidth',2);
xlabel('Time (s)'); ylabel('Temperature (K)');
legend('Mean','Max','Min','Location','best');
title('Temperature evolution in the domain');
grid on;

%% ---------------------------
%% Heat source function
function f = heatSource(location,state,params,style,scan_pattern,sqX,sqY)
    f = zeros(size(location.x));
    % Apply heat only inside square
    inSq = inpolygon(location.x, location.y, sqX, sqY);
    fvals = Heat_MultipleScans(location,state,params,style,scan_pattern);
    f(inSq) = fvals(inSq);
end
