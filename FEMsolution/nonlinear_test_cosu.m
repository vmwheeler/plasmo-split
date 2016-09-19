clear all;
%close all;
clc;
pathtogs4 = '/home/vmwheeler/Research/Writings/chapterxx_r/code';
%pathtogs4 = '/home/vmwheeler/Code/chapterxx_r/code';
addpath(strcat(pathtogs4,'/Base'));
addpath(strcat(pathtogs4,'/Elements'));
addpath(strcat(pathtogs4,'/ForceTerms/'));
addpath(strcat(pathtogs4,'/Extras'));
addpath('./NLTools')

%% Physical and numerical constants
numEle = 30;
numNodes = numEle+1;
rhoMax = 1.; rhoMin = 0.0; rhoEss = 0.0;
tEnd = 0.3;
numSteps = 10;
dt=tEnd/numSteps;

gs4 = GS4(rhoMax,rhoMin,rhoEss,dt,1);

%nonlinear tolerance
tol = 10e-6;


%% Initialize vectors
problemIC = zeros(numNodes,1);
xPlot = linspace(0,1,numNodes);
sol = zeros(numNodes,1);
xFlux = zeros(numNodes-1,1);
flux = zeros(numNodes-1,1);
fluxBD = zeros(numNodes-1,1);


%% Import (or create) mesh data

% generate nodal locations
nodeLocs = linspace(0,1,numNodes);

% initialize (set ICs) and number nodes
% also set guess for first timestep
yg = zeros(numNodes,1);
for i = 1:numNodes
    iv = problemIC(i);
    yg(i) = iv;
    nodes(i) = Node(i,nodeLocs(i),iv,0);
end


%% Create elements from mesh
% set the force term 
fhandle = @(x,t) 0;

% assign nodes to elements
for i = 1:numEle
    eleNodes = [nodes(i);nodes(i+1)];
    ele(i) = Sphere_1DL(i,eleNodes, [1,1,1,1], fhandle);
end

%% Set up system of equations
% initiallize
sysEQ = SystemEQ(nodes,ele);
sysEQ.dynamic_force = true;

%% set BCs
BC1 = BoundaryCondition(1,2,1,0,0.0,0.0);
sysEQ.addBC(BC1);
BC2 = BoundaryCondition(2,3,numNodes,0,1.0,0.0);
sysEQ.addBC(BC2);


%% Solve!
sysEQ.ready()
for n = 1:numSteps
    %manually crank gs4 forward
    gs4.tick()
    eps = 299999; %reset norm of residual
    ct = 1;
    while eps > tol
        res = EvalResidual_cosu(sysEQ,gs4,yg);
        J = EvalJacobian_cosu(sysEQ,gs4,yg);
        eps = norm(res);
        fprintf('************************\n')
        fprintf('iteration # %i\n', ct)
        fprintf('norm of residual = %f\n', eps)
        fprintf('************************\n')
        delta = - J \ res;
        yg = yg + delta;
        ct = ct + 1;
        if ct > 200
            fprintf('******\n iteration max met\n******\n')
            break
        end
    end
    fprintf('\n^^^^^^^^^^^^^^^^^^^^^^^^^\n')
    for i = 1:sysEQ.nNodes
        dely = yg(i)-sysEQ.nodes(i).y;
        sysEQ.nodes(i).y = yg(i);
        sysEQ.nodes(i).yd = (1-1/gs4.lam5)*sysEQ.nodes(i).yd + 1/gs4.lam5/gs4.dt*dely;
    end
end


%% Post-process

sol = yg

%sysEQ.computeDirs()
%for i = 1:numEle
%    xFlux(i) = (sysEQ.ele(i).nodes(2).loc + sysEQ.ele(i).nodes(1).loc)/2;
%    flux(i) = sysEQ.ele(i).yp;
%end

hT = figure();
plot(xPlot,sol)
set(gca,'fontsize',8)
ylabel('dimensionless temperature')
xlabel('dimensionless distance')

    
%hQ = figure();
%plot(xFlux,flux)
%set(gca,'fontsize',8)
%ylabel('dimensionless flux')
%xlabel('dimensionless distance')
