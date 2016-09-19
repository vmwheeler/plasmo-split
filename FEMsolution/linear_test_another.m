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
numEle = 10;
numNodes = numEle+1;
rhoMax = 1.; rhoMin = 0.0; rhoEss = 0.0;
tEnd = 50;
numSteps = 10;
dt=tEnd/numSteps;

gs4 = GS4(rhoMax,rhoMin,rhoEss,dt,1);

%% Initialize vectors
problemIC = zeros(numNodes,1);
xPlot = linspace(0,1,numNodes);
sol = zeros(numNodes,1);
xFlux = zeros(numNodes-1,1);
flux = zeros(numNodes-1,1);
fluxBD = zeros(numNodes-1,1);


%% Import (or create) mesh data

% generate nodal locations
nodeLocs = linspace(0.1,1,numNodes);

% initialize (set ICs) and number nodes
for i = 1:numNodes
    iv = problemIC(i);
    nodes(i) = Node(i,nodeLocs(i),iv,0);
end


%% Create elements from mesh
% set the force term 
fhandle = @force_term_0;

% assign nodes to elements
for i = 1:numEle
    eleNodes = [nodes(i);nodes(i+1)];
    ele(i) = Sphere_1DL(i,eleNodes, [1,1,1,1], fhandle);
    %ele(i) = Element1DL(i,eleNodes, [1,1], fhandle);
end

%% Set up system of equations
% initiallize
sysEQ = SystemEQ(nodes,ele);
%sysEQ.dynamic_force = true;

%% set BCs
BC1 = BoundaryCondition(1,2,1,0,0.0,0.0);
sysEQ.addBC(BC1);
BC2 = BoundaryCondition(2,3,numNodes,0,1.0,1.0);
sysEQ.addBC(BC2);



%% Solve!
sysEQ.ready()
sysEQ.bigK
for n = 1:numSteps
    fprintf('**********\nTimestep #%i\n**********\n',n)
    gs4.time_march(sysEQ);
end


%% Post-process

for j = 1:numNodes
    sol(j) = sysEQ.nodes(j).y;
end



sysEQ.computeDirs()
for i = 1:numEle
    xFlux(i) = (sysEQ.ele(i).nodes(2).loc + sysEQ.ele(i).nodes(1).loc)/2;
    flux(i) = sysEQ.ele(i).yp;
end

hT = figure();
plot(xPlot,sol)
set(gca,'fontsize',8)
ylabel('dimensionless temperature')
xlabel('dimensionless distance')

    
hQ = figure();
plot(xFlux,flux)
set(gca,'fontsize',8)
ylabel('dimensionless flux')
xlabel('dimensionless distance')



