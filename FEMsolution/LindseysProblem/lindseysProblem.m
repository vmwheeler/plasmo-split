clear;
close all;
clc; 
addpath('/home/vmwheeler/Code/chapterxx_r/code/Base',...
        '/home/vmwheeler/Code/chapterxx_r/code/Elements',...
        '/home/vmwheeler/Code/chapterxx_r/code/ForceTerms/',...
        '/home/vmwheeler/Code/chapterxx_r/code/Extras');

%% Physical and numerical constants
numEle = 29;
numNodes = numEle+1;
Beta = 1000.;
%Beta = 1000;
RR = 0.0025;
rhoMax = 1.; rhoMin = 0.0; rhoEss = 0.0;
tEnd = 10000.10;
numSteps = 10;
dt=tEnd/numSteps;

gs4 = GS4(rhoMax,rhoMin,rhoEss,dt,2);

%% Initialize vectors
problemIC = 0*ones(numNodes,1);
sol = zeros(numNodes,1);
flux = zeros(numNodes-1,1);


%% Import (or create) mesh data

% generate nodal locations
nodeLocs = linspace(0,RR,numNodes);

% initialize (set ICs) and number nodes
for i = 1:numNodes
    iv = problemIC(i);
    nodes(i) = Node(i,nodeLocs(i),iv,0);
end

%% Create elements from mesh
% set the force term 
fhandle = @(x,t)force_term_0(x,t);

% assign nodes to elements
for i = 1:numEle
    eleNodes = [nodes(i);nodes(i+1)];
    ele(i) = p1_1DL(i,eleNodes, [1,0.00001,3*Beta^2,1], fhandle);
end

% initiallize
sysEQ = SystemEQ(nodes,ele);
%sysEQ.dynamic_force = true;

% construct system
for q = 1:numEle
    sysEQ.addElement(ele(q));
end


q = 1.;
%% set BCs
BC1 = BoundaryCondition(1,2,1,0,0,0.0);
sysEQ.addBC(BC1);
BC2 = BoundaryCondition(2,3,numNodes,0,3/2*Beta*RR^2,4*3/2*Beta*RR^2*q);
%BC2 = BoundaryCondition(2,1,numNodes,0,0,1.0);
sysEQ.addBC(BC2);

sysEQ.ready()
%% Solve!
for n = 1:numSteps
    fprintf('**********\nTimestep #%i\n**********\n',n)
    gs4.time_march(sysEQ);
end


%% Post-process
xPlot = linspace(0,RR,numNodes);
for j = 1:numNodes
    x = nodes(j).loc;
    sol(j) = nodes(j).y;
end


for i = 1:numEle
    xFlux(i) = (sysEQ.ele(i).nodes(2).loc + sysEQ.ele(i).nodes(1).loc)/2;
    flux(i) = 1;
end



%% Show off your work
figure(3)
plot(xPlot,sol,'-rs')
ylabel('G')
xlabel('r')
%ylim([0,1])
%xlim([0,1])

% figure(2)
% plot(xFlux,flux,'-ob')
% ylabel('dimensionless flux')
% xlabel('dimensionless distance')
% xlim([0,1])
% ylim([-0.1,0.3])