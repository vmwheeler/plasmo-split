clear all;
close all;
clc;
pathtogs4 = '/home/vmwheeler/Research/Writings/chapterxx_r/code';
%pathtogs4 = '/home/vmwheeler/Code/chapterxx_r/code';
addpath(strcat(pathtogs4,'/Base'));
addpath(strcat(pathtogs4,'/Elements'));
addpath(strcat(pathtogs4,'/ForceTerms/'));
addpath(strcat(pathtogs4,'/Extras'));
addpath('./NLTools')

%% Physical and numerical constants
numEle = 20;
numNodes = numEle+1;
rhoMax = 1.; rhoMin = 0.0; rhoEss = 0.0;
tEnd = 100.0;
numSteps = 10;
% pick a tstar to nondimensionalize time and make fix numerical issues
% due to really really small numbers
tstar = 1.E-7;
dt = tEnd/numSteps;

% radii in nanometers
crad = 30.;
crad_nm = crad*1.E-9;
srad = 60.;
srad_nm = srad*1.E-9;

% Concentration ratio
CR = 50;

% Ceria stuff
% thermal conductivity in W/mK taken as low estimate from Khafisov.  Note
% that the nonlinear curve as a function of temp is given... incorporate
% later!
k_CeO2 = 6.0;
% molar mass 
M_CeO2 = 172.115; %grams per mol
M_CeO2 = M_CeO2 / 1E3; %kg per mol

% density (wikipedia... find a better source)
rho_CeO2 = 7.65; %grams per cubic cm
rho_CeO2 = rho_CeO2 / 1E3 * (1E2)^3; %kg per cubic m

% specific heat at constant pressure from Ricken1984.  Note that cp is a
% very non-trivial function of temp.  we are just picking a reasonable one
% for now
cp_CeO2 = 80.; %J per mol per Kelvin

rhocp_CeO2 = rho_CeO2/M_CeO2*cp_CeO2;

% Gold stuff... all taken from Incropera and Dewitt pg900 table A.1
rho_Au = 19300; %kg per cubic m
cp_Au = 129.; %K per kg per K
k_Au = 317.; %W per m per K

rhocp_Au = rho_Au*cp_Au;

% lets use the result from Feng instead of a traditional Nusselt number
% correlation.
lambda = 22.;
%h=0

%ambient fluid temp
%note that this should be consistent with the temp kf was chosen at
Tinf = 300; % K

gs4 = GS4(rhoMax,rhoMin,rhoEss,dt,1);

%nonlinear tolerance
tol = 10e-6;

%import the necessary data...
pathtoinputdata = '../bhfield/dataout/';

% get absorbed power and radii
uabsfname = strcat(pathtoinputdata,num2str(crad),'-',num2str(srad),'_','UabsVr.dat');
uabsdataraw = importdata(uabsfname,' ',3);
rads = uabsdataraw.data(:,1);
pabs = uabsdataraw.data(:,2);

% get emitted power
uemfname = strcat(pathtoinputdata,num2str(crad),'-',num2str(srad),'_','UemVr.dat');
uemdataraw = importdata(uemfname,' ',4);
pem = uemdataraw.data(:,2:end);
% get temperatures
testtemp = strsplit(uemdataraw.textdata{3});
temps = transpose(cellfun(@str2num,testtemp(2:end)));

% get emitted power dT
uemdTfname = strcat(pathtoinputdata,num2str(crad),'-',num2str(srad),'_','UemdTVr.dat');
uemdTdataraw = importdata(uemdTfname,' ',4);
pemdT = uemdTdataraw.data(:,2:end);

% make some interpolation functions
[rr,tt] = ndgrid(rads,temps);
pem_func = griddedInterpolant(rr,tt,pem);
pemdT_func = griddedInterpolant(rr,tt,pemdT);
pabs_func = griddedInterpolant(rads,pabs);
%pnet_func = @(r,t) pabs_func(r) - pem_func(r,t);

%{
% cute picture of the emission vs radius and temperature to verify
% interpolation function
[xx,yy] = meshgrid(temps,rads);
figure()
surface(xx,yy,pem)
%surface(xx,yy,pem_func(xx,yy))

% absorption too
figure()
plot(rads,pabs)
%CHECK
%}



%% Initialize vectors

xPlot = linspace(0,1,numNodes);
xFlux = zeros(numNodes-1,1);
flux = zeros(numNodes-1,1);
fluxBD = zeros(numNodes-1,1);


%% Import (or create) mesh data

% generate nodal locations
nodeLocs = linspace(0,srad_nm,numNodes);

%{
%check some orders of magnitude:
% transient term
(srad_nm)^2*rhocp_CeO2
% one chunk of the diffusion term
(srad_nm)^2*k_CeO2
% the other bigger chunk
2*(srad_nm)*k_CeO2
% the absorption term (yuge)
(srad_nm)^2*pabs_func(20*1E-9)
% emission term at 1000K
(srad_nm)^2*pem_func(20*1E-9,1000)
% convection loss at boundary assuming particle at 1000K and fluid at 100K
(srad_nm)^2*h*(1000-300)
%}

% initialize (set ICs) and number nodes
% also set guess for first timestep
initial_temp = 300;

problemIC = initial_temp*ones(numNodes,1);
yg = problemIC;
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
    ele(i) = CoreShell_1DL(i,eleNodes, [crad_nm,rhocp_Au/tstar,rhocp_CeO2/tstar,k_Au,k_CeO2], fhandle);
end

%% Set up system of equations
% initiallize
sysEQ = SystemEQ(nodes,ele);
sysEQ.dynamic_force = true;


%% set BCs
BC1 = BoundaryCondition(1,2,1,0,0.0,0.0);
sysEQ.addBC(BC1);
BC2 = BoundaryCondition(2,3,numNodes,0, srad_nm^2*h, srad_nm^2*h * Tinf);
sysEQ.addBC(BC2);

%first make a really crude jacobian approximation and don't change it
% I think this will have to be refined later
J = EvalJacobian(sysEQ,gs4);

nsnaps = 5;
times = linspace(0,tEnd,nsnaps);
ns = round(times/dt);
ns(1) = 1;
tids = ns*dt*tstar;
nct = 1;

%% Solve!
sysEQ.ready()
for n = 1:numSteps
    %manually crank gs4 forward
    gs4.tick()
    eps = 299999; %reset norm of residual
    ct = 1;
    while eps > tol
        res = EvalResidual_nls(sysEQ,gs4,yg,srad_nm,CR,pabs_func,pem_func);
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

    if n == ns(nct)
        disp('plucking solution')
        for i = 1:sysEQ.nNodes
            sol(i,nct) = yg(i);
        end
        nct = nct + 1;
    end
end


%% Post-process

hT = figure(3);
set(gca,'fontsize',8)
ylabel('temperature')
xlabel('distance')

    
%% Get EPRT solutions for comparison and plot each snapshot
for i = 1:length(ns)
    set(0,'CurrentFigure',hT)
    hold on;
    l(i)=plot(xPlot,sol(:,i), '-s', 'DisplayName',strcat('t=',num2str(tids(i),'%3.2e')));
end

set(0,'CurrentFigure',hT)
box on
legend(l)
 