function [Res] = EvalResidual_cosu(sys,gs4,yg)

Cee = sys.bigC;
Kay = sys.bigK;

fprintf('tn = %f\n', gs4.tn)
%set ic for yd if necessary
if gs4.tn == 0
    fprintf('===================\nSetting consistent IC for yd\n===================\n')
    yinit = zeros(sys.nNodes,1);
    for i = 1:sys.nNodes
        yinit(i) = sys.nodes(i).y;
    end
    %exF = cos(yinit);
    exF = zeros(sys.nNodes,1);
    ydinit = Cee \ (-exF-Kay*yinit);
    for i = 1:sys.nNodes
        sys.nodes(i).yd = ydinit(i);
    end
end

%on to the important business

%allocate space for solution vectors
ynpw1 = zeros(sys.nNodes,1);
ydnpl6w1 = zeros(sys.nNodes,1);

% set npw1 and npl6w1 vals
for i = 1:sys.nNodes
    yn = sys.nodes(i).y;
    ydn = sys.nodes(i).yd;
    ydg = (1-1/gs4.lam5)*ydn + 1/gs4.lam5/gs4.dt*(yg(i)-yn);
    %fprintf('%f, %f\n',yg(i),yn)
    ynpw1(i) = yn + gs4.w1*(yg(i) - yn);
    ydnpl6w1(i) = ydn + gs4.lam6w1*(ydg-ydn);
end

%do a linear inpterpolation of ynpw1 to get a function that can be
% integrated
locs = linspace(0,1,sys.nNodes);
ynpw1_func = griddedInterpolant(locs,ynpw1);

%ynpw1_func(locs)
sys.force = zeros(sys.nNodes,1);

%compute f at tnpw1
for i = 1:sys.nEle
    elmt = sys.ele(i);
    x1 = elmt.nodes(1).loc;
    x2 = elmt.nodes(2).loc;
    ff1 = @(x) (1-(x2-x)./elmt.dx).*x*x*cos(ynpw1_func(x));
    ff2 = @(x) (x2-x)./elmt.dx.*x*x*cos(ynpw1_func(x));
    elmt.force = [ integral(ff1,x1,x2,'ArrayValued', true);
        integral(ff2,x1,x2,'ArrayValued', true) ];
    
    %get node numbers for element to be assembled
    nn = zeros(elmt.nnpe,1);
    for j = 1:length(nn)
        nn(j) = elmt.nodes(j).num;
    end
    
    %external force assembly
    for m = 1:length(nn)
        sys.force(nn(m)) = sys.force(nn(m)) + elmt.force(m);
    end
end

exF = sys.force;

% set type 2 and 3 BC values in stiffness matrices

for i = 1:sys.nbc
    if sys.bcs(i).type == 2 || sys.bcs(i).type == 3
        [Cee,Kay,exF] = sys.bcs(i).applytype23(Cee,Kay,exF);
    end
end

Res = Cee*ydnpl6w1 + Kay*ynpw1 - exF;

%'enforce' dirichlet bcs
for i = 1:sys.nbc
    if sys.bcs(i).type == 1
        Res(sys.bcs(i).where) = 0;
    end
end

%Res

end
