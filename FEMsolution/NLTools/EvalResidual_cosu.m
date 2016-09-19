function [Res] = EvalResidual_cosu(sys,gs4,yg)

Cee = sys.bigC;
Kay = sys.bigK;


%set ic for yd if necessary
if gs4.tn == 0
    yinit = zeros(sys.nNodes,1);
    for i = 1:sys.nNodes
        yinit(i) = sys.nodes(i).y;
    end
    exF = cos(yinit);
    ydinit = Cee \ (-exF-Kay*yinit);
    for i = 1:sys.nNodes
        sys.nodes(i).yd = ydinit(i);
    end
end

%on to the important business
ynpw1 = zeros(sys.nNodes,1);
ydnpl6w1 = zeros(sys.nNodes,1);


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
locs
ynpw1
ynpw1_func = interp1(locs,ynpw1,method,'pp');
ynpw1_func
ynpw1_func(0.5)
ynpw1(end)

%compute f at tnpw1
for i = 1:sys.nEle
    elmt = sys.ele(i);
    x1 = elmt.nodes(1).loc;
    x2 = elmt.nodes(2).loc;
    ff1 = @(x) (1-(x2-x)./elmt.dx).*x*x*cos(ynpw1_func(x));
    ff2 = @(x) (x2-x)./elmt.dx.*x*x*cos(ynpw1_func(x));
    elmt.force = [ integral(ff1,x1,x2,'ArrayValued', true);
        integral(ff2,x1,x2,'ArrayValued', true) ];
end

elmt.force

wprkdd

% set type 2 and 3 BC values in stiffness matrices

for i = 1:sys.nbc
    if sys.bcs(i).type == 2 || sys.bcs(i).type == 3
        [Cee,Kay,exF] = sys.bcs(i).applytype23(Cee,Kay,exF);
    end
end

Res = Cee*ydnpl6w1 + Kay*ynpw1 + exF;

%'enforce' dirichlet bcs
for i = 1:sys.nbc
    if sys.bcs(i).type == 1
        Res(sys.bcs(i).where) = 0;
    end
end

%Res

end
