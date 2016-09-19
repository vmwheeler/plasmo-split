function [J] = EvalJacobian_cosu(sys,gs4,yg)

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
    
% set type 2 and 3 BC values in stiffness matrices
Cee = sys.bigC;
Kay = sys.bigK;
exF = sys.force;
for i = 1:sys.nbc
    if sys.bcs(i).type == 2 || sys.bcs(i).type == 3
        [Cee,Kay,exF] = sys.bcs(i).applytype23(Cee,Kay,exF);
    end
end

%This is the contribution from the nonlinear source term...
locs = linspace(0,1,sys.nNodes);
ynpw1_func = griddedInterpolant(locs,ynpw1);
Mp = zeros(sys.nNodes,sys.nNodes);
%compute f' at tnpw1
for i = 1:sys.nEle
    elmt = sys.ele(i);
    x1 = elmt.nodes(1).loc;
    x2 = elmt.nodes(2).loc;
    ff1 = @(x) -(1-(x2-x)./elmt.dx).*x*x*sin(ynpw1_func(x));
    ff2 = @(x) -(x2-x)./elmt.dx.*x*x*sin(ynpw1_func(x));
    elmt.force = [ integral(ff1,x1,x2,'ArrayValued', true);
        integral(ff2,x1,x2,'ArrayValued', true) ];
    
    %get node numbers for element to be assembled
    nn = zeros(elmt.nnpe,1);
    for j = 1:length(nn)
        nn(j) = elmt.nodes(j).num;
    end
    
    %external force assembly
    %for m = 1:length(nn)
    %    Mp(nn(m),nn(m)) = Mp(nn(m),nn(m)) - elmt.force(m);
    %end
end



%Here is J
J = gs4.lam6w1/gs4.lam5/gs4.dt*Cee +gs4.lam5w2/gs4.lam5*(Kay + Mp);

%apply dirichlet bc
for i = 1:sys.nbc
    bc = sys.bcs(i);
    no = bc.where;
    if sys.bcs(i).type == 1
        J(no,:) = zeros(1,length(J(no,:)));
        J(no,no) = 1;   
    end
end

%J

end
