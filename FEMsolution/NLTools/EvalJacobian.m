function [J] = EvalJacobian(sys,gs4)

for i = 1:sys.nNodes
    yn(i) = sys.nodes(i).y;
end
    
% set type 2 and 3 BC values in stiffness matrices
Cee = sys.bigC;
Kay = sys.bigK;
exF = sys.force;
for i = 1:sys.nbc
    if sys.bcs(i).type == 2 || sys.bcs(i) == 3
        [Cee,Kay,exF] = sys.bcs(i).applytype23(Cee,Kay,exF);
    end
end

%This is the contribution from the nonlinear source term...
%for now it is zero
Mp = zeros(sys.nNodes,sys.nNodes);

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

J

end
