function [Res] = EvalResidual(sys,gs4,ydg)

ynpw1 = zeros(sys.nNodes,1);
ydnpl6w1 = zeros(sys.nNodes,1);
for i = 1:sys.nNodes
    yn = sys.nodes(i).y;
    ydn = sys.nodes(i).yd;
    yg = yn + gs4.lam4*gs4.dt*ydn + gs4.lam5*gs4.dt*(ydg(i)-ydn);
    fprintf('%f, %f\n',yg,yn)
    ynpw1(i) = yn + gs4.w1*(yg - yn);
    ydnpl6w1(i) = ydn + gs4.lam6w1*(ydg(i)-ydn);
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

Res = Cee*ydnpl6w1 + Kay*ynpw1 + exF;

for i = 1:sys.nbc
    if sys.bcs(i).type == 1
        Res(sys.bcs(i).where) = 0;
    end
end

Res

end