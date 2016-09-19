function [Res] = EvalResidual_easy(sys,gs4,yg)

ynpw1 = zeros(sys.nNodes,1);
ydnpl6w1 = zeros(sys.nNodes,1);

Cee = sys.bigC;
Kay = sys.bigK;
sys.updateForce(gs4.tnpw1)
exF = sys.force;


if gs4.tn == 0
    yinit = zeros(sys.nNodes,1);
    for i = 1:sys.nNodes
        yinit(i) = sys.nodes(i).y;
    end
    ydinit = Cee \ (-exF-Kay*yinit);
    for i = 1:sys.nNodes
        sys.nodes(i).yd = ydinit(i);
    end
end

for i = 1:sys.nNodes
    yn = sys.nodes(i).y;
    ydn = sys.nodes(i).yd;
    ydg = (1-1/gs4.lam5)*ydn + 1/gs4.lam5/gs4.dt*(yg(i)-yn);
    %fprintf('%f, %f\n',yg(i),yn)
    ynpw1(i) = yn + gs4.w1*(yg(i) - yn);
    ydnpl6w1(i) = ydn + gs4.lam6w1*(ydg-ydn);
end

% set type 2 and 3 BC values in stiffness matrices
for i = 1:sys.nbc
    if sys.bcs(i).type == 2 || sys.bcs(i).type == 3
        [Cee,Kay,exF] = sys.bcs(i).applytype23(Cee,Kay,exF);
    end
end

Res = Cee*ydnpl6w1 + Kay*ynpw1 + exF;

for i = 1:sys.nbc
    if sys.bcs(i).type == 1
        Res(sys.bcs(i).where) = 0;
    end
end

%Res

end
