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

Res=sys.bigC*ydnpl6w1 + sys.bigK*ynpw1 + sys.force;


Res(sys.nNodes) = 0;

Res

end