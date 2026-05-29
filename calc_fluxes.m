function [f1,f2,f3] = calc_fluxes(rr,uu,pp,ee)

f1 = rr.*uu;
f2 = pp+rr.*uu.^2;
f3 = uu.*(pp+ee);

end