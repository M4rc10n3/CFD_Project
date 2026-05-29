function [f1,f2,f3] = decod(p,u,h)

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIGAMMA_PAR

%A seconda di quello che dice lo schema, possiamo modificare questa
%funzione per trovare le grandezze che ci interessano.

rho=p/h*ga;
e=p*gb+.5*rho*u^2;

f1=rho*u;
f2=p+rho*u^2;
f3=u*(p+e);

end