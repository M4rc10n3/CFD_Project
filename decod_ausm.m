function [phi1,phi2,phi3] = decod_ausm(p,u,h)

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIGAMMA_PAR

%A seconda di quello che dice lo schema, possiamo modificare questa
%funzione per trovare le grandezze che ci interessano.

rho=p/h*ga;
h_t = h + u^2/2;

phi1=rho*u;
phi2=p+rho*u^2;
phi3=rho*u*h_t;
end