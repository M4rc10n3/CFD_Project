function [f1,f2,f3] = decod_ausm(rho,u,ener)

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIGAMMA_PAR
global h %Nuove variabili globali

%A seconda di quello che dice lo schema, possiamo modificare questa
%funzione per trovare le grandezze che ci interessano.

p = (ener - .5*rho*u^2)./gb;
h = p / rho * ga; % ricordarsi che stiamo calcolando h ma non lo usiamo
%rho=p/h*ga;
%e=p*gb+.5*rho*u^2;

f1=rho*u;
f2=p+rho*u^2;
f3=u*(p+ener);

end