function [aa,mm,ss,ee,hh,p0,t0] = calc_other_vars(pp,rr,tt,uu)

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIGAMMA_PAR

%Tutte le variabili sono normalizzate rispetto alla velocità di riferimento
aa = gf*sqrt(tt); %Velocità del suono normalizzata
mm = uu./aa; %Numero di Mach
ss = log(pp)-gamma*log(rr); %Entropia
ee = gb*pp+.5d0.*rr.*uu.^2; %Energia interna
hh = ga*tt; %Entalpia
p0 = pp.*(1.+gd.*mm.^2).^ga; %Pressione totale
t0 = tt.*(1.+gd.*mm.^2); %Temperatura totale

end