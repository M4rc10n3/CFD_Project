function eno1
%ENO1 Summary of this function goes here
%   Detailed explanation goes here

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIgamma_PAR
global nc ncm k kout ka iord itest stab  % era USE NUM_PAR
global x                                          % era USE GEOM_VAR
global p t u s rho a e amach ptot ttot flow flht h htot % era USE VARS
global w1 w2 w3 f1 f2 f3 phi1 phi2 phi3                    % era USE VARS
global p t u s rho a e amach ptot ttot flow flht h htot % era USE VARS
global w1 w2 w3 f1 f2 f3 phi1 phi2 phi3                    % era USE VARS
global pxeno uxeno hxeno ppxeno hhxeno ppt ut hht  % era USE ENO

ncmm =nc-2;
ncmmm=nc-3;

for n=3:ncmm
    nm1=n-1; %Indice precedente
    n00=n; %Indice attuale
    np1=n+1; %Indice successivo
    
    %       ppxp=log(p(np1))-log(p(n00))
    %       ppxm=log(p(n00))-log(p(nm1))
    %       CALL minmod(ppxp,ppxm,ppxdum)
    
    %       uxp=u(np1)-u(n00)
    %       uxm=u(n00)-u(nm1)
    %       CALL minmod(uxp,uxm,uxdum)
    
    %       hhxp=log(h(np1))-log(h(n00))
    %       hhxm=log(h(n00))-log(h(nm1))
    %       CALL minmod(hhxp,hhxm,hhxdum)
    
    %n00 punto centrale; np1 punto a destra; nm1 punto a sinistra di n00
    goa    = gamma/a(n00);
    r1xp   = (log(p(np1))-goa*u(np1))-(log(p(n00))-goa*u(n00)); %p sta per "più" (+) ->
    %valori in np1 meno valori in n00

    r1xm   = (log(p(n00))-goa*u(n00))-(log(p(nm1))-goa*u(nm1)); %m sta per "meno" (-) ->
    %valori in n00 meno valori in nm1

    r1xdum = minmod(r1xp,r1xm); %operatore limite
    
    r2xp   = (log(h(np1))-gj*log(p(np1)))-(log(h(n00))-gj*log(p(n00)));
    r2xm   = (log(h(n00))-gj*log(p(n00)))-(log(h(nm1))-gj*log(p(nm1)));
    r2xdum = minmod(r2xp,r2xm);
    
    r3xp   = (log(p(np1))+goa*u(np1))-(log(p(n00))+goa*u(n00));
    r3xm   = (log(p(n00))+goa*u(n00))-(log(p(nm1))+goa*u(nm1));
    r3xdum = minmod(r3xp,r3xm);
    
    %Pendenze
    ppxdum  = .5*(r3xdum+r1xdum); %del logaritmo della pressione
    uxdum   = .5*(r3xdum-r1xdum)/goa; %della velocità
    hhxdum  = r2xdum+gj*ppxdum; %del logaritmo dell'entalpia
    
%   other possibility  - limiting the primitive variables directly
%    ppxp = log(p(np1))-log(p(n00));
%    ppxm = log(p(n00))-log(p(nm1));
%    ppxdum = minmod(ppxp,ppxm);

%    uxp = u(np1)-u(n00);
%    uxm = u(n00)-u(nm1);
%    uxdum = minmod(uxp,uxm);

%    hhxp = log(h(np1))-log(h(n00));
%    hhxm = log(h(n00))-log(h(nm1));
%    hhxdum = minmod(hhxp,hhxm);

    ppxeno(n) =  ppxdum;
    uxeno(n)  =   uxdum;
    hhxeno(n) =  hhxdum;    
end

n = 2;
ppxeno(n)= log(p(3))-log(p(2));
uxeno(n) = u(3)-u(2);
hhxeno(n)= log(h(3))-log(h(2));
n = ncm;
ppxeno(n)= log(p(ncm))-log(p(ncmm));
uxeno(n) = u(ncm)-u(ncmm);
hhxeno(n)= log(h(ncm))-log(h(ncmm));

% Le derivate temporali vengono calcolate a partire dalle pendenze, che sono 
% le derivate spaziali.
for n=2:ncm
    ppt(n) = -(u(n)*ppxeno(n)+gamma*uxeno(n));
    ut(n)  = -(u(n)*uxeno(n)+t(n)*ppxeno(n));
    hht(n) =   ppt(n)/ga - u(n)*(hhxeno(n)-ppxeno(n)/ga);
end

%Pezzo di codice che non serve a nulla 

% for n=2:ncm
%    pptt(n) = 0.0;
%    utt(n)  = 0.0;
%    hhtt(n) = 0.0;
%end

end

