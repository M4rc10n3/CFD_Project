function march

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIgamma_PAR
global nc ncm k kout ka iord itest stab % era USE NUM_PAR
global b c diverg c1 c2 x0 dx                     % era USE GEOM_PAR
global x                                          % era USE GEOM_VAR
global   time dt dtodx dtitest4 timemax timek    % era USE TIME_PAR
global p t u s rho a e amach ptot ttot flow flht h htot % era USE VARS
global w1 w2 w3 f1 f2 f3 phi1 phi2 phi3                    % era USE VARS
global nrt rtrms rtmax                    % era USE CONVERG

%Questa funzione non è da modificare

dtodx=dt/dx;

rtrms(1:3)=0.;
rtmax(1:3)=0.;
nrt(1:3)= 0;
nin  = 2;

for n=nin:ncm

    enu=dtodx;

    w1old = w1(n);
    w2old = w2(n);
    w3old = w3(n);

    %"phi" indica il flusso sulle diverse celle; è stato trovato in split()
    w1(n) = w1(n) - enu*(phi1(n)-phi1(n-1));
    w2(n) = w2(n) - enu*(phi2(n)-phi2(n-1));
    w3(n) = w3(n) - enu*(phi3(n)-phi3(n-1));
    
    rt(1) = (w1(n)-w1old);
    rt(2) = (w2(n)-w2old);
    rt(3) = (w3(n)-w3old);
    rtrms(1:3) = rtrms(1:3) + rt(1:3).^2;

    if (abs(rt(1)) >= abs(rtmax(1)))
        rtmax(1) = abs(rt(1));
        nrt(1) = n;
    end
    if (abs(rt(2)) >= abs(rtmax(2)))
        rtmax(2) = abs(rt(2));
        nrt(2) = n;
    end
    if (abs(rt(3)) >= abs(rtmax(3)))
        rtmax(3) = abs(rt(3));
        nrt(3) = n;
    end
end

rtrms(1:3) = sqrt(rtrms(1:3)/nc);

for n=nin:ncm
    %Probabilmente dobbiamo cambiare come si trovano queste variabili, 
    % perché nel nostro algoritmo w1, w2 e w3 non sono più queste cose
    rho(n)   = w1(n); %Questa giusta
    u(n)     = w2(n)/w1(n); %Anche questa giusta
    e(n)     = w3(n); %Questa è ERRATA, da sostituire con
    %e(n)     = w3(n)/w1(n);
    p(n)     = (e(n)-.5d0*rho(n)*u(n)*u(n))/gb;
    
    if (p(n) < 0.0)
        fprintf('pressione negativa in n==%i - p(n)=%f - p(n-1)=%f',n,time,p(n),p(n-1))
        
    end
    
    t(n)     = p(n)/rho(n);
    a(n)     = gf*sqrt(t(n));
    h(n)     = ga*t(n);
    amach(n) = u(n)/a(n);
    s(n)     = log(p(n))-gamma*log(rho(n));
    ptot(n)  = p(n)*(1.d0+gd*amach(n)^2)^ga;
    ttot(n)  = t(n)*(1.d0+gd*amach(n)^2);
    htot(n)  = h(n)+.5d0*u(n)^2;
    f1(n)    = rho(n)*u(n);
    f2(n)    = p(n)+rho(n)*u(n)*u(n);
    f3(n)    = u(n)*(p(n)+e(n));
    flow(n)  = rho(n)*u(n);
    flht(n)  = u(n)*(p(n)+e(n));
end
end