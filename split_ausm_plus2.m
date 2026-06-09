function split_ausm_plus2

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIgamma_PAR
global nc ncm k kout ka iord itest stab % era USE NUM_PAR
global b c diverg c1 c2 x0 dx                     % era USE GEOM_PAR
global   time dt dtodx dtitest4 timemax timek     % era USE TIME_PAR
global  xsh10 xsh20 p1 rho1 u1 p2 rho2 u2 p3 rho3 u3 p4 rho4 u4 % era USE INIT_PAR
global  vsh1 vsh2 vsh3 vsh4 vsh5 vsh3l vsh3r vsh4l vsh4r % era USE INIT_PAR
global  xsh1 xsh2 xsh3 xsh4 xsh5 xsh3l xsh3r xsh4l xsh4r % era USE INIT_PAR
global p t u s rho a e amach ptot ttot flow flht h htot % era USE VARS
global w1 w2 w3 f1 f2 f3 phi1 phi2 phi3                    % era USE VARS
global pxeno uxeno hxeno ppxeno hhxeno ppt ut hht % era USE ENO

enuo2 = 0.5*dt/dx;

%Variabili: 
% ncm: è il numero totale delle interfacce = n.celle-1 =nc-1
% ncmm: è il numero di iterazioni da compiere perché all'iterazione ncmm ci
% occupiamo anche dell'interfaccia ncmm+1=ncm, che è la sinistra
% dell'ultima cella

ncmm = ncm-1;
for n=2:ncmm
    nm = n;
    np = nm+1;
    
    pa = p(nm);
    pb = p(np);
    ua = u(nm);
    ub = u(np);
    ha = h(nm);
    hb = h(np);
    
    % Le condizioni al bordo qui sotto servono poi per i calcoli 
    % da riga 240 in poi; da capire quali implementare
    %% Condizioni al bordo di D'Ambrosio
    %{
    if n == 2
        p002   =  pa;
        u002   =  ua;
        h002   =  ha;
        rho002 =  p002/h002*ga;
        a002   =  sqrt(gamma*p002/rho002);
    end
    
    if n == ncmm
        pncm   =  pb;
        uncm   =  ub;
        hncm   =  hb;
        rhoncm =  pncm/hncm*ga;
        ancm   =  sqrt(gamma*pncm/rhoncm);
    end
    %}

    %% Condizioni al bordo nostre
    if n == 2
        p002   = pa; 
        u002   = ua;
        h002   = ha;
        rho002 = p002/h002*ga;

        % Nota: sx = dx, cioè pa = pb = p002, etc.

        %Entalpia specifica totale
        h_t002 = h002 + u002^2/2;
        
        %In questo caso non necessario calcolare entrambi i valori
        %a_starL002 = sqrt(2*(gamma-1)/(gamma+1)*h_t002);
        %a_starR002 = sqrt(2*(gamma-1)/(gamma+1)*h_t002);
        %a_tildeL002 = a_starL002^2/max(a_starL002,abs(u002)); %u002 = ua
        %a_tildeR002 = a_starR002^2/max(a_starR002,abs(ub));
        %a002 = min(a_tildeR002,a_tildeL002);

        a_star002 = sqrt(2*(gamma-1)/(gamma+1)*h_t002);
        % In questo caso ignoriamo il calcolo di a_tilde poiché 
        % a002 = a_tilde002
        a002 = a_star002^2/max(a_star002,abs(u002)); 
        M002 = u002/a002;
        
        if abs(M002) >= 1 
            M_cors_plus = 0.5*(M002+abs(M002));
            P_cors_plus = 0.5*(1+sign(M002));

            M_cors_minus = 0.5*(M002-abs(M002));
            P_cors_minus = 0.5*(1-sign(M002));
        else  % valori ottimali già settati beta=1/8, alpha=3/16
            M_cors_plus = 0.25*(M002+1)^2 + (1/8)*(M002^2-1)^2; 
            P_cors_plus = 0.25*(M002+1)^2*(2-M002) + (3/16)*M002*(M002^2-1)^2;

            M_cors_minus = -0.25*(M002-1)^2 - (1/8)*(M002^2-1)^2;
            P_cors_minus = 0.25*(M002-1)^2*(2+M002) - (3/16)*M002*(M002^2-1)^2;
        end

        m002 = M_cors_plus + M_cors_minus;
        p002 = P_cors_plus*p002 + P_cors_minus*p002;

        m_plus002 = 0.5*(m002+abs(m002));
        m_minus002 = 0.5*(m002-abs(m002));
    end
    
    if n == ncmm
        pncm   = pb;
        uncm   = ub;
        hncm   = hb;
        rhoncm = pncm/hncm*ga;

        % Nota: sx = dx, cioè pa = pb = pncm, etc.
        h_tncm = hncm + uncm^2/2;
        
        %Idem come sopra: non è necessario calcolare entrambi i valori.
        %a_starLncm = sqrt(2*(gamma-1)/(gamma+1)*h_tncm);
        %a_starRncm = sqrt(2*(gamma-1)/(gamma+1)*h_tncm);
        a_starncm = sqrt(2*(gamma-1)/(gamma+1)*h_tncm);
        %a_tildeLncm = a_starLncm^2/max(a_starLncm,abs(uncm));
        %a_tildeRncm = a_starRncm^2/max(a_starRncm,abs(uncm));
        %ancm = min(a_tildeRncm,a_tildeLncm);
        ancm = a_starncm^2/max(a_starncm,abs(uncm));
        Mncm = uncm/ancm;
        
        if abs(Mncm) >= 1 
            M_cors_plus = 0.5*(Mncm+abs(Mncm));
            P_cors_plus = 0.5*(1+sign(Mncm));

            M_cors_minus = 0.5*(Mncm-abs(Mncm));
            P_cors_minus = 0.5*(1-sign(Mncm));
        else  % valori ottimali già settati beta=1/8, alpha=3/16
            M_cors_plus = 0.25*(Mncm+1)^2 + (1/8)*(Mncm^2-1)^2;
            P_cors_plus = 0.25*(Mncm+1)^2*(2-Mncm) + (3/16)*Mncm*(Mncm^2-1)^2;

            M_cors_minus = -0.25*(Mncm-1)^2 - (1/8)*(Mncm^2-1)^2;
            P_cors_minus = 0.25*(Mncm-1)^2*(2+Mncm) - (3/16)*Mncm*(Mncm^2-1)^2;
        end

        mncm = M_cors_plus + M_cors_minus;
        pncm = P_cors_plus*pncm + P_cors_minus*pncm;

        m_plus_ncm = 0.5*(mncm+abs(mncm));
        m_minus_ncm = 0.5*(mncm-abs(mncm));
    end
    

    %% Per aumentare l'ordine dell'algoritmo
    if iord ~= 1
        ppa = log(pa);
        ppb = log(pb);
        hha = log(ha);
        hhb = log(hb);
        ppa  =  ppa + .5* ppxeno(n);
        ppb  =  ppb - .5* ppxeno(n+1);
        ua   =  ua  + .5* uxeno(n);
        ub   =  ub  - .5* uxeno(n+1);
        hha  =  hha + .5* hhxeno(n);
        hhb  =  hhb - .5* hhxeno(n+1);
        
        ppa = ppa + enuo2*ppt(nm);
        ua  = ua  + enuo2*ut(nm);
        hha = hha + enuo2*hht(nm);
        
        if n == 2
            pp002  =  log(p(nm));
            u002   =  u(nm);
            hh002  =  log(h(nm));
            pp002  =  pp002 - .5* ppxeno(n);
            u002   =  u002  - .5* uxeno(n);
            hh002  =  hh002 - .5* hhxeno(n);
            pp002  =  pp002 + enuo2*ppt(n);
            u002   =  u002  + enuo2*ut(n);
            hh002  =  hh002 + enuo2*hht(n);
            p002   =  exp(pp002);
            h002   =  exp(hh002);
            rho002 =  p002/h002*ga;
            a002   =  sqrt(gamma*p002/rho002);
        end
        
        ppb = ppb + enuo2*ppt(np);
        ub  = ub  + enuo2*ut(np);
        hhb = hhb + enuo2*hht(np);
        
        if n == ncmm
            ppncm  =  log(p(np));
            uncm   =  u(np);
            hhncm  =  log(h(np));
            ppncm  =  ppncm    + .5* ppxeno(n+1);
            uncm   =  uncm     + .5* uxeno(n+1);
            hhncm  =  hhncm    + .5* hhxeno(n+1);
            ppncm  =  ppncm    + enuo2*ppt(n+1);
            uncm   =  uncm     + enuo2*ut(n+1);
            hhncm  =  hhncm    + enuo2*hht(n+1);
            pncm   =  exp(ppncm);
            hncm   =  exp(hhncm);
            rhoncm =  pncm/hncm*ga;
            ancm   =  sqrt(gamma*pncm/rhoncm);
        end
        
        pa = exp(ppa);
        pb = exp(ppb);
        ha = exp(hha);
        hb = exp(hhb);
    end
    
    %% Tutte le interfacce che non sono la seconda o la penultima
    if n >= 2 && n <= ncmm

        % calcolo sul paper usando a per j e b per j+1
        rhoa = pa/ha*ga;
        rhob = pb/hb*ga;
    
        h_ta = ha + ua^2/2;
        h_tb = hb + ub^2/2;
        
        a_star_a = sqrt(2*(gamma-1)/(gamma+1)*h_ta);
        a_star_b = sqrt(2*(gamma-1)/(gamma+1)*h_tb);

        a_tilde_a = a_star_a^2/max(a_star_a,abs(ua));
        a_tilde_b = a_star_b^2/max(a_star_b,abs(ub));

        a_n = min(a_tilde_a,a_tilde_b);

        Ma = ua/a_n;
        Mb = ub/a_n;

        if abs(Ma) >= 1 % valori ottimali già settati beta=1/8, alpha=3/16
            M_cors_plus = 0.5*(Ma+abs(Ma));
            P_cors_plus = 0.5*(1+sign(Ma));
        else
            M_cors_plus = 0.25*(Ma+1)^2 + (1/8)*(Ma^2-1)^2;
            P_cors_plus = 0.25*(Ma+1)^2*(2-Ma) + (3/16)*Ma*(Ma^2-1)^2;
        end

        if abs(Mb) >= 1  
            M_cors_minus = 0.5*(Mb-abs(Mb));
            P_cors_minus = 0.5*(1-sign(Mb));
        else
            M_cors_minus = -0.25*(Mb-1)^2 - (1/8)*(Mb^2-1)^2;
            P_cors_minus = 0.25*(Mb-1)^2*(2+Mb) - (3/16)*Mb*(Mb^2-1)^2;
        end
    
        m_n = M_cors_plus + M_cors_minus;
        p_n = P_cors_plus*pa + P_cors_minus*pb;
        m_n_plus = 0.5*(m_n+abs(m_n));
        m_n_minus = 0.5*(m_n-abs(m_n));
    
    % Formula A3 del paper scomposta nelle 3 componenti:
    phi1(n) = a_n*(m_n_plus*rhoa + m_n_minus*rhob); 
    phi2(n) = a_n*(m_n_plus*rhoa*ua + m_n_minus*rhob*ub) + p_n;
    phi3(n) = a_n*(m_n_plus*rhoa*h_ta + m_n_minus*rhob*h_tb);
    
    end

    %{ 
    Dopo aver calcolato questi tre flussi, all'interno del file march.m si
    vanno a modificare come si trovano la quantità che ci interessano.
    Infatti, le variabili conservative dell'algoritmo AUSM+ sono diverse da quelle
    usate dal Professor D'Ambrosio per implementare il metodo di Euleuro a
    1 dimensione.
    %}

end %end of the for loop

%% Condizioni al bordo che differiscono a seconda del test svolto
%{
if(itest == 1)  %REFLECTING WALL B.C.
    %%DA RIVEDERE -> la formula (64) del paper potrebbe aiutare?
    r1dum  = p002-rho002*a002*u002;
    r2dum  = h002-p002/rho002;
    pin    = r1dum;
    uin    = 0.0;
    hin    = r2dum+pin/rho002;
    [phi1(1),phi2(1),phi3(1)] = decod_ausm(pin,uin,hin);
    r3dum  = pncm+rhoncm*ancm*uncm;
    r2dum  = hncm-pncm/rhoncm;
    pex    = r3dum;
    uex    = 0.0;
    hex    = r2dum+pex/rhoncm;
    [phi1(ncm),phi2(ncm),phi3(ncm)] = decod_ausm(pex,uex,hex);
end
%}

if(itest == 1)  %REFLECTING WALL B.C.
    %%DA RIVEDERE -> la formula (64) del paper potrebbe aiutare?    
    ain    = a002;
    pin    = p002;
    uin    = u002;
    hin    = h002;
    rhoin  = rho002;
    rho_2  = p(3)/h(3)*ga;
    h_tin  = hin + uin^2/2;
    h_t_2  = h(3) + u(3)^2/2;
    m_minus_in = m_minus002;
    m_plus_in  = m_plus002;
    phi1(1) = ain*(m_plus_in*rhoin + m_minus_in*rho_2);
    phi2(1) = ain*(m_plus_in*rhoin*uin + m_minus_in*rho_2*u(3)) + pin;
    phi3(1) = ain*(m_plus_in*rhoin*h_tin + m_minus_in*rho_2*h_t_2);

    aex    = ancm;
    pex    = pncm;
    uex    = uncm;
    hex    = hncm;
    rhoex  = rhoncm;
    rho_ncmm = p(ncmm)/h(ncmm)*ga;
    h_tex  = hex + uex^2/2;
    h_t_ncmm = h(ncmm) + u(ncmm)^2/2;
    m_minus_ex = m_minus_ncm;
    m_plus_ex  = m_plus_ncm;
    phi1(ncm) = aex*(m_plus_ex*rho_ncmm + m_minus_ex*rhoex);
    phi2(ncm) = aex*(m_plus_ex*rho_ncmm*u(ncmm) + m_minus_ex*rhoex*uex) + pex;
    phi3(ncm) = aex*(m_plus_ex*rho_ncmm*h_t_ncmm + m_minus_ex*rhoex*h_tex);
end

% Versione precedente non funzionante
%{
if(itest >= 2)  % REFLECTING WALL B.C.
    %Per noi inutili:
    %r1dum  = p002-rho002*a002*u002;
    %r2dum  = h002-p002/rho002;

    pin    = p002;
    uin    = u002;
    hin    = h002;
    [phi1(1),phi2(1),phi3(1)] = decod_ausm(pin,uin,hin);

    %Per noi inutili:
    %r3dum  = pncm+rhoncm*ancm*uncm;
    %r2dum  = hncm-pncm/rhoncm;

    pex    = pncm;
    uex    = uncm;
    hex    = hncm;
    [phi1(ncm),phi2(ncm),phi3(ncm)] = decod_ausm(pex,uex,hex);
end
%}

if(itest >= 2)  % REFLECTING WALL B.C.
    ain    = a002;
    pin    = p002;
    uin    = u002;
    hin    = h002;
    rhoin  = rho002;
    rho_2  = p(3)/h(3)*ga;
    h_tin  = hin + uin^2/2;
    h_t_2  = h(3) + u(3)^2/2;
    m_minus_in = m_minus002;
    m_plus_in  = m_plus002;
    phi1(1) = ain*(m_plus_in*rhoin + m_minus_in*rho_2);
    phi2(1) = ain*(m_plus_in*rhoin*uin + m_minus_in*rho_2*u(3)) + pin;
    phi3(1) = ain*(m_plus_in*rhoin*h_tin + m_minus_in*rho_2*h_t_2);

    aex    = ancm;
    pex    = pncm;
    uex    = uncm;
    hex    = hncm;
    rhoex  = rhoncm;
    rho_ncmm = p(ncmm)/h(ncmm)*ga;
    h_tex  = hex + uex^2/2;
    h_t_ncmm = h(ncmm) + u(ncmm)^2/2;
    m_minus_ex = m_minus_ncm;
    m_plus_ex  = m_plus_ncm;
    phi1(ncm) = aex*(m_plus_ex*rho_ncmm + m_minus_ex*rhoex);
    phi2(ncm) = aex*(m_plus_ex*rho_ncmm*u(ncmm) + m_minus_ex*rhoex*uex) + pex;
    phi3(ncm) = aex*(m_plus_ex*rho_ncmm*h_t_ncmm + m_minus_ex*rhoex*h_tex);
end

end %end of function