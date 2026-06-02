function split_ausm_plus   % QUESTA (FORSE) VA MODIFICATA PER SCRIVERE IL CODICE

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIgamma_PAR
global nc ncm k kout ka iord itest stab % era USE NUM_PAR
global b c diverg c1 c2 x0 dx                     % era USE GEOM_PAR
global   time dt dtodx dtitest4 timemax timek     % era USE TIME_PAR
global  xsh10 xsh20 p1 rho1 u1 p2 rho2 u2 p3 rho3 u3 p4 rho4 u4 % era USE INIT_PAR
global  vsh1 vsh2 vsh3 vsh4 vsh5 vsh3l vsh3r vsh4l vsh4r % era USE INIT_PAR
global  xsh1 xsh2 xsh3 xsh4 xsh5 xsh3l xsh3r xsh4l xsh4r % era USE INIT_PAR
global p t u s rho a e amach ptot ttot flow flht h htot % era USE VARS  <-- importanti
global w1 w2 w3 f1 f2 f3 phi1 phi2 phi3                    % era USE VARS
global pxeno uxeno hxeno ppxeno hhxeno ppt ut hht % era USE ENO
% LE VARIABILI PRIMITIVE POTREBBERO NON ESSERE p,u,h: ANDRANNO VALUTATE LE
% V.P. NEI PUNTI np E nm. 

enuo2 = 0.5*dt/dx; 

ncmm = ncm-1;
for n=2:ncmm
    nm=n;
    np=nm+1; % cerchiamo i flussi alle interfacce calcolando le grandezze nelle celle che ce l'hanno in comune
    
    pa  =  p(nm);  % LE VARIABILI PRIMITIVE POTREBBERO NON ESSERE p,u,h: ANDRANNO VALUTATE LE
    pb  =  p(np);  % V.P. NEI PUNTI np E nm.
    ua  =  u(nm);
    ub  =  u(np);
    ha  =  h(nm);
    hb  =  h(np);

    rhoa = pa/(ha/ga);
    rhob = pb/(hb/ga);

    aa = sqrt(gamma*pa/rhoa);
    ab = sqrt(gamma*pb/rhob);

    if n == 2 % A1 e A2
        p002   =  pa;
        u002   =  ua;
        h002   =  ha;
        rho002 =  p002/h002*ga;
        a002   =  sqrt(gamma*p002/rho002);
    end
    
    if n == 2
        a_starL002 = sqrt(2*(gamma+1)/(gamma-1)*h002);
        a_starR002 = sqrt(2*(gamma+1)/(gamma-1)*hb);
        a_tildeL002 = a_starL002^2/max(a_starL002,abs(u002));
        a_tildeR002 = a_starR002^2/max(a_starR002,abs(ub));
        a002 = min(a_tildeR002,a_tildeL002);
        M002 = u002/a002;
        Mb = ub/ab;
        a_n = a002;
        
        if abs(Ma0) >= 1 % valori di beta ottimale già settato a 1/8, alpha a 3/16
            M_cors_plus = 0.5*(Ma0+abs(Ma0));
            P_cors_plus = 0.5*(1+sign(Ma0));
        else
            M_cors_plus = 0.5*(Ma0+1)^2 + (1/8)*(M002^2-1)^2;
            P_cors_plus = 0.25*(M002+1)^2*(2-M002)+(3/16)*M002*(M002^2-1)^2;
        end

        if abs(Mb) >= 1  
            M_cors_minus = 0.5*(Mb-abs(Mb));
            P_cors_minus = 0.5*(1-sign(Mb));
        else
            M_cors_minus = 0.5*(Mb+1)^2 + (1/8)*(Mb^2-1)^2;
            P_cors_minus = 0.25*(Mb-1)^2*(2+Mb)-(3/16)*Mb*(Mb^2-1)^2;
        end

        m002 = M_cors_plus + M_cors_minus;
        p002 = P_cors_plus*p002 + P_cors_minus*pb;

        m_plus002 = 0.5*(m002 + abs(m002));
        m_minus002 = 0.5*(m002 - abs(m002));
        m_n_plus = m_plus002;
        m_n_minus = m_minus002;
        p_n = p002;
    end

    if n == ncmm
        pncm   =  pb;
        uncm   =  ub;
        hncm   =  hb;
        rhoncm =  pncm/hncm*ga;
        ancm   =  sqrt(gamma*pncm/rhoncm);
        Mncm   =  uncm/ancm;
        Ma     =  ua/aa;
    end

    if n == ncmm
        a_starLncm = sqrt(2*(gamma+1)/(gamma-1)*ha);
        a_starRncm = sqrt(2*(gamma+1)/(gamma-1)*hncm);
        a_tildeLncm = a_starLncm^2/max(a_starLncm,abs(ua));
        a_tildeRncm = a_starRncm^2/max(a_starRncm,abs(uncm));
        ancm = min(a_tildeRncm,a_tildeLncm);
        Mncm = uncm/ancm;
        a_n = ancm;
        
        if abs(Ma) >= 1 % valori di beta ottimale già settato a 1/8, alpha a 3/16
            M_cors_plus = 0.5*(Ma+abs(Ma));
            P_cors_plus = 0.5*(1+sign(Ma));
        else
            M_cors_plus = 0.5*(Ma+1)^2 + (1/8)*(Ma^2-1)^2;
            P_cors_plus = 0.25*(Ma+1)^2*(2-Ma)+(3/16)*Ma*(Ma^2-1)^2;
        end

        if abs(Mncm) >= 1  
            M_cors_minus = 0.5*(Mncm-abs(Mncm));
            P_cors_minus = 0.5*(1-sign(Mncm));
        else
            M_cors_minus = 0.5*(Mncm+1)^2 + (1/8)*(Mncm^2-1)^2;
            P_cors_minus = 0.25*(Mncm-1)^2*(2+Mncm)-(3/16)*Mncm*(Mncm^2-1)^2;
        end

        mncm = M_cors_plus + M_cors_minus;
        pncm = P_cors_plus*pa + P_cors_minus*pncm;

        m_plusncm = 0.5*(mncm + abs(mncm));
        m_minusncm = 0.5*(mncm - abs(mncm));
        m_n_plus = m_plusncm;
        m_n_minus = m_minusncm;
        p_n = pncm;
    end

% differenza significativa tra primo e secondo ordine: qui usiamo p,u,h
% come set di variabili. Al II ord, usiamo ln(p),u,ln(h): ciò consente che
% p e h non diventino mai negativi

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
    
    rhoa = pa/ha*ga;   % POTREMMO DOVER CALCOLARE ALTRE COSE  
    rhob = pb/hb*ga;
    % aa   = sqrt(gamma*pa/rhoa);
    % ab   = sqrt(gamma*pb/rhob);
    
    icalc=0;
    
    % r3a=pa+(rhoa*aa)*ua;  %--> r3a0r3c= circa a pc+(rhoa*aa)*uc 
    % r2a=ha-pa/rhoa;
    % r2b=hb-pb/rhob;
    % r1b=pb-(rhob*ab)*ub;
    % uc = (r3a-r1b)/((rhoa*aa)+(rhob*ab));
    % ud = uc;
    % pc = r3a-(rhoa*aa)*uc;  % approssimata, perché dovrebbe usare rhoc*ac
    % pd = pc;
    % hc = pc/rhoa+r2a; % approssimata, perché dovrebbe 
    % hd = pd/rhob+r2b;
    
    % if pc <= 0.0 || pd <= 0.0 || hc<= 0.0 || hd <= 0.0
    %     icalc=1;  % flag se qualche pressione diventa negativa
    % end
    
    % if icalc == 1  % in questo caso, switcha le variabili vecchie con quelle "logaritmiche" (tipiche del II ordine)
    %     ppa = log(pa);
    %     ppb = log(pb);
    %     hha = log(ha);
    %     hhb = log(hb);
    %     r3a = ppa + gamma/aa*ua;
    %     r2a = hha-ppa/ga;
    %     r2b = hhb-ppb/ga;
    %     r1b = ppb-gamma/ab*ub;
    %     uc = (r3a-r1b)/(gamma/aa+gamma/ab);
    %     ud = uc;
    %     ppc= r3a-gamma/aa*uc;
    %     ppd= ppc;
    %     hhc= ppc/ga+r2a;
    %     hhd= ppd/ga+r2b;
    %     pc = exp(ppc);
    %     pd = exp(ppd);
    %     hc = exp(hhc);
    %     hd = exp(hhd);
    %     fprintf('warning icalc=1 at k=%i and n=%i',k,n)
    % end
    
    % ac = sqrt(2.0*gd*hc);
    % ad = sqrt(2.0*gd*hd);
    % ala=ua-aa;
    % alc=uc-ac;
    % ald=ud+ad;
    % alb=ub+ab;
    % alx=uc;
    % [f1a,f2a,f3a] = decod(pa,ua,ha);  % mi preparo a calcolare i salti di flusso
    % [f1b,f2b,f3b] = decod(pb,ub,hb);
    % [f1c,f2c,f3c] = decod(pc,uc,hc);
    % [f1d,f2d,f3d] = decod(pd,ud,hd);  % A decod VANNO FORNITE LE VARIABILI  
    % GIUSTE, ED EVENTUALMENTE VANNO ESPLICITATI I FLUSSI IN FUNZIONE DELLA
    %  DATA TERNA DI VARIABILI
    % ...................................................

    % df1r = 0.0;
    % df2r = 0.0;
    % df3r = 0.0;
    % df1l = 0.0;
    % df2l = 0.0;
    % df3l = 0.0;
    % 
    % if (alx >=0) % c/d
    %     alam2r = 1.0;
    %     alam2l = 0.0;
    % else         % c\d
    %     alam2r = 0.0;
    %     alam2l = 1.0;
    % end
    % df1r =df1r + alam2r*(f1d-f1c);
    % df2r =df2r + alam2r*(f2d-f2c);
    % df3r =df3r + alam2r*(f3d-f3c);
    % df1l =df1l + alam2l*(f1d-f1c);
    % df2l =df2l + alam2l*(f2d-f2c);
    % df3l =df3l + alam2l*(f3d-f3c);
    % 
    % if (ala*alc >= 0)
    %     if (ala >=0) % a/c
    %         alam1r = 1.0;
    %         alam1l = 0.0;
    %     else         % a\c
    %         alam1r = 0.0;
    %         alam1l = 1.0;
    %     end
    %     df1r = df1r + alam1r*(f1c-f1a);
    %     df2r = df2r + alam1r*(f2c-f2a);
    %     df3r = df3r + alam1r*(f3c-f3a);
    %     df1l = df1l + alam1l*(f1c-f1a);
    %     df2l = df2l + alam1l*(f2c-f2a);
    %     df3l = df3l + alam1l*(f3c-f3a);
    % else % transonic case
    %     gdgd=2.0*gd;
    %     if (icalc == 0)
    %         sqrhst=0.5*(-sqrt(gdgd)*aa+sqrt(gdgd*aa*aa+4.0*(r3a/rhoa+r2a)));
    %         hst = sqrhst^2;
    %         ast = sqrt(gdgd)*sqrhst;
    %         ust = ast;
    %         pst = rhoa*(-r2a+hst);
    %     end
    %     if (icalc == 1)
    %         adum = 2.0*aa/gdgd^1.5;
    %         bdum = aa/gamma/sqrt(gdgd)*(r3a+ga*r2a);
    %         iterstar (adum,bdum,cdum)
    %         ast = sqrt(gdgd)*cdum;
    %         ust = ast;
    %         hst = gb*ast^2;
    %         pst= exp((log(hst)-r2a)*ga);
    %     end
    %     [f1st,f2st,f3st] = decod(pst,ust,hst);
    %     if (ala <=0) % a\*/c
    %         df1r = df1r + (f1c-f1st);
    %         df2r = df2r + (f2c-f2st);
    %         df3r = df3r + (f3c-f3st);
    %         df1l = df1l + (f1st-f1a);
    %         df2l = df2l + (f2st-f2a);
    %         df3l = df3l + (f3st-f3a);
    %     else          % a/*\c
    %         df1l = df1l + (f1c-f1st);
    %         df2l = df2l + (f2c-f2st);
    %         df3l = df3l + (f3c-f3st);
    %         df1r = df1r + (f1st-f1a);
    %         df2r = df2r + (f2st-f2a);
    %         df3r = df3r + (f3st-f3a);
    %     end
    % end
    % 
    % if (ald*alb >= 0)
    %     if (ald >=0) %d/b
    %         alam3r = 1.0;
    %         alam3l = 0.0;
    %     else         %d\b
    %         alam3r = 0.0;
    %         alam3l = 1.0;
    %     end
    %     df1r = df1r + alam3r*(f1b-f1d);
    %     df2r = df2r + alam3r*(f2b-f2d);
    %     df3r = df3r + alam3r*(f3b-f3d);
    %     df1l = df1l + alam3l*(f1b-f1d);
    %     df2l = df2l + alam3l*(f2b-f2d);
    %     df3l = df3l + alam3l*(f3b-f3d);
    % else
    %     gdgd=2.0*gd;
    %     if (icalc == 0)
    %         sqrhst=0.5*(-sqrt(gdgd)*ab+sqrt(gdgd*ab*ab+4.*(r1b/rhob+r2b)));
    %         hst = sqrhst^2;
    %         ast = sqrt(gdgd)*sqrhst;
    %         ust = -ast;
    %         pst = rhob*(-r2b+hst);
    %     end
    %     if (icalc == 1)
    %         adum = 2.0 * ab / gdgd^1.5;
    %         bdum = ab / gamma / sqrt(gdgd) * (r1b + ga * r2b);
    %         iterstar(adum, bdum, cdum);
    %         ast = sqrt(gdgd) * cdum;
    %         ust = -ast;
    %         hst = gb * ast^2;
    %         pst = exp((log(hst) - r2b) * ga);
    %     end
    %     [f1st,f2st,f3st] = decod(pst,ust,hst);
    %     if (ald <=0) %d\*/b
    %         df1r = df1r + (f1b-f1st);
    %         df2r = df2r + (f2b-f2st);
    %         df3r = df3r + (f3b-f3st);
    %         df1l = df1l + (f1st-f1d);
    %         df2l = df2l + (f2st-f2d);
    %         df3l = df3l + (f3st-f3d);
    %     else         %d/*\b
    %         df1l = df1l + (f1b-f1st);
    %         df2l = df2l + (f2b-f2st);
    %         df3l = df3l + (f3b-f3st);
    %         df1r = df1r + (f1st-f1d);
    %         df2r = df2r + (f2st-f2d);
    %         df3r = df3r + (f3st-f3d);
    %     end
    % end

    rhoa = pa/ha*ga;  
    rhob = pb/hb*ga;
    aa  = sqrt(gamma*pa/rhoa);
    ab  = sqrt(gamma*pb/rhob);
    Ma  =  ua/aa;
    Mb  =  ub/ab;
    
    a_star_a = sqrt(2*(gamma-1)/(gamma+1)*ha);
    a_star_b = sqrt(2*(gamma-1)/(gamma+1)*b);
    a_tilde_a = a_star_a^2/max(a_star_a,abs(ua));
    a_tilde_b = a_star_b^2/max(a_star_b,abs(ub));
    a_n = min(a_tilde_a,a_tilde_b);
    Ma = ua/a_n;
    Mb = ub/a_n;

    if n > 2 && n < ncmm
        if abs(Ma) >= 1 % valori di beta ottimale già settato a 1/8, alpha a 3/16
            M_cors_plus = 0.5*(Ma+abs(Ma));
            P_cors_plus = 0.5*(1+sign(Ma));
        else
            M_cors_plus = 0.5*(Ma+1)^2 + (1/8)*(Ma^2-1)^2;
            P_cors_plus = 0.25*(Ma+1)^2*(2-Ma)+(3/16)*Ma*(Ma^2-1)^2;
        end
    
        if abs(Mb) >= 1  
            M_cors_minus = 0.5*(Mb-abs(Mb));
            P_cors_minus = 0.5*(1-sign(Mb));
        else
            M_cors_minus = 0.5*(Mb+1)^2 + (1/8)*(Mb^2-1)^2;
            P_cors_minus = 0.25*(Mb-1)^2*(2+Mb)-(3/16)*Mb*(Mb^2-1)^2;
        end
    
        m_n = M_cors_plus + M_cors_minus;
        p_n = P_cors_plus*pa + P_cors_minus*pb;
        m_n_plus = 0.5*(m_n + abs(m_n));
        m_n_minus = 0.5*(m_n - abs(m_n));
    end

    phi1(n) = a_n*(m_n_plus*rhoa+m_n_minus*rhob); %A3
    phi2(n) = a_n*(m_n_plus*rhoa*ua+m_n_minus*rhob*ub)+p_n;
    phi3(n) = a_n*(m_n_plus*rhoa*ha+m_n_minus*rhob*hb);
end % end of the do loop

if(itest == 1)  %REFLECTING WALL B.C.
    r1dum  = p002-rho002*a002*u002;
    r2dum  = h002-p002/rho002;
    pin    = r1dum;
    uin    = 0.0;
    hin    = r2dum+pin/rho002;
    [phi1(1),phi2(1),phi3(1)] = decod(pin,uin,hin); 
    r3dum  = pncm+rhoncm*ancm*uncm;
    r2dum  = hncm-pncm/rhoncm;
    pex    = r3dum;
    uex    = 0.0;
    hex    = r2dum+pex/rhoncm;
    [phi1(ncm),phi2(ncm),phi3(ncm)] = decod(pex,uex,hex);
end
if(itest >= 2)  % REFLECTING WALL B.C.
    r1dum  = p002-rho002*a002*u002;
    r2dum  = h002-p002/rho002;
    ain    = a002;
    pin    = p002;
    uin    = u002;
    hin    = h002;
    rhoin  = rho002;
    rho_2 = p(3)/h(3)*ga;
    m_minus_in = m_minus002;
    m_plus_in = m_plus002;
    phi1(1) = ain*(m_plus_in*rhoin+m_minus_in*rho_2);
    phi2(1) = ain*(m_plus_in*rhoin*uin+m_minus_in*rho_2*u(3))+pin;
    phi3(1) = ain*(m_plus_in*rhoin*hin+m_minus_in*rho_2*h(3));
    r3dum  = pncm+rhoncm*ancm*uncm;
    r2dum  = hncm-pncm/rhoncm;
    aex    = ancm;
    pex    = pncm;
    uex    = uncm;
    hex    = hncm;
    rhoex = rhoncm;
    rho_ncmm = p(ncmm)/h(ncmm)*ga;
    m_minus_ex = m_minusncm;
    m_plus_ex = m_plusncm;
    phi1(ncm) = aex*(m_plus_ex*rho_ncmm+m_minus_ex*rhoex);
    phi2(ncm) = aex*(m_plus_ex*rho_ncmm*u(ncmm)+m_minus_ex*rhoex*uex)+pex;
    phi3(ncm) = aex*(m_plus_ex*rho_ncmm*h(ncmm)+m_minus_ex*rhoex*hex);
end
end