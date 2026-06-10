function split

%Input e output di questa fuznione sono tutte queste variabili che passiamo
%tra le funzioni che hanno queste righe come introduzione. 
% In particolare, a noi interessano p,u e h
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

%Questa è la funzione da variare (se ce n'è bisogno) per il nostro progetto

enuo2 = 0.5*dt/dx;

ncmm = ncm-1;
for n=2:ncmm

    % Prendiamo in considerazione la cella n-esima e quella successiva
    nm=n;
    np=nm+1; % np=n+1
    
    %Potremmo non dover calcolare p, u e h
    pa  =  p(nm);
    pb  =  p(np);
    ua  =  u(nm);
    ub  =  u(np);
    ha  =  h(nm);
    hb  =  h(np);
    
    % Ai contorni condizioni speciali
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
    
    if iord ~= 1  % si usa quando l'ordine è 2

        %Si fa un cambio di variabile passando ai logaritmi per evitare che
        %le variabili diventino negative
        ppa = log(pa);
        ppb = log(pb);
        hha = log(ha);
        hhb = log(hb);

        %Calcolo le grandezze alle interfacce delle celle
        ppa  =  ppa + .5* ppxeno(n);
        ppb  =  ppb - .5* ppxeno(n+1);
        ua   =  ua  + .5* uxeno(n);
        ub   =  ub  - .5* uxeno(n+1);
        hha  =  hha + .5* hhxeno(n);
        hhb  =  hhb - .5* hhxeno(n+1);

        %Aggiungo anche le derivate temporali
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
    
    %Potremmo non dover calcolare rho e a
    rhoa = pa/ha*ga;
    rhob = pb/hb*ga;
    aa   = sqrt(gamma*pa/rhoa);
    ab   = sqrt(gamma*pb/rhob);
    
    % Variabile che dà errore nel caso in cui qualcosa diventi negativo
    icalc=0;
    
    r3a=pa+(rhoa*aa)*ua;
    r2a=ha-pa/rhoa;
    r2b=hb-pb/rhob;
    r1b=pb-(rhob*ab)*ub;
    uc = (r3a-r1b)/((rhoa*aa)+(rhob*ab));
    ud = uc;
    pc = r3a-(rhoa*aa)*uc;
    pd = pc;
    hc = pc/rhoa+r2a;
    hd = pd/rhob+r2b;
    
    % Con espansioni molto forti, le quantità possono diventare negative e
    % ciò farebbe saltare l'algoritmo.
    if pc <= 0.0 || pd <= 0.0 || hc<= 0.0 || hd <= 0.0
        icalc=1;
    end
    
    % Per evitare che l'algoritmo salti, calcoliamo tutto attraverso i
    % logaritmi.
    if icalc == 1
        ppa = log(pa);
        ppb = log(pb);
        hha = log(ha);
        hhb = log(hb);
        r3a = ppa + gamma/aa*ua;
        r2a = hha-ppa/ga;
        r2b = hhb-ppb/ga;
        r1b = ppb-gamma/ab*ub;
        uc = (r3a-r1b)/(gamma/aa+gamma/ab);
        ud = uc;
        ppc= r3a-gamma/aa*uc;
        ppd= ppc;
        hhc= ppc/ga+r2a;
        hhd= ppd/ga+r2b;
        pc = exp(ppc);
        pd = exp(ppd);
        hc = exp(hhc);
        hd = exp(hhd);
        fprintf('warning icalc=1 at k=%i and n=%i \n',k,n)
    end
    
    ac = sqrt(2.0*gd*hc);
    ad = sqrt(2.0*gd*hd);
    ala=ua-aa; %lambda_1 nella zona a
    alc=uc-ac;
    ald=ud+ad;
    alb=ub+ab;
    alx=uc;
    % Si calcolano i flussi delle varie quantità nelle diverse zone,
    % passando da p,u e h alla terna che ci serve ( rho*u, p+rho*u^2 e
    % (E+p)*u )
    [f1a,f2a,f3a] = decod(pa,ua,ha);
    [f1b,f2b,f3b] = decod(pb,ub,hb);
    [f1c,f2c,f3c] = decod(pc,uc,hc);
    [f1d,f2d,f3d] = decod(pd,ud,hd);
    % ...................................................

    df1r = 0.0;
    df2r = 0.0;
    df3r = 0.0;
    df1l = 0.0;
    df2l = 0.0;
    df3l = 0.0;
    
    if (alx >=0) % c/d -> Inclinazione superficie di contatto tra zona c 
        % e zona d
        alam2r = 1.0;
        alam2l = 0.0;
    else % c\d
        alam2r = 0.0;
        alam2l = 1.0;
    end
    df1r =df1r + alam2r*(f1d-f1c);
    df2r =df2r + alam2r*(f2d-f2c);
    df3r =df3r + alam2r*(f3d-f3c);
    df1l =df1l + alam2l*(f1d-f1c);
    df2l =df2l + alam2l*(f2d-f2c);
    df3l =df3l + alam2l*(f3d-f3c);
    
    if (ala*alc >= 0)
        if (ala >=0) % a/c
            alam1r = 1.0;
            alam1l = 0.0;
        else         % a\c
            alam1r = 0.0;
            alam1l = 1.0;
        end
        df1r = df1r + alam1r*(f1c-f1a);
        df2r = df2r + alam1r*(f2c-f2a);
        df3r = df3r + alam1r*(f3c-f3a);
        df1l = df1l + alam1l*(f1c-f1a);
        df2l = df2l + alam1l*(f2c-f2a);
        df3l = df3l + alam1l*(f3c-f3a);
    else % transonic case
        gdgd=2.0*gd;
        if (icalc == 0)
            sqrhst=0.5*(-sqrt(gdgd)*aa+sqrt(gdgd*aa*aa+4.0*(r3a/rhoa+r2a)));
            hst = sqrhst^2;
            ast = sqrt(gdgd)*sqrhst;
            ust = ast;
            pst = rhoa*(-r2a+hst);
        end
        if (icalc == 1)
            adum = 2.0*aa/gdgd^1.5;
            bdum = aa/gamma/sqrt(gdgd)*(r3a+ga*r2a);
            cdum = iterstar(adum,bdum);
            ast = sqrt(gdgd)*cdum;
            ust = ast;
            hst = gb*ast^2;
            pst= exp((log(hst)-r2a)*ga);
        end
        [f1st,f2st,f3st] = decod(pst,ust,hst);

        if (ala <=0) % a\*/c
            df1r = df1r + (f1c-f1st);
            df2r = df2r + (f2c-f2st);
            df3r = df3r + (f3c-f3st);
            df1l = df1l + (f1st-f1a);
            df2l = df2l + (f2st-f2a);
            df3l = df3l + (f3st-f3a);
        else          % a/*\c
            df1l = df1l + (f1c-f1st);
            df2l = df2l + (f2c-f2st);
            df3l = df3l + (f3c-f3st);
            df1r = df1r + (f1st-f1a);
            df2r = df2r + (f2st-f2a);
            df3r = df3r + (f3st-f3a);
        end
    end
    
    if (ald*alb >= 0)
        if (ald >=0) %d/b
            alam3r = 1.0;
            alam3l = 0.0;
        else         %d\b
            alam3r = 0.0;
            alam3l = 1.0;
        end
        df1r = df1r + alam3r*(f1b-f1d);
        df2r = df2r + alam3r*(f2b-f2d);
        df3r = df3r + alam3r*(f3b-f3d);
        df1l = df1l + alam3l*(f1b-f1d);
        df2l = df2l + alam3l*(f2b-f2d);
        df3l = df3l + alam3l*(f3b-f3d);
    else
        gdgd=2.0*gd;
        if (icalc == 0)
            sqrhst=0.5*(-sqrt(gdgd)*ab+sqrt(gdgd*ab*ab+4.*(r1b/rhob+r2b)));
            hst = sqrhst^2;
            ast = sqrt(gdgd)*sqrhst;
            ust = -ast;
            pst = rhob*(-r2b+hst);
        end
        if (icalc == 1)
            adum = 2.0 * ab / gdgd^1.5;
            bdum = ab / gamma / sqrt(gdgd) * (r1b + ga * r2b);
            cdum = iterstar(adum, bdum);
            ast = sqrt(gdgd) * cdum;
            ust = -ast;
            hst = gb * ast^2;
            pst = exp((log(hst) - r2b) * ga);
        end
        [f1st,f2st,f3st] = decod(pst,ust,hst);
        if (ald <=0) %d\*/b
            df1r = df1r + (f1b-f1st);
            df2r = df2r + (f2b-f2st);
            df3r = df3r + (f3b-f3st);
            df1l = df1l + (f1st-f1d);
            df2l = df2l + (f2st-f2d);
            df3l = df3l + (f3st-f3d);
        else         %d/*\b
            df1l = df1l + (f1b-f1st);
            df2l = df2l + (f2b-f2st);
            df3l = df3l + (f3b-f3st);
            df1r = df1r + (f1st-f1d);
            df2r = df2r + (f2st-f2d);
            df3r = df3r + (f3st-f3d);
        end
    end

    %Queste sono la variabili che poi si usano per aggiornare le variabili
    %nel file march.m:
    phi1(n)=f1a+df1l;
    phi2(n)=f2a+df2l;
    phi3(n)=f3a+df3l;

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
    pin    = p002;
    uin    = u002;
    hin    = h002;
    [phi1(1),phi2(1),phi3(1)] = decod(pin,uin,hin);
    r3dum  = pncm+rhoncm*ancm*uncm;
    r2dum  = hncm-pncm/rhoncm;
    pex    = pncm;
    uex    = uncm;
    hex    = hncm;
    [phi1(ncm),phi2(ncm),phi3(ncm)] = decod(pex,uex,hex);
end
end