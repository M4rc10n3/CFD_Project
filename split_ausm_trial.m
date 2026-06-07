function split_ausm_trial

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

enuo2 = 0.5*dt/dx;
ncmm = ncm-1;

% ============================
% LOOP INTERNO
% ============================
for n = 2:ncmm

    nm = n;
    np = n+1;

    % -------------------------------------------------
    % 1. Stati base
    % -------------------------------------------------
    pa = p(nm);
    pb = p(np);
    ua = u(nm);
    ub = u(np);
    ha = h(nm);
    hb = h(np);

    % -------------------------------------------------
    % 2. ENO / MUSCL (se attivo)
    % -------------------------------------------------
    if iord ~= 1

        ppa = log(pa);
        ppb = log(pb);
        hha = log(ha);
        hhb = log(hb);

        ppa = ppa + 0.5*ppxeno(n);
        ppb = ppb - 0.5*ppxeno(n+1);
        ua  = ua  + 0.5*uxeno(n);
        ub  = ub  - 0.5*uxeno(n+1);
        hha = hha + 0.5*hhxeno(n);
        hhb = hhb - 0.5*hhxeno(n+1);

        ppa = ppa + enuo2*ppt(nm);
        ua  = ua  + enuo2*ut(nm);
        hha = hha + enuo2*hht(nm);

        ppb = ppb + enuo2*ppt(np);
        ub  = ub  + enuo2*ut(np);
        hhb = hhb + enuo2*hht(np);

        pa = exp(ppa);
        pb = exp(ppb);
        ha = exp(hha);
        hb = exp(hhb);
    end

    % -------------------------------------------------
    % 3. Variabili termodinamiche coerenti
    % -------------------------------------------------
    rhoa = pa/(ha/ga);
    rhob = pb/(hb/ga);

    aa = sqrt(gamma*pa/rhoa);
    ab = sqrt(gamma*pb/rhob);

    Ma = ua/aa;
    Mb = ub/ab;

    % -------------------------------------------------
    % 4. Velocità del suono all'interfaccia (AUSM+ robusto)
    % -------------------------------------------------
    a_star_a = sqrt(((gamma+1)/2) * (pa/rhoa + 0.5*ua^2));
    a_star_b = sqrt(((gamma+1)/2) * (pb/rhob + 0.5*ub^2));

    a_tilde_a = a_star_a^2 / max(a_star_a, abs(ua));
    a_tilde_b = a_star_b^2 / max(a_star_b, abs(ub));

    a_n = min(a_tilde_a, a_tilde_b);

    % safety
    if a_n <= 0 || ~isfinite(a_n)
        error('a_n non valido al nodo %d', n);
    end

    % -------------------------------------------------
    % 5. AUSM+ Mach splitting (Liou)
    % -------------------------------------------------
    if abs(Ma) >= 1
        Mplus_a = 0.5*(Ma + abs(Ma));
        Pplus_a = 0.5*(1 + sign(Ma));
    else
        Mplus_a = 0.5*(Ma+1)^2 + (1/8)*(Ma^2-1)^2;
        Pplus_a = 0.25*(Ma+1)^2*(2-Ma) + (3/16)*Ma*(Ma^2-1)^2;
    end

    if abs(Mb) >= 1
        Mminus_b = 0.5*(Mb - abs(Mb));
        Pminus_b = 0.5*(1 - sign(Mb));
    else
        Mminus_b = -0.5*(Mb-1)^2 + (1/8)*(Mb^2-1)^2;
        Pminus_b = 0.25*(Mb-1)^2*(2+Mb) - (3/16)*Mb*(Mb^2-1)^2;
    end

    % -------------------------------------------------
    % 6. Flux interface
    % -------------------------------------------------
    m_half = Mplus_a + Mminus_b;

    p_half = Pplus_a*pa + Pminus_b*pb;

    m_plus  = 0.5*(m_half + abs(m_half));
    m_minus = 0.5*(m_half - abs(m_half));

    % -------------------------------------------------
    % 7. Flux AUSM+
    % -------------------------------------------------
    phi1(n) = a_n * (m_plus*rhoa + m_minus*rhob);

    phi2(n) = a_n * (m_plus*rhoa*ua + m_minus*rhob*ub) + p_half;

    phi3(n) = a_n * (m_plus*rhoa*ha + m_minus*rhob*hb);

end

% ============================
% BORDO SINISTRO
% ============================
if itest >= 1

    % ghost cell sinistra (riflessione tipica)
    pa = p(2);
    pb = p(2);
    ua = -u(2);
    ub = u(2);
    ha = h(2);
    hb = h(2);

    rhoa = pa/(ha/ga);
    rhob = pb/(hb/ga);

    aa = sqrt(gamma*pa/rhoa);
    ab = sqrt(gamma*pb/rhob);

    Ma = ua/aa;
    Mb = ub/ab;

    a_n = aa;

    m_half = 0.5*(Ma + abs(Ma)) + 0.5*(Mb - abs(Mb));

    p_half = pa;

    m_plus  = 0.5*(m_half + abs(m_half));
    m_minus = 0.5*(m_half - abs(m_half));

    phi1(1) = a_n*(m_plus*rhoa + m_minus*rhob);
    phi2(1) = a_n*(m_plus*rhoa*ua + m_minus*rhob*ub) + p_half;
    phi3(1) = a_n*(m_plus*rhoa*ha + m_minus*rhob*hb);

end

% ============================
% BORDO DESTRO
% ============================
if itest >= 1

    pa = p(ncmm);
    pb = p(ncmm);
    ua = u(ncmm);
    ub = -u(ncmm);
    ha = h(ncmm);
    hb = h(ncmm);

    rhoa = pa/(ha/ga);
    rhob = pb/(hb/ga);

    aa = sqrt(gamma*pa/rhoa);
    ab = sqrt(gamma*pb/rhob);

    Ma = ua/aa;
    Mb = ub/ab;

    a_n = aa;

    m_half = 0.5*(Ma + abs(Ma)) + 0.5*(Mb - abs(Mb));

    p_half = pa;

    m_plus  = 0.5*(m_half + abs(m_half));
    m_minus = 0.5*(m_half - abs(m_half));

    phi1(ncm) = a_n*(m_plus*rhoa + m_minus*rhob);
    phi2(ncm) = a_n*(m_plus*rhoa*ua + m_minus*rhob*ub) + p_half;
    phi3(ncm) = a_n*(m_plus*rhoa*ha + m_minus*rhob*hb);

end

end