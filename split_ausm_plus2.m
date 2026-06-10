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
global ischeme

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
    % da riga 250 in poi
    %% Condizioni al bordo
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

        % In questo caso ignoriamo il calcolo di a_tilde poiché 
        % a002 = a_tilde002
        %a_tildeL002 = a_starL002^2/max(a_starL002,abs(u002)); %u002 = ua
        %a_tildeR002 = a_starR002^2/max(a_starR002,abs(ub));

        %a002 = min(a_tildeR002,a_tildeL002);

        a_star002 = sqrt(2*(gamma-1)/(gamma+1)*h_t002);
        a002 = a_star002^2/max(a_star002,abs(u002)); 

        M002 = u002/a002;

        M_cors_plus = M_beta(M002, 1/8, 'plus');
        M_cors_minus = M_beta(M002, 1/8, 'minus');
        
        P_cors_plus = P_alpha(M002, 3/16, 'plus');
        P_cors_minus = P_alpha(M002, 3/16, 'minus');

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

        a_starncm = sqrt(2*(gamma-1)/(gamma+1)*h_tncm);
        ancm = a_starncm^2/max(a_starncm,abs(uncm));

        Mncm = uncm/ancm;

        M_cors_plus = M_beta(Mncm, 1/8, 'plus');
        M_cors_minus = M_beta(Mncm, 1/8, 'minus');

        P_cors_plus = P_alpha(Mncm, 3/16, 'plus');
        P_cors_minus = P_alpha(Mncm, 3/16, 'minus');

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

        if itest == 1 || itest == 3
            a_n = 0.5*(a_star_a+a_star_b);
        else
            a_n = min(a_tilde_a,a_tilde_b);
        end

        Ma = ua/a_n;
        Mb = ub/a_n;
        
        M_cors_plus = M_beta(Ma, 1/8, 'plus');
        M_cors_minus = M_beta(Mb, 1/8, 'minus');

        P_cors_plus = P_alpha(Ma, 3/16, 'plus');
        P_cors_minus = P_alpha(Mb, 3/16, 'minus');
        

        m_n = M_cors_plus + M_cors_minus;
        p_n = P_cors_plus*pa + P_cors_minus*pb;
        m_n_plus = 0.5*(m_n+abs(m_n));
        m_n_minus = 0.5*(m_n-abs(m_n));
    
        if ischeme == 2
            %Caso AUSM+
            % Formula A3 del paper scomposta nelle 3 componenti:
            phi1(n) = a_n*(m_n_plus*rhoa + m_n_minus*rhob); 
            phi2(n) = a_n*(m_n_plus*rhoa*ua + m_n_minus*rhob*ub) + p_n;
            phi3(n) = a_n*(m_n_plus*rhoa*h_ta + m_n_minus*rhob*h_tb);
        elseif ischeme == 3
            %AUSMPW
            % P_cors_plus e P_cors_minus sono le P maiuscole dell'AUSMPW, la
            % formula a pagina 318 per p_s ha probabilmente un errore nell'ultima
            % P_R
            p_s = P_cors_plus * pa + P_cors_minus * pb;
            fa = f_limiter(Ma, ua, a_n, pa, pb, p_s, 'left');
            fb = f_limiter(Mb, ub, a_n, pa, pb, p_s, 'right');
            
            %Controlliamo il valore di m minuscolo e decidiamo come
            %procedere
            if m_n >= 0
                %Usiamo la prima formula delle (42) di pagina 318 del paper
                phi1(n) = a_n*((1+fa)*M_cors_plus* rhoa + ...
                    (1+fb)*M_cors_minus* pw(rhoa,rhob) );
                phi2(n) = a_n*((1+fa)*M_cors_plus* rhoa*ua + ...
                    (1+fb)*M_cors_minus* pw(rhoa*ua,rhob*ub) ) + ...
                    (P_cors_plus*pa + P_cors_minus*pb);
                phi3(n) = a_n*((1+fa)*M_cors_plus* rhoa*h_ta+ ...
                    (1+fb)*M_cors_minus* pw(rhoa,rhob)*h_ta );
            else
                %Usiamo la seconda formula delle (42) di pagina 318 del paper
                phi1(n) = a_n*((1+fa)*M_cors_plus* pw(rhob,rhoa) + ...
                    (1+fb)*M_cors_minus* rhob );
                phi2(n) = a_n*((1+fa)*M_cors_plus* pw(rhob*ub,rhoa*ua) + ...
                    (1+fb)*M_cors_minus* rhob*ub ) + ...
                    (P_cors_plus*pa + P_cors_minus*pb);
                phi3(n) = a_n*((1+fa)*M_cors_plus* pw(rhob,rhoa)*h_tb + ...
                    (1+fb)*M_cors_minus* rhob*h_tb );
            end
        else
            fprintf(['What just happened? How did you find a scheme ' ...
                'that we have not implemented?'])
        end
    end

end %end of the for loop

%% Calcolo delle variabili sul bordo
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

if ischeme == 2
    phi1(1) = ain*(m_plus_in*rhoin + m_minus_in*rho_2);
    phi2(1) = ain*(m_plus_in*rhoin*uin + m_minus_in*rho_2*u(3)) + pin;
    phi3(1) = ain*(m_plus_in*rhoin*h_tin + m_minus_in*rho_2*h_t_2);
elseif ischeme == 3
    %Ci serve nuovamente definire i valori di tutte le variabili, della
    %cella a sinistra (indicata con in) e di quella a destra (indicata
    %con _2)
    p_2 = p(3);
    u_2 = u(3);
    a_n; %Da definire
    M_2; %Da definire
    p_s = P_cors_plus002 * pin + P_cors_minus002 * p_2;
    f002 = f_limiter(M002, uin, a_n, pin, p_2, p_s, 'left');
    fb = f_limiter(M_2, u_2, a_n, pin, p_2, p_s, 'right');
    
    if m_minus_in >= 0
        phi1(n) = a_n*((1+f002)*M_cors_plus* rhoa + ...
            (1+fb)*M_cors_minus* pw(rhoa,rhob) );
        phi2(n) = a_n*((1+f002)*M_cors_plus* rhoa*ua + ...
            (1+fb)*M_cors_minus* pw(rhoa*ua,rhob*ub) ) + ...
            (P_cors_plus*pa + P_cors_minus*pb);
        phi3(n) = a_n*((1+f002)*M_cors_plus* rhoa*h_ta+ ...
            (1+fb)*M_cors_minus* pw(rhoa,rhob)*h_ta );
    else
        phi1(n) = a_n*((1+f002)*M_cors_plus* pw(rhob,rhoa) + ...
            (1+fb)*M_cors_minus* rhob );
        phi2(n) = a_n*((1+f002)*M_cors_plus* pw(rhob*ub,rhoa*ua) + ...
            (1+fb)*M_cors_minus* rhob*ub ) + ...
            (P_cors_plus*pa + P_cors_minus*pb);
        phi3(n) = a_n*((1+f002)*M_cors_plus* pw(rhob,rhoa)*h_tb + ...
            (1+fb)*M_cors_minus* rhob*h_tb );
    end
end

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

end %end of function

%% Implementation of some functions used for AUSM+ and AUSMPW

%Maybe H is the total enthalpy e not the h_t specific used above?

%pl computes the value of pl(x, y) - from formula (40) of page 318
%of the AUSMPW paper
function result = pl(x, y)
    m = min(x/y, y/x);
    if m >= 3/4 && m < 1
        result = 4 * m - 3;
        return
    elseif m < 3/4 && m >= 0
        result = 0.0;
        return
    else
        fprintf(['Something went wrong with the computation of pl; ' ...
            'x and y are probably the same number\n'])
        result = -1;
    end
end

%w computes the value of w(x, y) - from page 319
function result = w(x, y) 
    result = 1 - (min(x/y, y/x))^3;
end

%pw computes the value of pw(x, y) - from page 319
function result = pw(x,y) 
    result = (1 - w(x, y)) * x + w(x,y) * y;
end

%M_beta calculates the value of M_beta depending on which beta you
%give it as an input and which sign you wish your M_beta to have
function result = M_beta(M, beta, sgn)
    %First, we need to decide the sign of M_beta
    if strcmpi(sgn, 'plus') %strcmpi gives true (or 1) if the two strings 
        % are the same (it is case insensitive), otherwise it gives 
        % false (or 0) as a result

        if abs(M) <= 1
            result = 0.25*(M+1)^2 + beta*(M^2-1)^2;
            return
        else
            result = 0.5*(M+abs(M));
            return
        end

    elseif strcmpi(sgn, 'minus')

        if abs(M) <= 1
            result = -0.25*(M-1)^2 - beta*(M^2-1)^2;
            return
        else
            result = 0.5*(M-abs(M));
            return
        end
        
    else
        fprintf(['Something went wrong: have you spelt "plus" and "minus" ' ...
            'correctly in every usage of M_beta? \n'] )
    end
end

%f_limiter computes the value of the function f "similar to a limiter"
%introduced with the AUSMPW algorithm - formula (39) of page 318
%As inputs, it needs:
% - the Mach number and the velocity of the correct cell (left or right);
% - the value of the speed of sound on the interface;
% - the value of pressure on the two consecutives cell;
% - the values of p_s;
% - whether you want to calculate the 'left' or the 'right' one.
function result = f_limiter(M, u, a_n, p_l, p_r, p_s, position)
    if abs(M) < 1
        if strcmpi(position, 'left')
            result = (p_l/p_s - 1)*pl(p_l,p_r)*abs(M_beta(M,0,'plus')* ...
                min(1, (abs(u)/a_n)^(0.25)));
            return
        elseif strcmpi(position, 'right')
            result = (p_r/p_s - 1)*pl(p_r,p_l)*abs(M_beta(M,0,'minus')* ...
                min(1, (abs(u)/a_n)^(0.25)));
            return
        else
            fprintf(['Something went wrong: have you spelt "left" and "right" ' ...
            'correctly in every usage of f_limiter? \n'] )
        end
    else
        result = 0.0;
    end
end

function result = P_alpha(M, alpha, sgn)
    %First, we need to decide the sign of P_alpha
    if strcmpi(sgn, 'plus')

        if abs(M) <= 1
            result = 0.25*(M+1)^2*(2-M) + alpha*M*(M^2-1)^2;
            return
        else
            result = 0.5*(1+sign(M));
            return
        end

    elseif strcmpi(sgn, 'minus')

        if abs(M) <= 1
            result = 0.25*(M-1)^2*(2+M) - alpha*M*(M^2-1)^2;
            return
        else
            result = 0.5*(1-sign(M));
            return
        end
        
    else
        fprintf(['Something went wrong: have you spelt "plus" and "minus" ' ...
            'correctly in every usage of P_alpha? \n'] )
    end
end

%}