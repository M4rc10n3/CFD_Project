function init(rhol,ul,pl,rhor,ur,pr)

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIGAMMA_PAR
global nc ncm k kout ka iord itest stab % era USE NUM_PAR
global b c diverg c1 c2 x0 dx                     % era USE GEOM_PAR
global   time dt dtodx dtitest4 timemax timek     % era USE TIME_PAR
global  xsh10 xsh20 p1 rho1 u1 p2 rho2 u2 p3 rho3 u3 p4 rho4 u4 % era USE INIT_PAR
global  vsh1 vsh2 vsh3 vsh4 vsh5 vsh3l vsh3r vsh4l vsh4r % era USE INIT_PAR
global  xsh1 xsh2 xsh3 xsh4 xsh5 xsh3l xsh3r xsh4l xsh4r % era USE INIT_PAR
global  xcsf2                                            % era USE INIT_PAR
global x                                          % era USE GEOM_VAR
global p t u s rho a e amach ptot ttot flow flht h htot % era USE VARS
global w1 w2 w3 f1 f2 f3 phi1 phi2 phi3                    % era USE VARS
global pxeno uxeno hxeno ppxeno hhxeno ppt ut hht  % era USE ENO

%     *** Initialize deg-rad conversion and gamma combinations ***

% Gamma è un dato di input visibile nella finestra di dialogo
ga=gamma/(gamma-1.);
gb=1.d0/(gamma-1.);
gc=(gamma+1.)/(gamma-1.);
gd=(gamma-1.)/2.;
ge=(gamma+1.)/2.;
gf = sqrt(gamma);
gg=2./(gamma-1.);
gh=(gamma+1.)/(2.*gamma);
gi=(gamma-1.)/(2.*gamma);
gj=(gamma-1.)/gamma;

%                     *** Test case selection ***

if itest == 1 % shock tube
    fprintf('euler 1d - shock tube test case')
    dtodx = 0.411;  %default value for 100 points and 35 t-steps
    ka    = 0.35*nc+1;      %default value
   
    kout=ka;
    
    % A partire dagli input calcoliamo le prime variabili in modo
    % adimensionale
    timemax = 1.e6;
    ratrho = rhol/rhor; %Semplice rapporto delle densità
    ratp   = pl/pr; %Semplice rapporto delle pressioni
    pmin   =1.; %Pressione di riferimento
    rhomin =1.; %Densità di riferimento
    tmin   =pmin/rhomin; %Temperatura di riferimento
    umin   = 0.0; % La corrente è ferma, dunque la velocità è 0
    pmax   =pmin*ratp;
    rhomax =rhomin*ratrho;
    tmax   =pmax/rhomax;
    umax   = 0.0;

    % Dopo calcoliamo anche le altre variabili che ci servono.
    % Velocità suono, numero di Mach, entropia, energia interna, entalpia, 
    % pressione totale, temperatura totale
    % Tutte le variabili sono normalizzate rispetto alle quantità di
    % riferimento
    [amin,achmin,smin,emin,hmin,ptotmin,ttotmin] = calc_other_vars(pmin,rhomin,tmin,umin);
    [amax,achmax,smax,emax,hmax,ptotmax,ttotmax] = calc_other_vars(pmax,rhomax,tmax,umax);
    
    ncm = nc-1; %nc è il numero di celle => ncm è il numero di interfacce
    nhigh = nc/2+1; % Primo punto a destra della discontinuità
    diverg = 1.;
    grid(); % Crea la mesh computazionale, ovvero crea il vettore x  
    % in cui sono memorizzate le posizioni di tutti i centri di cella

    % call allocate_vars
    p(nhigh:nc)   = pmin;
    rho(nhigh:nc) = rhomin;
    t(nhigh:nc)   = tmin;
    p(1:nhigh-1)   = pmax;
    rho(1:nhigh-1) = rhomax;
    t(1:nhigh-1)   = tmax;
    u(1:nc)     = 0.0;
    
    % "calc_other_vars" prende come input anche vettori per ridarne in
    % output:
    [a,amach,s,e,h,ptot,ttot] = calc_other_vars(p,rho,t,u);
    [w1,w2,w3] = calc_cons(rho,u,e); % Calcola le variabili 
    % conservative -> presentazione 3 del corso
    [f1,f2,f3] = calc_fluxes(rho,u,p,e); % Calcola le variabili di flusso
    flow = f1;
    flht = f3;
    
elseif itest >=2 % GENERAL RIEMANN PROBLEMS
    
    tl = pl/rhol; %Temperatura a sinistra
    tr = pr/rhor; %Temperatura a destra
    
    [al,achl,sl,el,hl,ptotl,ttotl] = calc_other_vars(pl,rhol,tl,ul);
    [ar,achr,sr,er,hr,ptotr,ttotr] = calc_other_vars(pr,rhor,tr,ur);
    
    ncm=nc-1;
    nhigh=nc*x0+1;  %nhigh rappresenta l'indice della prima cella dopo la discontinuità iniziale
    diverg = 1.;
    kout = ka;
    grid
    % CALL ALLOCATE_VARS
    p(nhigh:nc)   = pr; %Finiscono in "r", sono a destra
    rho(nhigh:nc) = rhor;
    t(nhigh:nc)   = tr;
    u(nhigh:nc)   = ur;
    p(1:nhigh-1)   = pl; %Finiscono in "l", sono a sinistra
    rho(1:nhigh-1) = rhol;
    t(1:nhigh-1)   = tl;
    u(1:nhigh-1)   = ul;
    [a,amach,s,e,h,ptot,ttot] = calc_other_vars(p,rho,t,u);
    [w1,w2,w3] = calc_cons(rho,u,e);
    [f1,f2,f3] = calc_fluxes(rho,u,p,e);
    
    flow= f1;
    flht= f3;
end  % fine della condizione if

% inizializzazione a zero delle derivate spaziali e temporali di p, u, h
pxeno(1:nc) = 0.0d0;
uxeno(1:nc) = 0.0d0;
hxeno(1:nc) = 0.0d0;
ppxeno(1:nc)= 0.0d0;
hhxeno(1:nc)= 0.0d0;
ppt(1:nc)= 0.0d0;
ut(1:nc)= 0.0d0;
hht(1:nc)= 0.0d0;

end