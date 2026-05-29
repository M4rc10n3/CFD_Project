function [utha,uthb,uthc,uthd,ptha,pthb,pthc,pthd,atha,athb,athc,athd,rhotha,rhothb,rhothc,rhothd,stha,sthb,sthc,sthd,alama,alamb,alamc,alamd,alamx,vshl,vshr,el,er] = rpsolver()

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIgamma_PAR
global nc ncm k kout ka iord itest stab % era USE NUM_PAR
global b c diverg c1 c2 x0 dx                     % era USE GEOM_PAR
global x                                          % era USE GEOM_VAR
global   time dt dtodx dtitest4 timemax timek     % era USE TIME_PAR
global p t u s rho a e amach ptot ttot flow flht h htot % era USE VARS
global w1 w2 w3 f1 f2 f3 phi1 phi2 phi3                    % era USE VARS
global pex pexsup pexinf pexpp                    % era USE BC
global  xsh10 xsh20 p1 rho1 u1 p2 rho2 u2 p3 rho3 u3 p4 rho4 u4 % era USE INIT_PAR
global  vsh1 vsh2 vsh3 vsh4 vsh5 vsh3l vsh3r vsh4l vsh4r % era USE INIT_PAR
global  xsh1 xsh2 xsh3 xsh4 xsh5 xsh3l xsh3r xsh4l xsh4r % era USE INIT_PAR
global  xcsf2                                            % era USE INIT_PAR

global pth tth uth sth rhoth ath eth amachth ptotth ttotth
global flowth flhtth hth htotth r1th r2th r3th fplotth xth
global xshth fthshu fthshd pu pd rhou tu su sd uu ud amu amd
global htotu htotd flowu flowd flhtu flhtd
global r1u r1d r2u r2d r3u r3d x1 x2 x3 x4 x5
global p1l p1r u1l u1r s1l s1r a1l a1r t1l t1r rho1l rho1r e1l e1r
global p2l p2r u2l u2r s2l s2r a2l a2r t2l t2r rho2l rho2r e2l e2r
global p3l p3r u3l u3r s3l s3r a3l a3r t3l t3r rho3l rho3r e3l e3r
global p4l p4r u4l u4r s4l s4r a4l a4r t4l t4r rho4l rho4r e4l e4r

%                    *** soluzione del problema di Riemann
%       *** attenzione: nel caso in cui venissero raggiunti i contorni, bisognerebbe cambiare qualcosa
atha=a(2); % "th" all'interno delle variabili indica che è il valore teorico
ptha=p(2);
utha=u(2);
stha=s(2);
athb=a(nc-1);
pthb=p(nc-1);
uthb=u(nc-1);
sthb=s(nc-1);
alama = 0.0; % (u-a)_a  % sono le caratteristiche lambda delle varie zone "a, b, c, d" del piano delle fasi
% Il piano delle fasi di cui si parla è quello con i fasci di espansione,
% la superficie di contatto e, appunto, la divisione in zone.
alamb = 0.0; % (u+a)_b
alamc = 0.0; % (u-a)_c
alamd = 0.0; % (u+a)_d
% Non è calcollata la "lambda" della superficie di contatto

% Una volta che si è capito il tipo del problema di Riemann, si possono
% trovare queste quantità (una sola delle due, sono quantità eventuali):
vshl = 0.0; % shock speed (fam. I)
vshr = 0.0; % shock speed (fam. III)

[pcd,ucd,rhothc,rhothd,el,er] = riemexact(atha,athb,ptha,pthb,utha,uthb);

pthc = pcd;
pthd = pcd;
uthc = ucd;
uthd = ucd;
athc = sqrt(gamma*pthc/rhothc);
athd = sqrt(gamma*pthd/rhothd);
rhotha = gamma*ptha/atha/atha;
rhothb = gamma*pthb/athb/athb;

if (el)
    alama = utha-atha;
    alamc = uthc-athc;
    sthc = stha;
else
    if (rhothc-rhotha) == 0.0 % evanescent shock wave
        vshl = 0.0;
    else
        % Si calcola la velocità della propagazione della discontinuità
        vshl = (rhothc*uthc-rhotha*utha)/(rhothc-rhotha);
    end
    % Si calcola anche l'entropia
    sthc = log(pthc)-gamma*log(rhothc);
end
if (er)
    alamd = uthd+athd;
    alamb = uthb+athb;
    sthd = sthb;
else
    if (rhothd-rhothb) == 0.0 % evanescent shock wave
        vshr = 0.0;
    else
        vshr = (rhothd*uthd-rhothb*uthb)/(rhothd-rhothb);
    end
    sthd = log(pthd)-gamma*log(rhothd);
end

alamx=uthc;
end