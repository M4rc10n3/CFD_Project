function thesol(timestep,utha,uthb,uthc,uthd,ptha,pthb,pthc,pthd,atha,athb,athc,athd,rhotha,rhothb,rhothc,rhothd,stha,sthb,sthc,sthd,alama,alamb,alamc,alamd,alamx,vshl,vshr,el,er)

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
global htotu htotd flowu flowd  flhtd
global r1u r1d r2u r2d r3u r3d x1 x2 x3 x4 x5
global p1l p1r u1l u1r s1l s1r a1l a1r t1l t1r rho1l rho1r e1l e1r
global p2l p2r u2l u2r s2l s2r a2l a2r t2l t2r rho2l rho2r e2l e2r
global p3l p3r u3l u3r s3l s3r a3l a3r t3l t3r rho3l rho3r e3l e3r
global p4l p4r u4l u4r s4l s4r a4l a4r t4l t4r rho4l rho4r e4l e4r

%     *** soluzione del problema di Riemann
%       *** attenzione: nel caso in cui venissero raggiunti i contorni, bisognerebbe cambiare qualcosa

if (el)
    x1 = x0 + alama*timestep;
    x2 = x0 + alamc*timestep;
else
    x1 = x0 + vshl *timestep;
    x2 = x1;
end
x3 = x0 + alamx*timestep;
if (er)
    x4 = x0 + alamd*timestep;
    x5 = x0 + alamb*timestep;
else
    x4 = x0 + vshr *timestep;
    x5 = x4;
end

for n=1:nc
    if (x(n) < x1)
        pth(n)=ptha;
        uth(n)=utha;
        sth(n)=stha;
        ath(n)=atha;
    end
    if (x(n) >= x1 && x(n) < x2)  % left expansion
        alam=(x(n)-x0)/timestep;
        rie =gg*atha + utha;
        ath(n)=(rie-alam)/gc;
        uth(n)= alam + ath(n);
        sth(n)=stha;
        pth(n)=ptha*(ath(n)/atha)^(1.0/gi);
    end
    if (x(n) >= x2 && x(n) < x3)
        pth(n)=pthc;
        uth(n)=uthc;
        sth(n)=sthc;
        ath(n)=athc;
    end
    if (x(n) >= x3 && x(n) < x4)
        pth(n)=pthd;
        uth(n)=uthd;
        sth(n)=sthd;
        ath(n)=athd;
    end
    if (x(n) >= x4 && x(n) < x5) % right expansion
        alam=(x(n)-x0)/timestep;
        rie    = gg*athb - uthb;
        ath(n) = (rie+alam)/gc;
        uth(n) = alam - ath(n);
        sth(n) = sthb;
        pth(n) = pthb*(ath(n)/athb)^(1.0/gi);
    end
    if(x(n) >= x5)
        pth(n)=pthb;
        uth(n)=uthb;
        sth(n)=sthb;
        ath(n)=athb;
    end
    tth(n)=ath(n)^2/gamma;
    rhoth(n)=pth(n)/tth(n);
    eth(n)  =gb*pth(n)+.5*rhoth(n)*uth(n)^2;
    amachth(n) = uth(n)/ath(n);
end

p1l=ptha;
p1r=pthc;
u1l=utha;
u1r=uthc;
s1l=stha;
s1r=sthc;
a1l=atha;
a1r=athc;
t1l=a1l^2/gamma;
t1r=a1r^2/gamma;
rho1l=p1l/t1l;
rho1r=p1r/t1r;
e1l=gb*p1l+.5*rho1l*u1l^2;
e1r=gb*p1r+.5*rho1r*u1r^2;

p2l=pthc;
p2r=pthd;
u2l=uthc;
u2r=uthd;
s2l=sthc;
s2r=sthd;
a2l=athc;
a2r=athd;
t2l=a2l^2/gamma;
t2r=a2r^2/gamma;
rho2l=p2l/t2l;
rho2r=p2r/t2r;
e2l=gb*p2l+.5*rho2l*u2l^2;
e2r=gb*p2r+.5*rho2r*u2r^2;

p3l=pthd;
p3r=pthb;
u3l=uthd;
u3r=uthb;
s3l=sthd;
s3r=sthb;
a3l=athd;
a3r=athb;
t3l=a3l^2/gamma;
t3r=a3r^2/gamma;
rho3l=p3l/t3l;
rho3r=p3r/t3r;
e3l=gb*p3l+.5*rho3l*u3l^2;
e3r=gb*p3r+.5*rho3r*u3r^2;

if (timestep == 0)
    p1l=ptha;
    p1r=ptha;
    p2l=ptha;
    p2r=pthb;
    p3l=pthb;
    p3r=pthb;
    u1l=utha;
    u1r=utha;
    u2l=utha;
    u2r=uthb;
    u3l=uthb;
    u3r=uthb;
    s1l=stha;
    s1r=stha;
    s2l=stha;
    s2r=sthb;
    s3l=sthb;
    s3r=sthb;
    a1l=atha;
    a1r=atha;
    a2l=atha;
    a2r=athb;
    a3l=athb;
    a3r=athb;
    t1l=a1l^2/gamma;
    t1r=a1r^2/gamma;
    t2l=a2l^2/gamma;
    t2r=a2r^2/gamma;
    t3l=a3l^2/gamma;
    t3r=a3r^2/gamma;
    rho1l=p1l/t1l;
    rho1r=p1r/t1r;
    rho2l=p2l/t2l;
    rho2r=p2r/t2r;
    rho3l=p3l/t3l;
    rho3r=p3r/t3r;
    e1l=gb*p1l+.5*rho1l*u1l^2;
    e1r=gb*p1r+.5*rho1r*u1r^2;
    e2l=gb*p2l+.5*rho2l*u2l^2;
    e2r=gb*p2r+.5*rho2r*u2r^2;
    e3l=gb*p3l+.5*rho3l*u3l^2;
    e3r=gb*p3r+.5*rho3r*u3r^2;
end

for n=1:nc
    xth(n) = x(n);
end
%     *** cerca gli estremi dell'espansione/urto della prima famiglia (x1 e x2), della superficie di contatto (x3) e di espansione/urto della terza famiglia (x4 e x5)
ii = 0;
ik=1;
xth(1) = x(1);
pth(ik)     = p1l;
rhoth(ik)   = rho1l;
tth(ik)     = t1l;
uth(ik)     = u1l;
amachth(ik) = u1l/sqrt(gamma*t1l);
sth(ik)     = s1l;
for n=2:nc
    if (x(n) >= x1) && (x(n-1) < x1)
        ik=n+ii;
        ii = 1;
        for nn = nc+ii:-1:ik+1
            xth(nn)     = xth(nn-1);
            pth(nn)     = pth(nn-1);
            rhoth(nn)   = rhoth(nn-1);
            tth(nn)     = tth(nn-1);
            uth(nn)     = uth(nn-1);
            amachth(nn) = amachth(nn-1);
            sth(nn)     = sth(nn-1);
        end
        xth(ik)     = x1;
        pth(ik)     = p1l;
        rhoth(ik)   = rho1l;
        tth(ik)     = t1l;
        uth(ik)     = u1l;
        amachth(ik) = u1l/sqrt(gamma*t1l);
        sth(ik)     = s1l;
    end
    if (x(n) >= x2) && (x(n-1) < x2)
        ik=n+ii;
        ii = 2;
        for nn = nc+ii:-1:ik+1
            xth(nn)     = xth(nn-1);
            pth(nn)     = pth(nn-1);
            rhoth(nn)   = rhoth(nn-1);
            tth(nn)     = tth(nn-1);
            uth(nn)     = uth(nn-1);
            amachth(nn) = amachth(nn-1);
            sth(nn)     = sth(nn-1);
        end
        xth(ik)     = x2;
        pth(ik)     = p1r;
        rhoth(ik)   = rho1r;
        tth(ik)     = t1r;
        uth(ik)     = u1r;
        amachth(ik) = u1r/sqrt(gamma*t1r);
        sth(ik)     = s1r;
    end
    if (x(n-1) < x3) && (x(n) >= x3)
        ik=n+ii;
        ii = 3;
        for nn = nc+ii:-1:ik+1
            xth(nn)     = xth(nn-1);
            pth(nn)     = pth(nn-1);
            rhoth(nn)   = rhoth(nn-1);
            tth(nn)     = tth(nn-1);
            uth(nn)     = uth(nn-1);
            amachth(nn) = amachth(nn-1);
            sth(nn)     = sth(nn-1);
        end
        xth(ik)     = x3;
        pth(ik)     = p2l;
        rhoth(ik)   = rho2l;
        tth(ik)     = t2l;
        uth(ik)     = u2l;
        amachth(ik) = u2l/sqrt(gamma*t2l);
        sth(ik)     = s2l;
        ik=n+ii;
        ii = 4;
        for nn = nc+ii:-1:ik+1
            xth(nn)     = xth(nn-1);
            pth(nn)     = pth(nn-1);
            rhoth(nn)   = rhoth(nn-1);
            tth(nn)     = tth(nn-1);
            uth(nn)     = uth(nn-1);
            amachth(nn) = amachth(nn-1);
            sth(nn)     = sth(nn-1);
        end
        xth(ik)     = x3;
        pth(ik)     = p2r;
        rhoth(ik)   = rho2r;
        tth(ik)     = t2r;
        uth(ik)     = u2r;
        amachth(ik) = u2r/sqrt(gamma*t2r);
        sth(ik)     = s2r;
    end
    if (x(n-1) < x4) && (x(n) >= x4)
        ik=n+ii;
        ii = 5;
        for nn = nc+ii:-1:ik+1
            xth(nn)     = xth(nn-1);
            pth(nn)     = pth(nn-1);
            rhoth(nn)   = rhoth(nn-1);
            tth(nn)     = tth(nn-1);
            uth(nn)     = uth(nn-1);
            amachth(nn) = amachth(nn-1);
            sth(nn)     = sth(nn-1);
        end
        xth(ik)     = x4;
        pth(ik)     = p3l;
        rhoth(ik)   = rho3l;
        tth(ik)     = t3l;
        uth(ik)     = u3l;
        amachth(ik) = u3l/sqrt(gamma*t3l);
        sth(ik)     = s3l;
    end
    if (x(n-1) < x5) && (x(n) >= x5)
        ik=n+ii;
        ii = 6;
        for nn = nc+ii:-1:ik+1
            xth(nn)     = xth(nn-1);
            pth(nn)     = pth(nn-1);
            rhoth(nn)   = rhoth(nn-1);
            tth(nn)     = tth(nn-1);
            uth(nn)     = uth(nn-1);
            amachth(nn) = amachth(nn-1);
            sth(nn)     = sth(nn-1);
        end
        xth(ik)     = x5;
        pth(ik)     = p3r;
        rhoth(ik)   = rho3r;
        tth(ik)     = t3r;
        uth(ik)     = u3r;
        amachth(ik) = u3r/sqrt(gamma*t3r);
        sth(ik)     = s3r;
    end
end

end