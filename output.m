
function output

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIgamma_PAR
global nc ncm k kout ka iord itest stab  % era USE NUM_PAR
global b c diverg c1 c2 x0 dx            %        era USE GEOM_PAR
global x                                          % era USE GEOM_VAR
global   time dt dtodx dtitest4 timemax timek    % era USE TIME_PAR
global p t u s rho a e amach ptot ttot flow flht h htot % era USE VARS
global w1 w2 w3 f1 f2 f3 phi1 phi2 phi3                    % era USE VARS
global p t u s rho a e amach ptot ttot flow flht h htot % era USE VARS
global w1 w2 w3 f1 f2 f3 phi1 phi2 phi3                    % era USE VARS
global  xsh10 xsh20 p1 rho1 u1 p2 rho2 u2 p3 rho3 u3 p4 rho4 u4 % era USE INIT_PAR
global  vsh1 vsh2 vsh3 vsh4 vsh5 vsh3l vsh3r vsh4l vsh4r % era USE INIT_PAR
global  xsh1 xsh2 xsh3 xsh4 xsh5 xsh3l xsh3r xsh4l xsh4r % era USE INIT_PAR
global pth tth uth sth rhoth ath eth amachth ptotth ttotth flowth flhtth hth % era USE THEORETICAL
global htotth r1th r2th r3th fplotth xth % era USE THEORETICAL
global xshth fthshu fthshd pu pd rhou tu su sd uu ud amu amd % era USE THEORETICAL
global htotu htotd flowu flowd flhtu flhtd r1u r1d r2u r2d r3u r3d % era USE THEORETICAL
global x1 x2 x3 x4 x5 % era USE THEORETICAL
global p1l p1r u1l u1r s1l s1r a1l a1r t1l t1r rho1l rho1r e1l e1r % era USE THEORETICAL
global p2l p2r u2l u2r s2l s2r a2l a2r t2l t2r rho2l rho2r e2l e2r % era USE THEORETICAL
global p3l p3r u3l u3r s3l s3r a3l a3r t3l t3r rho3l rho3r e3l e3r % era USE THEORETICAL
global p4l p4r u4l u4r s4l s4r a4l a4r t4l t4r rho4l rho4r e4l e4r % era USE THEORETICAL
global fplot % era USE OUTPUT_VAR

%     **********************
%10      CONTINUE

for n=1:nc
xth(n) = x(n);
end

for k=1:kout
    thesol(timek(k));
    if (itest <= 15)
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
                ik=n+ii
                ii = 1
                for nn = nc+ii:ik+1:-1
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
                for nn = nc+ii:ik+1:-1
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
                ik=n+ii
                ii = 3
                for nn = nc+ii:ik+1:-1
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
                for nn = nc+ii:ik+1:-1
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
                for nn = nc+ii:ik+1:-1
                    xth(nn)     = xth(nn-1);
                    pth(nn)     = pth(nn-1);
                    rhoth(nn)   = rhoth(nn-1);
                    tth(nn)     = tth(nn-1);
                    uth(nn)     = uth(nn-1);
                    amachth(nn) = amachth(nn-1);
                    sth(nn)     = sth(nn-1);
                end
                xth(ik)     = x4
                pth(ik)     = p3l
                rhoth(ik)   = rho3l
                tth(ik)     = t3l
                uth(ik)     = u3l
                amachth(ik) = u3l/sqrt(gamma*t3l)
                sth(ik)     = s3l
            end
            if (x(n-1) < x5) && (x(n) >= x5)
                ik=n+ii
                ii = 6
                for nn = nc+ii:ik+1:-1
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
        
% scrittura soluzione esatta, da inserire in seguito
%        if (timek(k) == 0.0)
%            OPEN(UNIT=10,FILE='theo_sol.dat',SHARED,STATUS='UNKNOWN')
%        else
%            OPEN(UNIT=10,FILE='theo_sol.dat',SHARED,STATUS='OLD',
%            &        ACCESS='APPEND')
%        end
%        WRITE(10,*)'TITLE="ITEST=',ITEST,'"'
%        WRITE(10,*)'VARIABLES="x" "p" "rho" "T" "u" "M" "S"'
%        WRITE(10,*)'ZONE T="t=',TIMEK(K),'" I=',NC+6,' F=POINT'
%        for n=1:nc+6
%            WRITE( 10,'(7E19.10)')XTH(N),PTH(N),RHOTH(N),TTH(N),UTH(N),
%            &                    AMACHTH(N),STH(N)
%        end
%        WRITE(10,*)
%        CLOSE(10)
    end
end
fprintf('INPUT:         1: Run');
fprintf('               2: Stop');
prompt = 'Run or Stop?'
inp = input(prompt);

if (inp == 1) 
    return
elseif (inp == 2)
    stop
end

%GO TO 10
end
        