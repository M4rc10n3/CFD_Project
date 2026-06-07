function [p,u,rhoc,rhod,el,er] = riemexact(aa,ab,pa,pb,ua,ub)

global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIgamma_PAR

rhoa = gamma*pa/aa^2;
capaa = 1.0/ge/rhoa;
rhob = gamma*pb/ab^2;
capab = 1.0/ge/rhob;

%Booleani che servono per capire se a sinistra e destra vi è un'onda d'urto
%o un'espansione
elold = false;
erold = false;
el    = true; 
er    = true; 

pold = 0.5*(pa+pb);

while ((el ~= elold) || (er ~= erold)) % "|" è l' "or" di MatLab
    elold = el;
    erold = er;
    it = 0;
    f =1.0;
    pnew = pold *2.0;
    while (abs(f) > 1.e-15) && (abs((pnew-pold)/pnew) > 1.e-15)
        if (it > 0)
            pold = pnew;
        end
        it = it+1;
        pp = pold + pold*1.e-5;
        pm = pold - pold*1.e-5;
        % Si usa un metodo tra i 4 tipi di problemi di Riemann che
        % possono verificarsi a seconda dei valori di "el" e "er".
        if (el && er)
            %fs e fe stanno per "funzione shock" e "funzione espansione",
            %si trovano a riga 112
            fp = fel(pp,ua,aa,pa)      - fer(pp,ub,ab,pb);
            fm = fel(pm,ua,aa,pa)      - fer(pm,ub,ab,pb);
            f = fel(pold,ua,aa,pa)     - fer(pold,ub,ab,pb);
        elseif (el)
            fp = fel(pp,ua,aa,pa)      - fsr(pp,ub,capab,pb);
            fm = fel(pm,ua,aa,pa)      - fsr(pm,ub,capab,pb);
            f  = fel(pold,ua,aa,pa)    - fsr(pold,ub,capab,pb);
        elseif (er)
            fp = fsl(pp,ua,aa,pa)      - fer(pp,ub,ab,pb);
            fm = fsl(pm,ua,capaa,pa)   - fer(pm,ub,ab,pb);
            f  = fsl(pold,ua,capaa,pa) - fer(pold,ub,ab,pb);
        else
            fp = fsl(pp,ua,capaa,pa)   - fsr(pp,ub,capab,pb);
            fm = fsl(pm,ua,capaa,pa)   - fsr(pm,ub,capab,pb);
            f  = fsl(pold,ua,capaa,pa) - fsr(pold,ub,capab,pb);
        end
        
        % L'algoritmo che si usa è il Newton-Raphson con la differenza
        % delle 2 curve di cui vogliamo trovare l'intersezione.
        ff = (fp-fm)/(pp-pm); % è la pendenza della retta differenza tra le curve
        pnew = pold - f/ff; % nuovo punto trovato con l'algoritmo
         
        if (pnew < 0.0d0)
            pnew = 1.e-6;
        end
    end
    
    p = pnew;
    if (el)
        u = ua + gg*aa*(1.d0-(p/pa)^gi);
    elseif (~el)
        u = ua - (p-pa)*sqrt(capaa/(p+pa/gc));
    end
    pc = p;
    pd = p;
    uc = u;
    ud = u;
    
    if (el)
        rhoc = rhoa*(pc/pa)^(1.0/gamma);
    elseif (~el)
        rhoc = rhoa*(1.0/gc+pc/pa)/(pc/pa/gc+1.0);
    end
    if (er)
        rhod = rhob*(pd/pb)^(1.0/gamma);
    elseif (~er)
        rhod = rhob*(1.0/gc+pd/pb)/(pd/pb/gc+1.0);
    end

    % Ora si fa un controllo sui valori trovati; se non è soddisfatto, sono cambiati 
    % i valori di "el" e "er" e il "while" continua.
    ac = sqrt(gamma*pc/rhoc);
    ad = sqrt(gamma*pd/rhod);
    
    ala = ua-aa;
    alc = uc-ac;
    
    ald = ud+ad;
    alb = ub+ab;
    
    el = false;
    er = false;
    
    if (ala <= alc)
        el = true;
    end
    if (ald <= alb)
        er = true;
    end
    
    % Questo "while" si ripeterà al massimo 1 volta, perché dopo la prima
    % iterazione si capisce se si tratta di fasci di espansione o onde
    % d'urto.

    pold = pc;
end
    %Le prossime funzioni provengono da pag. 29, 30 e 31 della
    %presentazione 6 del corso.
    function fel = fel(p,ua,aa,pa) %Across a left rarefaction
        %      use pigamma_par
        fel = ua - gg*aa*((p/pa)^gi-1.0);
    end

    function fer = fer(p,ub,ab,pb) %Across the right rarefaction
        %      use pigamma_par
        fer = ub + gg*ab*((p/pb)^gi-1.0);
    end

    function fsl = fsl(p,ua,capaa,pa) %Across a left shock
        %      use pigamma_par
        fsl = ua - (p-pa)*sqrt(capaa/(p+pa/gc));
    end

    function fsr = fsr(p,ub,capab,pb) %Across a right shock
        %      use pigamma_par
        fsr = ub + (p-pb)*sqrt(capab/(p+pb/gc));
    end

end
% ************************************************************************      