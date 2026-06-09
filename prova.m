%% Prova di implementazione funzione nel file script

fprintf('Il minimo tra 5 e 7 è %i \n', min(5,7))

function c = min(a, b)
    if a < b
        c = a;
    else
        c = b;
    end
end