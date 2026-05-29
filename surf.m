function surf(xx,itest)

global b c diverg c1 c2 x0 dx                     % era USE GEOM_PAR

if itest ~= 3
    surf = 1. + (diverg-1.)*xx/(c-b)
else
    c1 = 2.5
    c2 = 0.3
    x0 = 0.1
    surf = c1*(xx+x0)+c2/(xx+x0)
end

end