function c = minmod(a,b)
%MINMOD Summary of this function goes here
%   Detailed explanation goes here
        c = 0.0;
        if (a*b > 0)
            c = sign(a)*min(abs(a),abs(b));
        end
end

