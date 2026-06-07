function c = minmod(a,b)
% MINMOD clculates the value of the "minmod" slope limiter for the couple 
% of variables a and b.
        c = 0.0;
        if (a*b > 0)
            c = sign(a)*min(abs(a),abs(b));
        end
end

