function c = iterstar(a, b)
% ITERSTAR  Iterative solver used in the transonic case of split.m
%   Solves: c = exp((b - c) / a)  by fixed-point iteration.
%   Input:  a, b  (scalar doubles)
%   Output: c     (converged solution)

global gb  % gamma combination: 1/(gamma-1)

try1 = 0.5 * a / gb;

for kip = 1:500
    try2 = exp((b - try1) / a);
    if abs(try2 - try1) < 1.e-5
        c = try2;
        return
    end
    try1 = try2;
end

error('ITERSTAR: iteration failed to converge after 500 steps')

end
