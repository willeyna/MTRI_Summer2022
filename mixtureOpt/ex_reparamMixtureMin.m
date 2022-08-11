%-------------------------------------------------------------------------------
%   Low dimensionality example of using simplex reparametrization to optimize
%       over micture constraints
%   Objective function is defined to be the distance between a user-defined
%   point tgt and a point on the regular simplex of the size of tgt. 
%
%   See reparamMixtureMin for definition of inputs
%-------------------------------------------------------------------------------

% user-defined example parameters -------------------
tgt = [1;1.1;1;.1; .2];
% overall normalization of x0 is not important due to parametrization via angles
x0 = rand(size(tgt));
objFun = @(a) pointDist(a, tgt);
gradTol = sqrt(eps);
nIter = 100;
paramIter = 10;
%----------------------------------------------------

[x, f, g, output] = reparamMixtureMin(objFun, x0, gradTol, nIter, paramIter);
disp("Optima:")
disp(x);
disp("||Gradient|| of Reparametrized Function at Optima:")
disp(norm(g));

function [f, g] = pointDist(x, tgt)
    f = norm(x-tgt)^2;
    g = 2*(x-tgt);
end