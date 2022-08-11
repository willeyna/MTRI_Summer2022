% reparamMixtureMin- Uses simplex reparametrization to find optima of mixture
% problem
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   x = reparamMixtureMin(fun, x0)
%   [x, f, g, output] = reparamMixtureMin(fun, x0, gradTol, nIter, paramIter)
%
% PURPOSE:
%   Minimize a mixture problem (min f(x)| sum(x)=1 x_i >= 0 for all i)
%   Minimizes over spherical coordinates, which are remapping element-wise
%   by Pade splines s.t. the search in spherical coordinates occurs near the 
%   spherical-to-cartesian transformation's isometry point, improving
%   conditioning
%
% INPUT:
%   fun       - Lambda function f(x) to be minimized where x is defined to be on 
%               the standard simplex
%   x0        - [N 1] Initial guess of the minimizer
%   gradTol   - Absolute tolerance of gradient at which to end the linesearch
%   maxIter     - Number of iterations that the L-BFGS routine will take in total
                between all reparametrizations
%   paramIter - Number of iterations the minimization will take
%               before reparametrizing the simplex (and thus lose Hessian
%               information)
%
% OUTPUT:
%   x         - Location of located optima on simplex 
%   f         - Function evaluation at optima
%   g         - Gradient of the objective function in terms of re-mapped
%               spherical coordinates (Note: NOT dfdx; dfdx !=0 necessirly at 
%               optima due to mixture constraints)
%   output    - Structure storing basic per-iteration details 
%
% ASSUMPTIONS: 
%   fun is of the form [f, g] = fun(x) x adhering to mixture constraints
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order.
%-------------------------------------------------------------------------------
%}
function [xi, fAlpha, dfdAlpha, output] = reparamMixtureMin(fun, x0, varargin)

% begin timing as soon as the script is called
T0 = tic();

[gradTol, maxIter, paramIter] = parseInputs(varargin{:});

% point of isometry from sphere to simplex 
N = numel(x0);
xIsom = ones([N 1])/N;
alphaIsom = simplex2sphere(xIsom);

% remove radial compoenent to reduce minimizer computations
alphaIsom = alphaIsom( 2:end );

% LBFGS options to be pre-defined
opts = LBFGSOptions;
opts.maxIter = paramIter;
opts.maxEval = 2*paramIter;
opts.verbose = true;
opts.gradTol = gradTol;
opts.debug = true;

% initialize variables for loop
dfdAlpha = 1;
xi = x0;

% basic per-iteration output for diagnostics
output = struct('x', [], ...
                'f', [], ...
                'g', [], ...
                'nIter', 0, ...
                'f0', fun(x0),  ...
                'T', [], ...
                'T0', T0);


% first order minimization loop
objFun = @(a)simplexParam(a, fun, xIsom, xi);
disp(xi)
while norm(dfdAlpha) > gradTol && output.nIter < maxIter

    [alpha, fAlpha, dfdAlpha, ~, iterOutput] = LBFGS(objFun, alphaIsom, opts);

    % find true point on simplex of current step
    xi = simplexRemap(sphere2simplex(alpha), xIsom, xi);
    % change obj param s.t. isometry pt maps to current step location
    objFun = @(a)simplexParam(a, fun, xIsom, xi);

    % Update output information %
    output.T = [output.T, iterOutput.iter.T];
    output.x = [output.x, iterOutput.iter.x];
    output.f = [output.f, iterOutput.iter.f];
    output.g = [output.g, iterOutput.iter.g];
    output.nIter = output.nIter + iterOutput.numIter;
end


end

% evaluates objective function at a point in spherical coordinates using 
%   sphere2simplex and the element-wise Pade mapping defined by (xc, yc)
%   to produce x(alpha)
function [fAlpha, dfdAlpha, kappa] =  simplexParam(alpha, fun, xc, yc)

    [yAlpha, dydAlpha] = sphere2simplex(alpha);

    [x, dxdY] = simplexRemap(yAlpha, xc, yc);

    [fAlpha, dfdx] = fun(x);

    dxdAlpha = dxdY * dydAlpha;
    dfdAlpha = (dfdx.' * dxdAlpha).';

    % condition number of mapping from sphere to shifted simplex point
    kappa = cond(dxdAlpha);
end

% Hadamard product squared of sphere2cart
function [yAlpha, dydAlpha] = sphere2simplex(alpha)

    [mAlpha, dMdAlpha] = sphere2cart([1;alpha]);
    dMdAlpha(:,1) = []; % removes dxdr since dr==0

    yAlpha = mAlpha.^2; 
    dydAlpha = 2*mAlpha .* dMdAlpha;
end


% remap point on the standard simplex using coordinate-wise Pade splines
function [x, dxdY] = simplexRemap(y, xc, yc)

    N = numel(y);

    % -- AD-HOC FIX FOR PADE PINNING --

%     % stops remapping too close to yc=1 where Pade splines are ill-conditioned
    yc = roundTo(yc, .9);
    yc = roundTo(yc, .1, false);

    % -------------------------- 

    % uses pade spline to shift isometry point to xc
    [P, dpdY] = padeVectorized(y, [xc, yc]);

    % normalize to stay on the simplex
    b = sum(P);
    x = P/b;

    dxdY = (1/b^2 * (b* eye(N) - P.*ones(N))) * (dpdY.*eye(N));
end

% evaluate pade splines and derivative row-wise on inputs
function [S,G] = padeVectorized(x, controls)

d = numel(x);
S = zeros([d 1]);
G = zeros([d 1]);

for i = 1:d
    [S(i,:), G(i,:)] = oneNodePadeSplines(x(i,:), [controls(i,:)]);
end
end

% element-wise sqrt sends point on simplex to a point in cartesian coordinates 
%   with a radius of one
function alpha = simplex2sphere(x)
    alpha = cart2sphere(sqrt(x));
end

function [gradTol, maxIter, paramIter] = parseInputs(varargin)
    gradTol = sqrt(eps);
    maxIter = 100;
    paramIter = 10;

switch nargin
    case 0
    case 1
        gradTol = varargin{1};
    case 2 
        gradTol = varargin{1};
        maxIter = varargin{2};
    case 3
        gradTol = varargin{1};
        maxIter = varargin{2};
        paramIter = varargin{3};
    otherwise
        error([mfilename,':invalidInput'],...
            'Unexpected number of inputs');
end

end