%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   fun = monotonicRegression(x, y, ord, knots)
%   fun = monotonicRegression(x, y, ord, knots, true)
%   fun = monotonicRegression(x, y, ord, knots, false, penaltyOrder, lambda)
%   [yhat, ghat] = fun(x)
%
% PURPOSE:
%   Fits a monotonically increasing regression spline of specified order to the
%       data (x,y) using I-Splines under an L2 loss
%   Allows for traditional monotone quadratic regression spline methods with
%       carefully selected knots or a more flexibile approach using finite
%       difference penalties (P-spline) and a smoothness parameter
%   Minimizes coefficients using LBFGS routine
%
% INPUT:
%   x             - [N 1] x values of input data
%   y             - [N 1] y values of input data
%   ord           - Order of I-splines to be fit 
%   knots         - Sequence of knots [a, t_1, ..., t_r, b] where a,b are the 
%                 values at the boundary of the fit region 
%   cdf           - [Boolean] Specifies whether the spline should adhere to
%                   cumulative density function restrictions (range=[0,1])
%                   Default false
%   penaltyOrder  - [Integer] Order of finite difference penalty applied to
%                   regression 
%   lambda        - Hyperparameter determining strength of penalty in the least
%                   squares problem
%
% OUTPUT:
%   fun     - Function that evaluates the spline and it's derivative at x
%   yhat    - Ordinates of spline at x
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%-------------------------------------------------------------------------------
%}
function fun = monotonicRegression(x, y, ord, knots, varargin)

% MATLAB's Bspline calculator expects monotonic x values
[x,sortIdx] = sort(x);
y = y(sortIdx);

% spcol expects x to be unique values if only calling 0th derivative
xUnique = unique(x);
% stores the original x multiplicities 
Q = histcounts(x, [xUnique; Inf]);

[cdf, dOrd, lambda, nIter] = parseInputs(varargin{:});

% get endpoints
a = knots(1);
b = knots(end);
% get order of B-splines
deg = ord - 1;

% a,b multiplicity ord = deg+1; "extended partition" boundary conditions
knots = [a*ones([1, deg]), knots, b*ones([1, deg])];

% iSplines evaluated at the unique x points
iUnique =  calculate_iSpline(knots, ord, xUnique, cdf);

% builds the full evaluation of all i Splines including x values with mult > 1
startIdx = [0, cumsum(Q)] + 1;
for i = 1:length(xUnique)
    I(startIdx(i):startIdx(i) + Q(i)-1, :) = repmat(iUnique(i,:), Q(i), 1);
end

% length of coefficient vector
N = size(I,2);

%%% Minimize for spline coefficeints %%%

P = differenceMatrix(N, dOrd);
if cdf == false
    % makes sure that y-intercept isn't penalized
    P(end, :) = [];
end

% minimize penalized least squares
objFun = @(a) penalizedLSTSQ(a, I, y, P, lambda);
beta0 = ones([N 1]);
beta0 = beta0/sum(beta0);

if cdf
    % minimize over the probability simplex with the manOpt package
    problem.M = multinomialfactory(N);
    problem.cost = @(beta) norm(I*beta - y, 2).^2 + ... 
                                lambda*(beta.' * (P'*P) * beta);
    problem.egrad = @(beta) 2.*(I.' * I * beta - I.'*y + lambda.*(P.'*P)*beta);
    
    options = struct("maxiter", nIter, "verbosity", 0, "Delta_bar", 1/N);
    [beta, ~, info] = trustregions(problem, beta0, options);
    gradients = [info.gradnorm];
    g = gradients(1,end);

    disp("Minimization complete with |grad|:")
    disp(g)
else
    % minimize with positive boundary constraint
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
    lb = [zeros(N-1, 1); -Inf];

    beta = fmincon(objFun, beta0, [], [], [], [], lb, [], [], options);
end

fun = @(a) evaluate_iSpline(beta, knots, ord, a, cdf);

end

% returns I-Spline basis + constant matrix evaluated at x
function [I, M] = calculate_iSpline(knots, ord, x, cdf)

% uses cumsum of B splines in De Boor (1976) (4.11) to calculate I-Splines
B = spcol(knots, ord, x);
I = 1 - cumsum(B,2);

% gets rid of redundant 0-spline
if cdf
    I(:,end) = [];
else
    % adds constant basis function to allow spline to have non-zero y-intercept
    I(:,end) = 1;
end

% B splines are an order lower and so need one less boundary knot
B_knots = knots(2:end-1);
% create B splines and transform into M-splines (derivative of iSpline by ftoc)
B = spcol(B_knots, ord-1, x);
% Leeuw (2017) equation (3)
integral_normalization = (ord-1)./(B_knots(ord:end) - B_knots(1:end-(ord-1)));
M = integral_normalization .* B;

end

% returns I-Spline basis + constant matrix evaluated at x
function [f, g] = evaluate_iSpline(beta, knots, ord, x, cdf)

    % MATLAB's Bspline calculator expects monotonic x values
    [x,sortIdx] = sort(x);

    [I, M] = calculate_iSpline(knots, ord, x, cdf);
    f = I*beta;
    if cdf
        g = M*beta;
    else
        % get rid of constant coefficient for derivative
        g = M*beta(1:end-1);
    end
, 2000
    % re-sort outputs to match original input
    f(sortIdx) = f;
    g(sortIdx) = g;
end

% contrained and penalizaed LSTSQ s.t. each coeff is non-negative
function [f,g] = penalizedLSTSQ(beta, N, y, P, lambda)
    f = norm(N*beta - y, 2).^2 + lambda*(beta.' * (P'*P) * beta);
    % wrong and not actually being used yet due to reparamMixtureMin grad broken
    g = 2.*(N.' * N * beta - N.'*y + lambda.*(P.'*P)*beta);
end

% recursively calculate higher order difference matrices from the Identity
function P = differenceMatrix(nparam, ord)
    
    % by default 0th order gives Identity matrix as Penalty 
    P = eye(nparam);
    for k = 1:ord
        for i = 1:(nparam-k)
            P(i,:) = P(i,:) - P(i+1,:);
        end
        P(end, :) = [];
    end

end

% this function should be handled better with an options strucutre 
%   but I ran out of time to add it
function [cdf, dOrd, lambda, nIter] = parseInputs(varargin)

cdf = false;
dOrd = 0;
lambda = 0;
nIter = 1000;

switch nargin
    case 1
        cdf = varargin{1};
    case 3
        cdf = varargin{1};
        dOrd = varargin{2};
        lambda = varargin{3};
    case 4
        cdf = varargin{1};
        dOrd = varargin{2};
        lambda = varargin{3};
        nIter = varargin{4};
    otherwise
        error([mfilename,':invalidInput'],...
            'Unexpected number of inputs');
end



end
