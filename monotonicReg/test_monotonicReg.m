% test_monotonicRegression - Test suite
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   This code provides the test suite.  It can be run through the MOxUnit
%   unit testing framework.
%
% PURPOSE:
%   This function provides the unit tests.  New bug reports for the code should
%   not be closed without ensuring that a test in the suite both fails before
%   the fix, and passes after it.
%
%-------------------------------------------------------------------------------
%}

function test_monotonicRegression

try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions = localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;

end

% makes sure calculated gradient matches with the numerical gradient
function test_gradient
      
    fun = fitTestSpline(100, 0.1, 4, 10);

    xlin = linspace(0,1,101);  
    [~, gDerived] = fun(xlin);
    gNumeric = numericalGradient(fun, xlin, 1E-6).';
    % collpase diagonal matrix
    gNumeric = sum(gNumeric,2);

    % numeric gradient is an approximation also so this is not a ground truth
    % comparison, thus removing the end points from the comparison should remove 
    % the points with the greatest error 
    assertElementsAlmostEqual(gDerived(2:end-1),gNumeric(2:end-1),'relative',1e-3)
end

% returns a calculated spline on a normal cdf with mean 0.5 and std 0.1
%   training data has N pts and normal random error around cdf of size eps
%   spline has order ord and nKnots=number of knots ignoring multiplicity 
function [fun, x, y] = fitTestSpline(N, eps, ord, nKnots)

    % test data generation
    x = rand(N, 1);
    knots = linspace(0,1,nKnots);
    model_fun = @(a) cdf('Normal', a, 0.5, 0.1);
    y = model_fun(x) + normrnd(0, eps, [N, 1]);
    
    fun = monotonicRegression(x, y, ord, knots, true);

end