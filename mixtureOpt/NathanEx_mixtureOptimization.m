% MattsEx_mixtureOptimization - Mixture optimization via reparameterization demo
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   NathanEx_mixtureOptimization()
%
% PURPOSE:
%   Modified version of
%   forge/mjwilder/SPARTA/MixtureOptimization/MattsEx_mixtureOptimization.m that
%   includes the Pade remapping variant + slightly updated code / comments
%
%   Provide a time comparison to ex_MixtureOptimization.  This function is
%   trying to identify points closes to tgt in a squared-error sense under the
%   constraint that they are non-negative and have coordinates that sum to one
%   (on the L1 ball).
%
%       f(x) = (x - tgt)^2 : ||x||_1 = r and x_i >= 0
%  
% INPUT:
%   None
% 
% OUTPUT:
%   Display to the terminal
%
% NOTES:
%   The default target point has been chosen such that the condition number of 
%   the Jacobian (a measure of the maximum function deviation) is poor relative 
%   to the isometric point of the nonlinear mapping.  This results in a more
%   poorly conditioned local structure that decreases optimization performance.
%   Setting the constant in front of beta to 1 shows how to correct this issue
%   with a translation of the domain.  Such translations must maintatin the
%   embedding, and I've added code to warn about constraint violations.
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid (if applicable),
%   and given in the correct order.
%-------------------------------------------------------------------------------
%}

N = 5;
%% 
tgt = abs(rand(N,1));

xIsom = sqrt(ones([N 1])/N);

%x0 is on the L2 unit ball
x0 = rand(N, 1);
x0 = x0/norm(x0,2);

beta = 0;      % Must induce non-negative points on L1 ball 

%%% Solve on the sphere through manOpt
problem.M = spherefactory(N);
problem.cost = @(x) costFunction(x,tgt,beta);
problem.egrad = @(x) gradFunction(x,tgt,beta);

[x1, xcost1, info1, options1] = rlbfgs(problem,x0);
dydx = diag(2*x1);
y = (x1.^2) + beta;
% condition number of mapping from sphere to simplex
disp("manOpt conditioning:")
disp(cond(dydx))


%%% Solve directly on the simplex through manOpt
problem.M = multinomialfactory(N);
problem.cost = @(x) norm(x-tgt,2)^2;
problem.egrad = @(x) 2*(x-tgt);

[x2, xcost2, info2, options2] = rlbfgs(problem,x0.^2);


% Set up spherical coordinates
alpha0 = cart2sphere( x0 );
alpha0 = alpha0(2:end);
objFun = @(a)objectiveFunction(a,tgt,beta);
opts = LBFGSOptions();
opts.verbose = true;
opts.debug = true;

%%% Solve with spherical coordinates and no shift on simplex (or beta--constant)

[alpha3,f3,g3,exitFlag3,output3] = LBFGS(objFun, alpha0, opts);
x3 = sphere2cart( [1;alpha3] ).^2;
x3 = x3 + beta;
disp(x3);

% Look at the conditioning of the internal mapping
[~,~, dxdAlpha,~] = objFun(alpha3);
disp("MTRI1 Conditioning:")
disp(cond(dxdAlpha));


%%% Solve with spherical coordinates + pade transformation

x0 = zeros([N 1]);
nIter = 200;
paramIter = nIter;

% objFun defined on the simplex as a function 
objFun = @(a) pointDist(a, tgt);
[x4, f4, g4, output4] = reparamMixtureMin(objFun, x0.^2, 1e-6, nIter, paramIter);
disp(x4);

% Display the per-iteration result
figure(1);
N1 = numel(output3.iter);
N2 = numel(output4.f);
semilogy(...
    [info1.iter],[info1.cost],'--bo',...
    [info2.iter],[info2.cost],'--mo', ... 
    0:N1, [output3.f0 [output3.iter.f]],'--ro', ...
    0:N2, [output4.f0, output4.f],'--ko');
xlabel('Iteration');
ylabel('Objective');
legend('Manopt Sphere', 'Manopt Simplex','MTRI Sphere', 'MTRI Pade');
xlim([0, 10])
prepareFigure();

% Display the per-iteration timing results
figure(2);
iterTimes3 = toc(output3.config.T0) - arrayfun(@toc, [output3.iter.T]);
iterTimes4 = toc(output4.T0) - arrayfun(@toc, [output4.T]);
loglog(...
    [info2.time], [info2.cost],'--mo',...
    [0 iterTimes3], [output3.f0 [output3.iter.f]],'--ro', ...
    [0 iterTimes4], [output4.f0 [output4.f]],'--ko');
xlabel('Time');
ylabel('Objective');
legend('manOpt', 'MTRI Spherical', 'MTRI Spherical w/ Pade');
grid('on');
prepareFigure();



% Squared error of x from tgt on the L1 shell [From spherical coordinates]
function [f,g,dxdAlpha,dMdAlpha] = objectiveFunction(alpha,tgt,beta)

% Default Values
P = numel(alpha)+1;

% Compute the weights at the current point
[MAlpha,dMdAlpha] = sphere2cart([1;alpha]);
x = MAlpha.^2;
x = x + beta;

% Evaluate the objective function in the original space
res = tgt-x;
f = norm(res)^2;
g = 2*res;

% Check for L1 ball and non-negativity
if any(x<0) || abs(sum(x)-1) > eps
    disp('Constraint Violation!');
end

% Propogate the gradient back into the transformed space
dMdAlpha(:,1) = [];          % Radius is 1 by definition
dxdAlpha = 2*MAlpha(:,ones(1,P-1)).*dMdAlpha;
g = -( g.'*dxdAlpha ).';

end


% Squared error of x from tgt on the L1 shell [From x on sphere]
function f = costFunction(x,tgt,beta)

% Evaluate the objective function in the original space
y = (x.^2) + beta;

res = tgt-y;
f = norm(res)^2;

% Check for L1 ball and non-negativity
if any(y<0) || abs(sum(y)-1) > eps
    disp('Constraint Violation!');
end

end

function g = gradFunction(x,tgt,beta)

y = (x.^2) + beta;

res = tgt-y;

g = -4*res.*x;

end

function [f, g] = pointDist(x, tgt)
    f = norm(x-tgt)^2;
    g = 2*(x-tgt);
end
