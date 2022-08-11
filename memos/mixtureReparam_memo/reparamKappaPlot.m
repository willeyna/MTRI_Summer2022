% Primitive script to plot a 1D slice of the parametrization Jacobian's
%   condition number 

function reparamKappaPlot

res = 1000;

% point of isometry from sphere to simplex 
N = 3;
xIsom = ones([N 1]) / N;

xlin = linspace(0.01, 1-0.01, res);

tgt = xIsom;
fun = @(x)dealVec({1,numel(x)}, [norm(x-tgt)^2 ; 2*(x-tgt)] );  


xi = zeros([N 1]);
kappa = zeros(res, 1);

xpin = rand([N 1]);
xpin = xpin/sum(xpin);

xpin = xIsom;

for i = 1:res

    % 1-d slice along [x1;a;a;...;a] on S_N ; y
    yi = [xlin(i) ; (1-xlin(i))/(N-1) * ones([N-1 1])];

    t = simplexRemap(yi, xIsom, xpin);
    % store the first element of the transformed simplex coordinate 
    xi(i) = t(1);

    % simplexParam needs alpha as an input (rather than yi)
    alpha = cart2sphere(sqrt(yi));        
    alpha = alpha(2:end);

    [~, ~, k] = simplexParam(alpha, fun, xIsom, xpin);
    kappa(i) = k;
    
end

plot(xlin, kappa);
title("Conditioning of Map from Spherical onto Reparametrized Simplex")
ylabel("Condition Number of Jacobian")
xlabel("x(1)")
hold on 
% plot(xi, kappa);
% xline(xpin(1),'--r');
% xline(xIsom(1),'--b');

end


%%  Bringing in functions directly from reparamMixtureMin.m 

function [fAlpha, dfdAlpha, kappa] =  simplexParam(alpha, fun, xc, yc)

    [yAlpha, dydAlpha] = sphere2simplex(alpha);

    [x, dxdY] = simplexRemap(yAlpha, xc, yc);

    [fAlpha, dfdx] = fun(x);

    dxdAlpha = dxdY * dydAlpha;
    dfdAlpha = (dfdx.' * dxdAlpha).';

    % condition number of mapping from sphere to shifted simplex point
    Skappa = cond(dydAlpha);
    kappa = cond(dxdAlpha);

end

% Hadamard product squared of sphere2cart
function [yAlpha, dydAlpha] = sphere2simplex(alpha)

    [MAlpha, dMdAlpha] = sphere2cart([1;alpha]);
    dMdAlpha(:,1) = []; % gets rid of dxdr 

    yAlpha = MAlpha.^2; 
    dydAlpha = 2*MAlpha .* dMdAlpha;
end


% remap point on the standard simplex using coordinate-wise Pade splines
function [x, dxdY] = simplexRemap(y, xc, yc)

    N = numel(y);

    % uses pade spline to shift isometry point to xc
    [P, dpdY] = padeVectorized(y, [xc, yc]);

    % normalize to stay on the simplex
    b = sum(P);
    x = P/b;

    dxdY = (1/b^2 * (b* eye(N) - P.*ones(N))) * (dpdY.*eye(N));
end

% evaluate pade splines row-wise on inputs
function [S,G] = padeVectorized(x, controls)

d = numel(x);
S = zeros([d 1]);
G = zeros([d 1]);

for i = 1:d
    [S(i,:), G(i,:)] = oneNodePadeSplines(x(i,:), controls(i,:));
end

end