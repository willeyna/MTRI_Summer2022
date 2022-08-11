% function to fit a single [1/1] degree Pade approximant using a 6 dim 
% parametrization (x0,y0,xc,yc,xf,yf) using parametrizedPade

function fun = fitPade(x, y) 

objFun = @(a) sum((evaluatePade(x,a) - y).^2);

lb = [-Inf, -Inf, 0, 0, 0, 0];
beta0 = [min(x), min(y), 1, 0.5, max(x) - min(x), max(y) - min(y)];
disp(beta0)
beta = fmincon(objFun, beta0, [], [], [], [], lb, [], []);

fun = @(a) evaluatePade(x,beta);

end

function yhat = evaluatePade(x, beta)

% beta := (x0, y0, xc-x0, yc-y0, xf-x0, yf-y0) s.t. beta(2:end) > 0
yhat = parametrizedPade(x, beta(1) + beta(3), beta(2) + beta(4), ...
       beta(1), beta(2), beta(1) + beta(5), beta(2) + beta(6));

end
