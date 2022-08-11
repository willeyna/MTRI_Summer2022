% example script for fitPade. Fits sigmoidal data using a [1/1] Pade approx

N = 100;
eps = 0.2;

x = -3 + 6*rand([N 1]);
x = sort(x);
mu = 1;
sigma = .5;

% model_fun = @(a) cdf('Normal', a, mu, sigma);
model_fun = @(a) sigmoid(a, mu, sigma);
y = model_fun(x) + normrnd(0,eps,[N, 1]);

fun = fitPade(x,y);

xlin = linspace(min(x), max(x), 1000);
ylin = model_fun(xlin);

hold on
scatter(x,y);
plot(x, fun(x));
plot(xlin, ylin, 'Color', 'Black' ,'LineStyle', '--')