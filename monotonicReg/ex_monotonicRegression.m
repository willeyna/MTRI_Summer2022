function ex_monotonicRegression()

% x selection and variation
N = 75;
eps = 0.3;
ord = 4;
% for plotting
xlin = linspace(0,1,1000);

% for non-cdf fits
x = rand(N, 1);
%% sin function fitting (test of monotonicity and regression spline)
knots = linspace(0,1,5);

sinfig = figure('Name', 'Sin(3pi/4 x)');
deriv_sinfig = figure('Name', 'd/dx Sin(3pi/4 x)');

model_fun = @(a) sin((3*pi/4)*a) + 0.5;
model_deriv = @(a) (3*pi/4) * cos((3*pi/4)*a);

y = model_fun(x) + normrnd(0,eps,[N, 1]);

fun = monotonicRegression(x, y, ord, knots, false);

[f,g]=fun(xlin);

figure(deriv_sinfig)
hold on
plot(xlin, model_deriv(xlin), 'Color', 'red', 'LineStyle', '--');
plot(xlin, g);

figure(sinfig)
hold on
scatter(x, y);
plot(xlin, model_fun(xlin), 'Color', 'black', 'LineStyle', '--');
plot(xlin, f);

%% normal cdf function fitting (test of cdf fitting)
knots = linspace(-1,1,50);

% for plotting
xlin = linspace(-1,1,10000);

y = rand([N 1]);

normfig1 = figure('Name', 'Normal cdf mu=.5, sigma=.1');
deriv_normfig1 = figure('Name', 'Normal pdf mu=.5, sigma=.1');

mu = 0;
sigma = .3;

model_fun = @(a) cdf('Normal', a, mu, sigma);
model_deriv = @(a) pdf('Normal', a, mu, sigma);

% pull from inverse cdf and add Gaussian noise
x = icdf('Normal', y, mu, sigma);
y = y + normrnd(0,eps,[N, 1]);

fun = monotonicRegression(x, y, ord, knots, true, 1, 1e3);
[f,g]=fun(xlin);

figure(deriv_normfig1)
hold on
plot(xlin, model_deriv(xlin), 'Color', 'red', 'LineStyle', '--');
plot(xlin, g);

figure(normfig1)
hold on
scatter(x, y);
plot(xlin, model_fun(xlin), 'Color', 'black', 'LineStyle', '--');
plot(xlin, f);

%% normal cdf function fitting (steeper ascent)
% example where number of knots > n 
knots = linspace(-.05,.05,100);

% for plotting
xlin = linspace(-.05,.05,10000);

y = rand([N 1]);

normfig2 = figure('Name', 'Normal cdf mu=.5, sigma=.1');
deriv_normfig2 = figure('Name', 'Normal pdf mu=.5, sigma=.1');

mu = 0;
sigma = .01;

model_fun = @(a) cdf('Normal', a, mu, sigma);
model_deriv = @(a) pdf('Normal', a, mu, sigma);

% pull from inverse cdf and add Gaussian noise
x = icdf('Normal', y, mu, sigma);
y = y + normrnd(0,eps,[N, 1]);

fun = monotonicRegression(x, y, ord, knots, true, 2, 1e5);
[f,g]=fun(xlin);

figure(deriv_normfig2)
hold on
plot(xlin, model_deriv(xlin), 'Color', 'red', 'LineStyle', '--');
plot(xlin, g);

figure(normfig2)
hold on
scatter(x, y);
plot(xlin, model_fun(xlin), 'Color', 'black', 'LineStyle', '--');
plot(xlin, f);


%% exponential dist cdf
knots = linspace(0,10,100);
% for plotting
xlin = linspace(0,10,1000);

y = rand([N 1]);

expfig = figure('Name', 'Exponential cdf, sigma=.1');
deriv_expfig = figure('Name', 'Exponential pdf, sigma=.1');

lambda = 1.5;

model_fun = @(a) cdf('exp', a, lambda);
model_deriv = @(a) pdf('exp', a, lambda);

x = icdf('exp', y, lambda);
y = y + normrnd(0,eps,[N, 1]);

fun = monotonicRegression(x, y, ord, knots, true, 2, 5e5);

[f,g]=fun(xlin);
plot(xlin, f);
plot(xlin, g);

figure(deriv_expfig)
hold on
plot(xlin, model_deriv(xlin), 'Color', 'red', 'LineStyle', '--');
plot(xlin, g);

figure(expfig)
hold on
scatter(x, y);
plot(xlin, model_fun(xlin), 'Color', 'black', 'LineStyle', '--');
plot(xlin, f);

%% fitting ASA data exposition data as in Wang and Li (2008) fig 5
car_fig = figure('Name', 'Difference from max mpg vs hp');
carData = importdata("isotonic_ex.txt");
x = carData(:,1);
y = carData(:,2);

xmin = min(x);
xmax = max(x);

knots = linspace(xmin, xmax, 300);
xlin = linspace(xmin, xmax, 1000);

fun = monotonicRegression(x, y, ord, knots, false, 1, 1e5);

f = fun(xlin);

figure(car_fig)
hold on
scatter(x,y);
alpha(.3)
plot(xlin, f);
end