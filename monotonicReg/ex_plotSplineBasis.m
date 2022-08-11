function ex_plotSplineBasis()

%% INPUT
% interior knots include all but the first and last element of knots
%    implicit here is that the end points are inputted as knots
knots = linspace(5,7,5);
% desired degree/ order of the *B-Splines* (1+deg is degree of I-Splines)
deg = 2;
%%

xlin = linspace(knots(1),knots(end),1000);
% get endpoints
a = knots(1);
b = knots(end);
% get order of B-splines
ord = deg + 1;

% number of interior knots-- consistent with notation in Leeuw
r = numel(knots) - 2;

% a,b multiplicity m (= deg+1); "extended partition"
B_knots = [a*ones([1, deg]), knots, b*ones([1, deg])];
% gets B/N splines s.t. sum splines = 1 for a given x (for visualization)
B = spcol(B_knots, ord, xlin);

% Leeuw (2017) equation (3)
integral_normalization = ord./(B_knots(1+ord:end) - B_knots(1:end-(ord)));
M = integral_normalization .* B; 

% clamped boundary conditions enforced by multiplicity of a and b
B1p_knots = [a, B_knots, b];
% calculates higher order splines for calculation of I-Splines
B1p = spcol(B1p_knots, ord+1, xlin);

% uses cumsum described in De Boor (1976) (4.11) to calculate I-Splines
I = 1 - cumsum(B1p,2);
%gets rid of useless 0-spline
I = I(1:end,1:end-1);

bfig = figure('Name', 'B-Splines');
mfig = figure('Name', 'M-Splines');
ifig = figure('Name', 'I-Splines');


figure(bfig)
%plot dashed lines at knots
for i = 1:r
    plot([knots(i+1), knots(i+1)], [0, 1], 'Color', 'k', 'LineStyle', '--');
    hold on
end
title("B-Spline Basis")
plot(xlin, B)

figure(mfig)
%plot dashed lines at knots
for i = 1:r
    plot([knots(i+1), knots(i+1)], [0, 1], 'Color', 'k', 'LineStyle', '--');
    hold on
end
title("Integral-Normalized B-Spline Basis")
plot(xlin, M)


figure(ifig)
%plot dashed lines at knots
for i = 1:r
    plot([knots(i+1), knots(i+1)], [0, 1], 'Color', 'k', 'LineStyle', '--');
    hold on
end
title("I-Spline Basis")
plot(xlin, I)
   
end