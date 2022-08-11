% unstable version of pade splining that connects [2/1] degree pade approximants
%   over a set of knots with user-inputted (x_k, y_k), (x_k, y'_k)

% [7-15-2022] Current Issues: Discontinuity in the 2nd derivative, unsure how to
% prevent this with current implementation. Furthermore, certain (most..)
% combinations of parameters give Pade approximants that are DISCONTINUOUS
% between the knot locations in order to ensure f, f' conditions

x = linspace(1,5, 1000);
N = numel(x);
y = zeros([N 1]);

% Parameter format
% [y_1, y'_1
%  y_2, y'_2
%  y_3, y'_3]

params = [1,2; ... 
          2,1/10; ... 
          3,1];

knots = [1,2,5];

% compute each approximants coefficients and evaluates into y(k)
for k = 1:N
    % index(k) in knots of left knot
    i = roundTo(sum(x(k) >= knots), numel(knots)-1);
    
    a00 = params(i, 1);
    
    a01 = 1/params(i, 2);
    
    a10 = params(i,2) * (knots(i+1) - knots(i)) * (params(i+1, 1) - params(i, 1)) ...
        / ((params(i,2) * (knots(i+1) - knots(i)) - (params(i+1, 1) - params(i, 1))));
    
    a11 = ((params(i+1, 1) - params(i, 1)) - params(i,2)*(knots(i+1) - knots(i))).^2 ...
        / (params(i,2)* (params(i,2)*params(i+1,2)*((knots(i+1) - knots(i)).^2) ...
            - (params(i+1, 1) - params(i, 1)).^2 ));
    
    y(k) = ((x(k) - knots(i))*(x(k)-knots(i+1)) + a00*a01*(x(k)-knots(i+1)) + a10*a11*(x(k)-knots(i)) ... 
            + a00*a11*(x(k)-knots(i)) + a00*a01*a10*a11) ...
            / (a01*(x(k)-knots(i+1)) + a11*(x(k)-knots(i)) + a01*a10*a11); 
end

plot(x,y)
hold on
scatter(knots, params(:,1))