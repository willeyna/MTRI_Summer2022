% Calculations for the pade monotonicity and reformulation memo


syms x_c y_c x0 y0 x_f y_f a01 a11 b11 a02 a12 b12 y x m_c

%% With fixed endpoints and cont. s'' plugging in mc gives ONE Pade approximant

% (1.32)
mcfixed = (y_c*(1-y_c))/(x_c*(1-x_c));

% (1.29)
s1_fixed = (y_c^2 * x)/ (m_c * x_c^2 + (y_c - m_c*x_c)*x);
s1_fixed = simplify(subs(s1_fixed, m_c, mcfixed), 'Steps', 3);

% (1.30)
s2_fixed = (y_c*(1-y_c)- m_c*(1-x_c)*x_c + (m_c*(1-x_c) - y_c*(1 - y_c))*x) ...
                / (1-y_c - m_c*(1-x_c)*x_c + (m_c*(1-x_c) - (1-y_c))*x );
s2_fixed = simplify(subs(s2_fixed, m_c, mcfixed), 'Steps', 3);

disp(latex(s1_fixed))
disp(latex(s2_fixed))

%% Without fixed endpoints solve for both sets of Pade coefficients
% midpoint and endpoint contraints on both spline pieces
e1 = y0 == (a01 + a11*x0)/(1+b11*x0);
e2 = y_f == (a02 + a12*x_f)/(1+b12*x_f);
e3 = y_c == (a01 + a11*x_c)/(1+b11*x_c);
e4 = y_c == (a02 + a12*x_c)/(1+b12*x_c);
e5 = (a11 - a01*b11)/((1 + b11*x_c)^2) == (a12 - a02*b12)/((1 + b12*x_c)^2);
e6 = -2*((b11 * (a11 - a01*b11))/((1+b11*x_c)^3)) == -2*((b12 * (a12 - a02*b12))/(1+b12*x_c)^3);

eqns = [e1, e2, e3, e4, e5, e6];
solvefor = [a01, a11, b11, a02, a12, b12];
% solve in terms of mid and end points
sol = solve(eqns, solvefor);
% Subsitute in terms of endpts and mdpt to see both pieces are the same again
s1 =  (a01 + a11*x)/(1 + b11*x);
s2 =  (a02 + a12*x)/(1 + b12*x);

s1 = subs(s1, [a01, a11, b11], [sol.a01, sol.a11, sol.b11]);
s2 = subs(s2, [a02, a12, b12], [sol.a02, sol.a12, sol.b12]);

% final "simpified" version of s1,s2 as a function of endpoints and midpoint
s1 = simplify(s1);
s2 = simplify(s2);

s1_prime = simplify(diff(s1));
s2_prime = simplify(diff(s2));

disp(isequal(s1,s2));

disp(latex(s1))
disp(latex(s2))
disp(latex(s1_prime))
disp(latex(s2_prime))
%% Without fixed endpoints solve for ONE SET of Pade coefficients (proves same as before)
% midpoint and endpoint contraints on the spline
e1 = y0 == (a01 + a11*x0)/(1+b11*x0);
e2 = y_f == (a01 + a11*x_f)/(1+b11*x_f);
e3 = y_c == (a01 + a11*x_c)/(1+b11*x_c);

eqns = [e1, e2, e3];
solvefor = [a01, a11, b11];
% solve in terms of mid and end points
sol = solve(eqns, solvefor);

s1 =  (a01 + a11*x)/(1 + b11*x);

s1 = subs(s1, [a01, a11, b11], [sol.a01, sol.a11, sol.b11]);

% final "simpified" version of s1,s2 as a function of endpoints and midpoint
s1 = simplify(s1);

s1_prime = simplify(diff(s1));

disp(latex(s1))
%disp(latex(s1_prime))
%% For two Pade approximants solve for derivative conditions

e1 = (a01 + a11*x_c)/(1 + b11*x_c) == (a02 + a12*x_c)/(1 + b12*x_c);
e2 = (a11 - a01*b11)/(1 + b11*x_c)^2 == (a12 - a02*b12)/(1 + b12*x_c)^2;
e3 = (2*b11*(a01*b11-a11))/(b11*x_c+1)^3 == (2*b12*(a02*b12-a12))/(b12*x_c+1)^3;

sol = solve(e1,e2,e3,a01,a11,b11);
disp(isequal(sol.a01, a02))
disp(isequal(sol.a11, a12))
disp(isequal(sol.b11, b12))

