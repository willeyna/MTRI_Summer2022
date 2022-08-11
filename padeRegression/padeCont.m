syms x0 x1 x2 y0 y1 y2 p0 p1 p2 a00 a01 a10 a11 b00 b01 b10 b11 R1 R2 x

a00 = y0;
a01 = 1/p0;
a10 = p0*(x1-x0)*(y1-y0)/(p0*(x1-x0) - (y1-y0));
a11 = (y1 - y0 - p0*(x1-x0))^2 / (p0 * (p0*p1*(x1-x0)^2 - (y1-y0)^2));

b00 = y1;
b01 = 1/p1;
b10 = p1*(x2-x1)*(y2-y1)/(p1*(x2-x1) - (y2-y1));
b11 = (y2 - y1 - p1*(x2-x1))^2 / (p1 * (p1*p2*(x2-x1)^2 - (y2-y1)^2));

R1 = ((x-x0)*(x-x1) + a00*a01*(x-x1) + a10*a11*(x-x0) + a00*a11*(x-x0) + a00*a01*a10*a11) ... 
     / (a01*(x-x1) + a11*(x-x0) + a01*a10*a11);

R2 = ((x-x1)*(x-x2) + b00*b01*(x-x2) + b10*b11*(x-x1) + b00*b11*(x-x1) + b00*b01*b10*b11) ... 
     / (b01*(x-x2) + b11*(x-x1) + b01*b10*b11);

R1pp = simplify(diff(R1, 2), 5);
R2pp = simplify(diff(R2, 2), 5);

latex(R1pp)
latex(R2pp)

R1pp = subs(R1pp, x, x1);
R2pp = subs(R2pp, x, x1);

res = simplify(solve(R1pp == R2pp, p1),10);