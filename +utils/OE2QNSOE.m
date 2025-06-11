function QNS = OE2QNSOE(x)
% input: [a, e, i, Omega, omega, nu] (all in radians except a and e)

a     = x(1);
e     = x(2);
inc   = x(3);
RAAN  = x(4);
argp  = x(5);
nu    = x(6);

% All computations in radians
E     = 2 * atan( sqrt((1-e)/(1+e)) * tan(nu/2) );
M     = E - e * sin(E);
lam   = M + argp + RAAN * cos(inc);
ex    = e * cos(argp);
ey    = e * sin(argp);
ix    = inc;
iy    = RAAN * sin(inc);

QNS   = [a; lam; ex; ey; ix; iy];
end
