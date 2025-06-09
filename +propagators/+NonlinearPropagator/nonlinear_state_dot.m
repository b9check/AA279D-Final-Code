function [statedot] = nonlinear_state_dot(t, state, mu, Re, J2)
statedot = zeros(12, 1);

% CHIEF ECI INTEGRATOR --------------------------
% Unpack state
x0 = state(1);
y0 = state(2);
z0 = state(3);
vx0 = state(4);
vy0 = state(5);
vz0 = state(6);
r0_vec = [x0; y0; z0];
v0_vec = [vx0; vy0; vz0];
r0 = norm(r0_vec);

% Get perturbations
z2 = z0^2;
r2 = r0^2;
factor = -1.5 * J2 * mu * Re^2 / r0^7;
perturbations = factor * [ ...
    x0 * (1 - 5*z2/r2);
    y0 * (1 - 5*z2/r2);
    z0 * (3 - 5*z2/r2) ];

% Calculate stateerivative
statedot(1:3) = state(4:6);     
accel_kepler = -mu / r0^3 * r0_vec;
statedot(4:6) = accel_kepler + perturbations; 

% DEPUTY RTN INTEGRATOR
% Unpack state
x1 = state(7);
y1 = state(8);
z1 = state(9);
x1dot = state(10);
y1dot = state(11);
z1dot = state(12);

% Positional derivative = velocity
statedot(7:9) = state(10:12);

% Compute chief theta derivatives and r derivatives
h = cross(r0_vec, v0_vec);
theta0dot = norm(h) / r0^2;
r0dot = dot(r0_vec, v0_vec) / r0;
theta0dot2 = (-2*r0dot*theta0dot) / r0;

% Calculate theta/derivates and r/derivatives from chief state
statedot(10) = 2*theta0dot*y1dot + theta0dot2*y1 + theta0dot^2*x1 - (mu*(r0+x1)) / ((r0+x1)^2 + y1^2 + z1^2)^1.5 + mu/r0^2;
statedot(11) = -2*theta0dot*x1dot - theta0dot2*x1 + theta0dot^2*y1 - (mu*y1) / ((r0+x1)^2 + y1^2 + z1^2)^1.5;
statedot(12) = -(mu*z1) / ((r0+x1)^2 + y1^2 + z1^2)^1.5;

end