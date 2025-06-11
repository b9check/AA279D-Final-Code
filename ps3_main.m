%% 2B
alpha0 = [6771; 0.1005; deg2rad(51.64); deg2rad(257); deg2rad(0); deg2rad(30)]; % Chief
alpha1 = [6771; 0.1006; deg2rad(51.69); deg2rad(257.05); deg2rad(0.05); deg2rad(29.95)]; % Deputy

% Calculate initial state of chief
mu = 398600.4418;
chief_ECI0 = utils.OE2ECI(alpha0, mu);
r0 = chief_ECI0(1:3);
v0 = chief_ECI0(4:6);

% Calculate initial state [RTN position and velocity] of deputy
deputy_ECI0 = utils.OE2ECI(alpha1, mu);
r1 = deputy_ECI0(1:3);
v1 = deputy_ECI0(4:6);
RTN_0 = utils.ECI2RTN(chief_ECI0, deputy_ECI0);
rRTN = RTN_0(1:3);
vRTN = RTN_0(4:6);

% Solve for constants
a = alpha0(1);         
e = alpha0(2);
f = alpha0(6);
mu = 398600.4418;    
n = sqrt(mu/a^3);
h = norm(cross(r0,v0));

% Calculate values needed to get constants
fdot = h/norm(r0)^2;
eta = sqrt(1-e^2);
k = 1 + e*cos(f);
c = k*cos(f);
s = k*sin(f);

% Initial state
xbar0 = [ ...
   rRTN(1)/a;     
   vRTN(1)/(fdot*a); 
   rRTN(2)/a;  
   vRTN(2)/(fdot*a);    
   rRTN(3)/a;         
   vRTN(3)/(fdot*a); 
];


% Get constants
phi_inv = (1/eta^2) * [
    -3*s*(k+e^2)/k^2,            c - 2*e,     0,  -s*(k+1)/k,    0,               0;
    -3*(e + c/k),                -s,          0,  -(c*(k+1)/k + e), 0,            0;
     3*k - eta^2,                e*s,         0,   k^2,            0,               0;
    -3*e*s*(k+1)/k^2,           -2 + e*c,     eta^2, -e*s*(k+1)/k,  0,               0;
     0,                          0,           0,   0,              eta^2*cos(f),   -eta^2*sin(f);
     0,                          0,           0,   0,              eta^2*sin(f),    eta^2*cos(f)
];

Constants = phi_inv * xbar0


%% 2C
% Propogate
step_size = 1*(pi/180);
f_change = 15*2*pi;
[f_vals, states] = get_YA_states(alpha0, Constants, step_size, f_change, n, mu, h);

x = states(:, 1);
vx = states(:, 2);
y = states(:, 3);
vy = states(:, 4);
z = states(:, 5);
vz = states(:, 6);

% Remove normalization
x = x*a;
vx = vx*a*fdot;
y = y*a;
vy = vy*a*fdot;
z = z*a;
vz = vz*a*fdot;

% Figure 1: RTN Position
figure(1)
subplot(2,2,1)
plot(x, y)
xlabel('R (m)')
ylabel('T (m)')
grid on

subplot(2,2,2)
plot(x, z)
xlabel('R (m)')
ylabel('N (m)')
grid on

subplot(2,2,3)
plot(y, z)
xlabel('T (m)')
ylabel('N (m)')
grid on

subplot(2,2,4)
plot3(x, y, z)
xlabel('R (m)')
ylabel('T (m)')
zlabel('N (m)')
grid on

% Figure 2: RTN Velocity
figure(2)
subplot(2,2,1)
plot(vx, vy)
xlabel('v_R (m/s)')
ylabel('v_T (m/s)')
grid on

subplot(2,2,2)
plot(vx, vz)
xlabel('v_R (m/s)')
ylabel('v_N (m/s)')
grid on

subplot(2,2,3)
plot(vy, vz)
xlabel('v_T (m/s)')
ylabel('v_N (m/s)')
grid on

subplot(2,2,4)
plot3(vx, vy, vz)
xlabel('v_R (m/s)')
ylabel('v_T (m/s)')
zlabel('v_N (m/s)')
grid on


function [f_vals, states] = get_YA_states(alpha0, constants, step_size, f_change, n, mu, h)
%GET_YA_STATES  Propagate YA state transition matrix instead of manual phi
%   alpha0: 6×1 [ρ0; e; … ; f0]
%   constants: 6×1 initial YA state
%   step_size: step size in true anomaly (rad)
%   f_change: total true‐anomaly change (rad)
%   n: mean motion (rad/s)
%   mu: gravitational parameter (m^3/s^2)
%   h: specific angular momentum (m^2/s)

    % extract elements
    e         = alpha0(2);
    f_initial = alpha0(6);
    f_final   = f_initial + f_change;
    f_vals    = f_initial:step_size:f_final;
    states    = zeros(length(f_vals), 6);

    % compute mean‐anomaly time history
    E = 2*atan( sqrt((1-e)./(1+e)) .* tan(f_vals/2) );
    E = unwrap(E);
    M = E - e.*sin(E);
    M0 = M(1);
    t  = (M - M0) ./ n;       % time since f0

    % recover semi-major axis from n: n = sqrt(mu / a^3)
    a = (mu / n^2)^(1/3);

    % build STM function handle using STMYA
    Options.f0 = f_initial;
    Options.e  = e;
    Options.a  = a;
    Options.mu = mu;
    Phi_fun    = STMYA([], Options);

    t0 = 0;
    for i = 1:length(f_vals)
        tf = t(i);
        Phi = Phi_fun(tf, t0);
        state_i = Phi * constants;
        states(i, :) = state_i';
    end
end


function A = EOMYA(System, Options)

% Scaled Version:
% A = @(f) [0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1;3/(1+Options.e*cos(f)) 0 0 0 -2 0;0 0 0 2 0 0;0 0 -1 0 0 0];

% e = Options.e;
% Unscaled Version (independent variable is true anomaly):
% A = @(f) [zeros(3,3) eye(3,3);
%           (3+e*cos(f))/(1+e*cos(f)) -2*e*sin(f)/(1+e*cos(f)) 0 2*e*sin(f)/(1+e*cos(f)) 2 0;
%           2*e*sin(f)/(1+e*cos(f)) e*cos(f)/(1+e*cos(f)) 0 -2 2*e*sin(f)/(1+e*cos(f)) 0;
%           0 0 -1/(1+e*cos(f)) 0 0 2*e*sin(f)/(1+e*cos(f))];

% Unscaled Version (independent variable is time):

function A = EOM(t,f0,e,a,mu)

M0 = TrueToMeanAnomaly(f0,e);
n = sqrt(mu/a^3);
f = MeanToTrueAnomaly(t*n+M0,e);
h = sqrt(mu*a*(1-e^2));
fdot = sqrt(mu/a^3)*(1+e*cos(f))^2/(1-e^2)^(3/2);
fddot = -2*mu*e*sin(f)*(1+e*cos(f))^3/(a*(1-e^2))^3;
k = mu/h^(3/2);

A = [zeros(3,3) eye(3,3);
     (2*k*fdot^(3/2)+fdot^2)  fddot 0 0 2*fdot 0;
     -fddot (-k*fdot^(3/2)+fdot^2) 0 -2*fdot 0 0;
     0 0 -k*fdot^(3/2) 0 0 0];

end

A = @(t) EOM(t,Options.f0,Options.e,Options.a,Options.mu);

end

function Phi = STMYA(System, Options)

function Phi = STM(tf,t0,f0,e,a,mu)

M0 = utils.TrueToMeanAnomaly(f0,e);
n = sqrt(mu/a^3);
ff = utils.MeanToTrueAnomaly((tf-t0)*n+M0,e);
rhof = 1+e*cos(ff);
rho0 = 1+e*cos(f0);
sf = rhof*sin(ff);
s0 = rho0*sin(f0);
cf = rhof*cos(ff);
c0 = rho0*cos(f0);
spf = cos(ff) + e*cos(2*ff);
sp0 = cos(f0) + e*cos(2*f0);
cpf = -sin(ff) - e*sin(2*ff);
cp0 = -sin(f0) - e*sin(2*f0);
k2 = sqrt(mu)/sqrt(a*(1-e^2))^3;
% M = TrueToMeanAnomaly(ff,e);
% M0 = TrueToMeanAnomaly(f0,e);
% J = (M-M0)/sqrt((1-e^2)^3);
J = k2*(tf-t0);

Phixy = [sf 0 (2-3*e*sf*J) -cf;
       cf*(1+1/rhof) 1 -3*rhof^2*J sf*(1+1/rhof);
       spf -0 -3*e*(spf*J+sf/rhof^2) -cpf;
       -2*sf 0 -3*(1-2*e*sf*J) 2*cf-e]*...
      [-3*s0*(1/rho0+e^2/rho0^2) 0 c0-2*e -s0*(1+1/rho0);
       -3*e*s0*(1/rho0+1/rho0^2) 1-e^2 e*c0-2 -e*s0*(1+1/rho0);
       3*rho0+e^2-1 0 e*s0 rho0^2;
       3*(c0/rho0+e) 0 s0 c0*(1+1/rho0)+e]*...
       1/(1-e^2);

rhof0 = 1+e*cos(ff-f0);
cf0 = rhof0*cos(ff-f0);
sf0 = rhof0*sin(ff-f0);
Phiz = 1/rhof0*[cf0 sf0;-sf0 cf0];

Phi = [Phixy(1:2,1:2) zeros(2,1) Phixy(1:2,3:4) zeros(2,1);
       zeros(1,2) Phiz(1,1) zeros(1,2) Phiz(1,2);
       Phixy(3:4,1:2) zeros(2,1) Phixy(3:4,3:4) zeros(2,1);
       zeros(1,2) Phiz(2,1) zeros(1,2) Phiz(2,2)];

Phi = [1/rhof*eye(3,3) zeros(3,3);
       k2*e*sin(ff)*eye(3,3) k2*rhof*eye(3,3)]*...
       Phi*...
      [rho0*eye(3,3) zeros(3,3);
       -e*sin(f0)*eye(3,3) 1/(k2*rho0)*eye(3,3)];

end

Phi = @(tf,t0) STM(tf,t0,Options.f0,Options.e,Options.a,Options.mu);

end




%% 2e
[a_0,e_0,i_0,W_0,w_0,f_0] = deal( ...
    6771, ...                      % a1
    0.1005, ...                    % e1
    deg2rad(51.64), ...            % i1
    deg2rad(257), ...              % O1
    deg2rad(0), ...                % w1
    deg2rad(30) ...                % f1
);
[a_1,e_1,i_1,W_1,w_1,f_1] = deal( ...
    6771, ...
    0.1006, ...
    deg2rad(51.69), ...
    deg2rad(257.05), ...
    deg2rad(0.05), ...
    deg2rad(29.95) ...
);


E_0 = 2*atan2(sqrt((1 - e_0)/(1 + e_0))*tan(f_0/2), 1);
M_0 = E_0 - e_0*sin(E_0);

E_1 = 2*atan2(sqrt((1 - e_1)/(1 + e_1))*tan(f_1/2), 1);
M_1 = E_1 - e_1*sin(E_1);

delta_a = (a_1 - a_0)/a_0;
delta_lambda = (M_1 + w_1) - (M_0 + w_0) + (W_1 - W_0)*cos(i_0);
delta_ex = e_1*cos(w_1) - e_0*cos(w_0);
delta_ey = e_1*sin(w_1) - e_0*sin(w_0);
delta_ix = i_1 - i_0;
delta_iy = (W_1 - W_0)*sin(i_0);
qns_init = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy];



%% 2f
quasi_elements = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy]';
alpha0 = [6771; 0.1005; deg2rad(51.64); deg2rad(257); deg2rad(0); deg2rad(30)]; % Chief

% [a, ex, ey, u, ix, iy

[f_vals2, states] = linear_ecc_mapping(quasi_elements, alpha0, 15*pi*2, 1*(pi/180));

x2 = states(:, 1);
y2 = states(:, 2);
z2 = states(:, 3);
vx2 = states(:, 4);
vy2 = states(:, 5);
vz2 = states(:, 6);


function [f_vals, states] = linear_ecc_mapping(quasi_elements, alpha0, f_change, step_size)
% Solve for constants
a = alpha0(1);         
e = alpha0(2);
w = alpha0(5);
i = alpha0(3);
mu = 398600.4418;    
n = sqrt(mu/a^3);

% Calculate values needed to get constants
eta = sqrt(1-e^2);
ex = e*cos(w);
ey = e*sin(w);

% F vals
f_initial = alpha0(6);
f_final = f_initial + f_change;
f_vals = f_initial:step_size:f_final;

% Calculate t matrix
E = 2*atan2(sqrt((1 - e)./(1 + e)).*tan(f_vals./2), 1);
E = unwrap(E);
M = E - e.*sin(E);               
M0 = M(1);                       
ts = (M - M0) ./ n;    

% Initialize states
states = zeros(length(f_vals), 6);

for j = 1:length(f_vals)
% assume scalars k, eta, n, t, ex, ey, u, i are already defined
f = f_vals(j);
k = 1 + e*cos(f);
kp = -e*sin(f);
u = f + w;
t = ts(j);

% Compute B matrix entries (x component shown as example)
bx1 = 1 / k + (3 / 2) * kp * n / eta^3 * t;
bx2 = -kp / eta^3;
bx3 = (1 / eta^3) * (ex * (k - 1) / (1 + eta) - cos(u));
bx4 = (1 / eta^3) * (ey * (k - 1) / (1 + eta) - sin(u));
bx6 = kp / eta^3 * cot(i);
by1 = -(3 / 2) * k * n / eta^3 * t;
by2 = k / eta^3;
by3 = (1 / eta^2) * ((1 + 1 / k) * sin(u) + ey / k + k / eta * (ey / (1 + eta)));
by4 = -(1 / eta^2) * ((1 + 1 / k) * cos(u) + ex / k + k / eta * (ex / (1 + eta)));
by6 = (1 / k - k / eta^3) * cot(i);
bz5 = (1 / k) * sin(u);
bz6 = -(1 / k) * cos(u);
% Velocities
bxd1 = kp / 2 + (3 / 2) * k^2 * (1 - k) * n / eta^3 * t;
bxd2 = k^2 / eta^3 * (k - 1);
bxd3 = k^2 / eta^3 * (eta * sin(u) + ey * (k - 1) / (1 + eta));
bxd4 = -k^2 / eta^3 * (eta * cos(u) + ex * (k - 1) / (1 + eta));
bxd6 = -k^2 / eta^3 * (k - 1) * cot(i);
byd1 = -(3 / 2) * k * (1 + k * kp * n / eta^3 * t);
byd2 = k^2 / eta^3 * kp;
byd3 = (1 + k^2 / eta^3) * cos(u) + ex / eta^2 * k  * (1 + k / eta * (1 - k) / (1 + eta));
byd4 = (1 + k^2 / eta^3) * sin(u) + ey / eta^2 * k  * (1 + k / eta * (1 - k) / (1 + eta));
byd6 = -(1 + k^2 / eta^3) * kp * cot(i);
bzd5 = cos(u)+ex;
bzd6 = sin(u)+ey;
% Assemble B matrix
B = [bx1, bx2, bx3, bx4, 0, bx6;
     by1, by2, by3, by4, 0, by6;
     0,    0,    0,    0, bz5, bz6;
     bxd1, bxd2, bxd3, bxd4, 0, bxd6;
     byd1, byd2, byd3, byd4, 0, byd6;
     0,    0,    0,    0, bzd5, bzd6];
M = [ ...
    a*eta^2*eye(3),    zeros(3); ...
    zeros(3),        (a*n)/eta*eye(3) ...
];

Phi = M * B;

state = Phi*quasi_elements;
states(j,:) = state;

end

end

