%% Problem set 1

%% 2b
% Constants
mu = 3.986004418e14;  
R = 6378137;          
J2 = 1.08262668e-3;   

% Convert OE to ECI
chief_OE = [6771e3, 0.0005, deg2rad(51.64), deg2rad(257), deg2rad(0), deg2rad(30)];
deputy_OE = [6771e3, 0.0015, deg2rad(52.14), deg2rad(258), deg2rad(0.5), deg2rad(25)];
chief_ECI = utils.OE2ECI(chief_OE, mu);
deputy_ECI = utils.OE2ECI(deputy_OE, mu);

% Display
disp(chief_ECI(1:3))
disp(chief_ECI(4:6))
disp(deputy_ECI(1:3))
disp(deputy_ECI(4:6))


%% 2c
% Setup
T = 2*pi*sqrt(chief_OE(1)^3 / mu); % Period
tspan = 0:1:10*T;  
state0 = chief_ECI;

% Call ode113 without J2 perturbations
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t, state, mu, R, J2, perturbated);
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
[t_out1, x_out1] = ode113(odefun, tspan, state0, opts);
r_mag1 = vecnorm(x_out1(:,1:3), 2, 2);
v_mag1 = vecnorm(x_out1(:,4:6), 2, 2);

% Call ode113 with J2 perturbations
perturbated = true;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t, state, mu, R, J2, perturbated);
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
[t_out2, x_out2] = ode113(odefun, tspan, state0, opts);
r_mag2 = vecnorm(x_out2(:,1:3), 2, 2);
v_mag2 = vecnorm(x_out2(:,4:6), 2, 2);

% Plot
figure(1)
plot3(x_out1(:,1), x_out1(:,2), x_out1(:,3))
xlabel('X_{ECI} (m)')
ylabel('Y_{ECI} (m)')
zlabel('X_{ECI} (m)')

figure(2)
plot3(x_out2(:,1), x_out2(:,2), x_out2(:,3))
xlabel('X_{ECI} (m)')
ylabel('Y_{ECI} (m)')
zlabel('X_{ECI} (m)')

% Position magnitude plots
figure(3)
plot(t_out1/T, r_mag1, 'b', 'DisplayName', 'No Perturbation');
hold on
plot(t_out2/T, r_mag2, 'r--', 'DisplayName', 'With Perturbation');
xlabel('Orbital Periods')
ylabel('||r|| (m)')
legend
grid on

% Velocity magnitude plots
figure(4)
plot(t_out1/T, v_mag1, 'b', 'DisplayName', 'No Perturbation');
hold on
plot(t_out2/T, v_mag2, 'r--', 'DisplayName', 'With Perturbation');
xlabel('Orbital Periods')
ylabel('||v|| (m/s)')
legend
grid on

% Separate subplots for clarity
figure(5)
subplot(2,1,1)
plot(t_out1/T, r_mag1)
xlabel('Orbital Periods')
ylabel('||r|| (m)')
title('Position Magnitude (No Perturbation)')
grid on

subplot(2,1,2)
plot(t_out1/T, v_mag1)
xlabel('Orbital Periods')
ylabel('||v|| (m/s)')
title('Velocity Magnitude (No Perturbation)')
grid on

figure(6)
subplot(2,1,1)
plot(t_out2/T, r_mag2)
xlabel('Orbital Periods')
ylabel('||r|| (m)')
title('Position Magnitude (With Perturbation)')
grid on

subplot(2,1,2)
plot(t_out2/T, v_mag2)
xlabel('Orbital Periods')
ylabel('||v|| (m/s)')
title('Velocity Magnitude (With Perturbation)')
grid on


%% 2d
% Propagate Keplerian
T = 2*pi*sqrt(chief_OE(1)^3 / mu); % Period
tspan = 0:1:10*T;  
state0 = chief_OE;

% Call ode113 without J2 perturbations
perturbated = false;
odefun = @(t, state) propagators.KeplerianPropagator.getStatedot(t, state, mu, R, J2, perturbated);
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
[t_out3, x_out3] = ode113(odefun, tspan, state0, opts);

x_RTNs = zeros(length(tspan), 6);
for i = 1:length(tspan)
    ECI_1 = x_out1(i, :);
    OE_2 = x_out3(i, :);
    ECI_2 = utils.OE2ECI(OE_2, mu)';
    x_RTN = utils.ECI2RTN(ECI_1, ECI_2)';
    x_RTNs(i, :) = x_RTN;
end


% Figure 7: Position error RTN projections
figure(7)
subplot(2,2,1)
plot(x_RTNs(:,1), x_RTNs(:,2), '.')
xlabel('R error (m)')
ylabel('T error (m)')
grid on

subplot(2,2,2)
plot(x_RTNs(:,1), x_RTNs(:,3), '.')
xlabel('R error (m)')
ylabel('N error (m)')
grid on

subplot(2,2,3)
plot(x_RTNs(:,2), x_RTNs(:,3), '.')
xlabel('T error (m)')
ylabel('N error (m)')
grid on

subplot(2,2,4)
plot3(x_RTNs(:,1), x_RTNs(:,2), x_RTNs(:,3), '.')
xlabel('R error (m)')
ylabel('T error (m)')
zlabel('N error (m)')
grid on

% Figure 8: Velocity error RTN projections
figure(8)
subplot(2,2,1)
plot(x_RTNs(:,4), x_RTNs(:,5), '.')
xlabel('R_{dot} error (m/s)')
ylabel('T_{dot} error (m/s)')
grid on

subplot(2,2,2)
plot(x_RTNs(:,4), x_RTNs(:,6), '.')
xlabel('R_{dot} error (m/s)')
ylabel('N_{dot} error (m/s)')
grid on

subplot(2,2,3)
plot(x_RTNs(:,5), x_RTNs(:,6), '.')
xlabel('T_{dot} error (m/s)')
ylabel('N_{dot} error (m/s)')
grid on

subplot(2,2,4)
plot3(x_RTNs(:,4), x_RTNs(:,5), x_RTNs(:,6), '.')
xlabel('R_{dot} error (m/s)')
ylabel('T_{dot} error (m/s)')
zlabel('N_{dot} error (m/s)')
grid on



% Figure 9: RTN Position Error vs Time
figure(9)
subplot(3,1,1)
plot(t_out3/T, x_RTNs(:,1))
ylabel('R error (m)')
grid on

subplot(3,1,2)
plot(t_out3/T, x_RTNs(:,2))
ylabel('T error (m)')
grid on

subplot(3,1,3)
plot(t_out3/T, x_RTNs(:,3))
xlabel('Orbital Periods')
ylabel('N error (m)')
grid on

% Figure 10: RTN Velocity Error vs Time
figure(10)
subplot(3,1,1)
plot(t_out3/T, x_RTNs(:,4))
xlabel('Time (s)')
ylabel('R_{dot} error (m/s)')
grid on

subplot(3,1,2)
plot(t_out3/T, x_RTNs(:,5))
ylabel('T_{dot} error (m/s)')
grid on

subplot(3,1,3)
plot(t_out3/T, x_RTNs(:,6))
xlabel('Orbital Periods')
ylabel('N_{dot} error (m/s)')
grid on


%% 2e
% PROPAGATE PERTURBED KEPLERIAN ELEMENTS
T = 2*pi*sqrt(chief_OE(1)^3 / mu); % Period
tspan = 0:1:10*T;  
state0 = chief_OE;

% Call ode113 without J2 perturbations
perturbated = true;
odefun = @(t, state) propagators.KeplerianPropagator.getStatedot(t, state, mu, R, J2, perturbated);
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
[t_out4, x_out4] = ode113(odefun, tspan, state0, opts);

% CALCULATE DESIRED QUANTITIES
% Preallocate arrays
n = size(x_out3,1);
e_vecs1     = zeros(n,3);
h_vecs1     = zeros(n,3);
energy_vec1 = zeros(n,1);
e_vecs2     = zeros(n,3);
h_vecs2     = zeros(n,3);
energy_vec2 = zeros(n,1);

% Loop ove*r time
for k = 1:n
    % --- First dataset (x_out3) ---
    oe1 = x_out3(k,:);
    state1 = utils.OE2ECI(oe1, mu)';   % 6x1
    r1 = state1(1:3);
    v1 = state1(4:6);
    e_vecs1(k,:)     = (1/mu)*((norm(v1)^2 - mu/norm(r1))*r1 - dot(r1,v1)*v1);
    h_vecs1(k,:)     = cross(r1, v1);
    energy_vec1(k)   = norm(v1)^2/2 - mu/norm(r1);
    
    % --- Second dataset (x_out4) ---
    oe2 = x_out4(k,:);
    state2 = utils.OE2ECI(oe2, mu)';   % 6x1
    r2 = state2(1:3);
    v2 = state2(4:6);
    e_vecs2(k,:)     = (1/mu)*((norm(v2)^2 - mu/norm(r2))*r2 - dot(r2,v2)*v2);
    h_vecs2(k,:)     = cross(r2, v2);
    energy_vec2(k)   = norm(v2)^2/2 - mu/norm(r2);
end

% Compute magnitudes
e_mag1 = vecnorm(e_vecs1, 2, 2);
h_mag1 = vecnorm(h_vecs1, 2, 2);
e_mag2 = vecnorm(e_vecs2, 2, 2);
h_mag2 = vecnorm(h_vecs2, 2, 2);


% PLOT
% Figure 11: Keplerian elements over time
figure(11)
subplot(6,1,1)
plot(t_out3/T, x_out3(:,1))
ylabel('a (m)')
grid on

subplot(6,1,2)
plot(t_out3/T, x_out3(:,2))
ylabel('e')
grid on

subplot(6,1,3)
plot(t_out3/T, rad2deg(x_out3(:,3)))
ylabel('i (deg)')
grid on

subplot(6,1,4)
plot(t_out3/T, rad2deg(x_out3(:,4)))
ylabel('\Omega (deg)')
grid on

subplot(6,1,5)
plot(t_out3/T, rad2deg(x_out3(:,5)))
ylabel('\omega (deg)')
grid on

subplot(6,1,6)
plot(t_out3/T, rad2deg(mod(x_out3(:,6), 2*pi)))
ylabel('\nu (deg)')
xlabel('Orbital Periods')
grid on

% Figure 12: Perturbed Keplerian elements over time
figure(12)
subplot(6,1,1)
plot(t_out4/T, x_out4(:,1))
ylabel('a (m)')
grid on

subplot(6,1,2)
plot(t_out4/T, x_out4(:,2))
ylabel('e')
grid on

subplot(6,1,3)
plot(t_out4/T, rad2deg(x_out4(:,3)))
ylabel('i (deg)')
grid on

subplot(6,1,4)
plot(t_out4/T, rad2deg(x_out4(:,4)))
ylabel('\Omega (deg)')
grid on

subplot(6,1,5)
plot(t_out4/T, rad2deg(x_out4(:,5)))
ylabel('\omega (deg)')
grid on

subplot(6,1,6)
plot(t_out4/T, rad2deg(mod(x_out4(:,6), 2*pi)))
ylabel('\nu (deg)')
xlabel('Orbital Periods')
grid on

% Figure 13: unperturbed case
figure(13)
subplot(3,1,1)
hold on
plot(t_out3/T, e_mag1)
plot(t_out4/T, e_mag2)
legend('Unperterbed', 'J2 Perturbed')
ylabel('|e|')
grid on

subplot(3,1,2)
hold on
plot(t_out3/T, h_mag1)
plot(t_out4/T, h_mag2)
ylabel('|h| (m^2/s)')
grid on

subplot(3,1,3)
hold on
plot(t_out3/T, energy_vec1)
plot(t_out4/T, energy_vec2)
ylabel('\epsilon (J/kg)')
xlabel('Orbital Periods')
grid on


%% 2f
% PROPAGATE PERTURBED KEPLERIAN ELEMENTS
T = 2*pi*sqrt(chief_OE(1)^3 / mu); % Period
tspan = 0:1:10*T;  
state0 = chief_OE;

% Call ode113 without J2 perturbations
odefun = @(t, state) propagators.KeplerianPropagator.getMeanStatedot(t, state, mu, R, J2);
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
[t_out5, x_out5] = ode113(odefun, tspan, state0, opts);


% FIGURE 14: Unperturbed vs. J2 Mean
figure(14)
subplot(6,1,1)
plot(t_out3/T, x_out3(:,1), 'g', ...
     t_out5/T, x_out5(:,1), 'r--')
ylabel('a (m)')
legend('Unperturbed','J2 Mean')
grid on
set(gca, 'XTickLabel', [])

subplot(6,1,2)
plot(t_out3/T, x_out3(:,2), 'g', ...
     t_out5/T, x_out5(:,2), 'r--')
ylabel('e')
grid on
set(gca, 'XTickLabel', [])

subplot(6,1,3)
plot(t_out3/T, rad2deg(x_out3(:,3)), 'g', ...
     t_out5/T, rad2deg(x_out5(:,3)), 'r--')
ylabel('i (deg)')
grid on
set(gca, 'XTickLabel', [])

subplot(6,1,4)
plot(t_out3/T, rad2deg(mod(x_out3(:,4), 2*pi)), 'g', ...
     t_out5/T, rad2deg(mod(x_out5(:,4), 2*pi)), 'r--')
ylabel('\Omega (deg)')
grid on
set(gca, 'XTickLabel', [])

subplot(6,1,5)
plot(t_out3/T, rad2deg(mod(x_out3(:,5), 2*pi)), 'g', ...
     t_out5/T, rad2deg(mod(x_out5(:,5), 2*pi)), 'r--')
ylabel('\omega (deg)')
grid on
set(gca, 'XTickLabel', [])

subplot(6,1,6)
plot(t_out3/T, rad2deg(mod(x_out3(:,6), 2*pi)), 'g', ...
     t_out5/T, rad2deg(mod(x_out5(:,6), 2*pi)), 'r--')
ylabel('\nu (deg)')
xlabel('Orbital Periods')
grid on


% FIGURE 15: Comparison of all 3 Methods
figure(15)
subplot(6,1,1)
plot(t_out3/T, x_out3(:,1), 'g', ...
     t_out4/T, x_out4(:,1), 'b', ...
     t_out5/T, x_out5(:,1), 'r--')
ylabel('a (m)')
legend('Unperturbed','J2 Osc','J2 Mean')
grid on
set(gca, 'XTickLabel', [])

subplot(6,1,2)
plot(t_out3/T, x_out3(:,2), 'g', ...
     t_out4/T, x_out4(:,2), 'b', ...
     t_out5/T, x_out5(:,2), 'r--')
ylabel('e')
grid on
set(gca, 'XTickLabel', [])

subplot(6,1,3)
plot(t_out3/T, rad2deg(x_out3(:,3)), 'g', ...
     t_out4/T, rad2deg(x_out4(:,3)), 'b', ...
     t_out5/T, rad2deg(x_out5(:,3)), 'r--')
ylabel('i (deg)')
grid on
set(gca, 'XTickLabel', [])

subplot(6,1,4)
plot(t_out3/T, rad2deg(mod(x_out3(:,4), 2*pi)), 'g', ...
     t_out4/T, rad2deg(mod(x_out4(:,4), 2*pi)), 'b', ...
     t_out5/T, rad2deg(mod(x_out5(:,4), 2*pi)), 'r--')
ylabel('\Omega (deg)')
grid on
set(gca, 'XTickLabel', [])

subplot(6,1,5)
plot(t_out3/T, rad2deg(mod(x_out3(:,5), 2*pi)), 'g', ...
     t_out4/T, rad2deg(mod(x_out4(:,5), 2*pi)), 'b', ...
     t_out5/T, rad2deg(mod(x_out5(:,5), 2*pi)), 'r--')
ylabel('\omega (deg)')
grid on
set(gca, 'XTickLabel', [])

subplot(6,1,6)
plot(t_out3/T, rad2deg(mod(x_out3(:,6), 2*pi)), 'g', ...
     t_out4/T, rad2deg(mod(x_out4(:,6), 2*pi)), 'b', ...
     t_out5/T, rad2deg(mod(x_out5(:,6), 2*pi)), 'r--')
ylabel('\nu (deg)')
xlabel('Orbital Periods')
grid on

