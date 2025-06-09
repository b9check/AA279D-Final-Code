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
tspan = 0:10:10*T;  
state0 = chief_ECI;

% Call ode113 without J2 perturbations
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t, state, mu, R, J2, perturbated);
[t_out1, x_out1] = ode113(odefun, tspan, state0);

% Call ode113 with J2 perturbations
perturbated = true;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t, state, mu, R, J2, perturbated);
[t_out2, x_out2] = ode113(odefun, tspan, state0);

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


%% 2d

