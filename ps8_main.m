%% - Reset

clc; clear; close all;

%% Initialization 

% Orbited body
body = 'earth';
const = utils.getConstants({body});

% Chief initial orbital elements
a_0 = 6771000 ;         % semi-major axis [m]
e_0 = 0.0005;           % eccentricity [-]
i_0 = deg2rad(51.64);   % inclination [rad]
W_0 = deg2rad(257);     % RAAN [rad]
w_0 = deg2rad(45);      % argument of perigee [rad]
f_0 = deg2rad(30);      % true anomaly [rad]

% Chief initial state in OE & ECI
initial_state_0_OE = [a_0, e_0, i_0, W_0, w_0, f_0];
initial_state_0_ECI = utils.OE2ECI(initial_state_0_OE, const, body);
r0_init = initial_state_0_ECI(1:3);
v0_init = initial_state_0_ECI(4:6);

% Deputy initial relative quasi-non-singular OE
% initial_state_1_rQNS_OE = [0; deg2rad(-0.05); 0; 0; 0; 0];
%initial_state_1_rQNS_OE = [0; deg2rad(-0.5); 0; 0; 0; 0];
delta_a = 0;
delta_lambda = deg2rad(-3.3773);
delta_e_x = 0.0007;
delta_e_y = 0.0007;
delta_i_x = deg2rad(0.4989);
delta_i_y = deg2rad(0.7850);
initial_state_1_rQNS_OE = [delta_a; delta_lambda; delta_e_x; delta_e_y; delta_i_x; delta_i_y];

% Deputy initial state in OE & ECI
initial_state_1_OE = rQNSOE2OE(initial_state_0_OE,initial_state_1_rQNS_OE);
initial_state_1_ECI = utils.OE2ECI(initial_state_1_OE, const, body);
r1_init = initial_state_1_ECI(1:3);
v1_init = initial_state_1_ECI(4:6);


% Compute the distance between the two spacecraft
distance = norm(r1_init - r0_init);

% Display the distance
fprintf('Distance between the two spacecraft: %.2f meters\n', distance);


%% - Generate Ground Truth State

% Orbital period
mu = const.(body).mu;
T = 2 * pi * sqrt(a_0^3 / mu); % Orbital period [s]

% Simulation parameters
tstart = 0;
tint = 1;
tend = 15*T;
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Run the perturbed FODE propagator for chief in ECI
perturbated = true;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t, state_ECI_0_p] = ode113(odefun, (tstart:tint:tend)', initial_state_0_ECI, options);


% Run the perturbed FODE for deputy in ECI
perturbated = true;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t, state_ECI_1_p] = ode113(odefun, (tstart:tint:tend)', initial_state_1_ECI, options);


% Consolidated state
x_c_gt = state_ECI_0_p;
x_d_gt = ECI2rQNSOE(state_ECI_0_p, state_ECI_1_p, const, body);
x_d_gt(:,2) = unwrap(x_d_gt(:,2));
x_gt = [x_c_gt, x_d_gt];

%% - Generate Measurements

% Number of time steps
N = length(t);

% Initialize measurement array
z_k = zeros(N, 9); 

% Measurement noise covariance (diagonal, SI units)
sigma_pos = 0.001; % m, GNSS position (1 mm)
sigma_vel = 0.001; % m/s, GNSS velocity (1 mm/s)
sigma_range = 0.001; % m, laser range (1 mm)
sigma_angles = 4.848e-3; % rad, camera angles (~1000 arcsecond)

R_k = diag([sigma_pos^2, sigma_pos^2, sigma_pos^2, ...
            sigma_vel^2, sigma_vel^2, sigma_vel^2, ...
            sigma_range^2, sigma_angles^2, sigma_angles^2]);

% Square root of covariance
sqrt_R_k = sqrtm(R_k);

for k = 1:N
    % Ground truth state
    x_k_gt = x_gt(k,:)'; % 12 x 1
    
    % Compute noiseless measurement
    z_k_true = h_k(x_k_gt, zeros(9,1), const, body);
    
    % Generate noise
    v_k = sqrt_R_k * randn(9,1);
    
    % Add noise
    z_k(k,:) = z_k_true + v_k;
end


%% - UKF

% UKF Parameters
n = 12; % State dimension
m = 9;  % Measurement dimension
alpha = 6*1e-1;
kappa = 0;
beta = 2;
lambda = alpha^2 * (n + kappa) - n;

% Weights
Wm = zeros(2*n+1,1);
Wc = zeros(2*n+1,1);
Wm(1) = lambda / (n + lambda);
Wc(1) = lambda / (n + lambda) + (1 - alpha^2 + beta);
for i = 2:(2*n+1)
    Wm(i) = 1 / (2 * (n + lambda));
    Wc(i) = 1 / (2 * (n + lambda));
end

% Process noise covariance
sigma_pos_ECI = 1e-2; % m, position noise for chief ECI
sigma_vel_ECI = 1e-3; % m/s, velocity noise for chief ECI
sigma_rQNS_dimless = 1e-6; % dimensionless, for delta_a, delta_ex, delta_ey
sigma_rQNS_angle = 1e-6; % rad, for delta_lambda, delta_ix, delta_iy
Q_k = diag([sigma_pos_ECI^2 * ones(3,1); ... % ECI position
            sigma_vel_ECI^2 * ones(3,1); ... % ECI velocity
            sigma_rQNS_dimless^2; ... % delta_a
            sigma_rQNS_angle^2; ... % delta_lambda
            sigma_rQNS_dimless^2 * ones(2,1); ... % delta_ex, delta_ey
            sigma_rQNS_angle^2 * ones(2,1)]); % delta_ix, delta_iy

% Measurement noise covariance (diagonal, SI units)
sigma_pos = 0.001; % m, GNSS position (1 mm)
sigma_vel = 0.001; % m/s, GNSS velocity (1 mm/s)
sigma_range = 0.001; % m, laser range (1 mm)
sigma_angles = 4.848e-3; % rad, camera angles (~1000 arcsecond)

R_k = diag([sigma_pos^2, sigma_pos^2, sigma_pos^2, ...
            sigma_vel^2, sigma_vel^2, sigma_vel^2, ...
            sigma_range^2, sigma_angles^2, sigma_angles^2]);

% Initial state and covariance
x_hat = [initial_state_0_ECI; initial_state_1_rQNS_OE];
P_c = 1e-1 * eye(6);
P_d = 1e-2 * eye(6);
P = [P_c, zeros(6);
    zeros(6), P_d];

% Storage
x_hat_store = zeros(N, n);
P_store = zeros(n, n, N);

z_hat_store = zeros(N, m);
z_post_store = zeros(N, m);

for k = 1:N
    % Time step
    delta_t = tint;
    if k == 1
        delta_t = 0;
    end

    % Generate sigma points
    sqrt_P = chol((n + lambda) * P, 'lower');
    chi = zeros(n, 2*n+1);
    chi(:,1) = x_hat;
    for i = 1:n
        chi(:,i+1) = x_hat + sqrt_P(:,i);
        chi(:,i+n+1) = x_hat - sqrt_P(:,i);
    end

    % Prediction step
    chi_pred = zeros(n, 2*n+1);
    for i = 1:(2*n+1)
        chi_pred(:,i) = f_k(chi(:,i), zeros(6,1), zeros(n,1), delta_t, const, body);
    end
    x_hat_pred = zeros(n,1);
    for i = 1:(2*n+1)
        x_hat_pred = x_hat_pred + Wm(i) * chi_pred(:,i);
    end
    P_pred = Q_k;
    for i = 1:(2*n+1)
        diff = chi_pred(:,i) - x_hat_pred;
        P_pred = P_pred + Wc(i) * (diff * diff');
    end

    % Update step
    sqrt_P_pred = chol((n + lambda) * P_pred, 'lower');
    chi_pred_new = zeros(n, 2*n+1);
    chi_pred_new(:,1) = x_hat_pred;
    for i = 1:n
        chi_pred_new(:,i+1) = x_hat_pred + sqrt_P_pred(:,i);
        chi_pred_new(:,i+n+1) = x_hat_pred - sqrt_P_pred(:,i);
    end

    % Measurement sigma points
    Z_pred = zeros(m, 2*n+1);
    for i = 1:(2*n+1)
        Z_pred(:,i) = h_k(chi_pred_new(:,i), zeros(m,1), const, body);
    end
    z_hat = zeros(m,1);
    for i = 1:(2*n+1)
        z_hat = z_hat + Wm(i) * Z_pred(:,i);
    end

   
    % Measurement covariance and cross-covariance
    S_k = R_k;
    P_xz = zeros(n,m);
    for i = 1:(2*n+1)
        z_diff = Z_pred(:,i) - z_hat;
        S_k = S_k + Wc(i) * (z_diff * z_diff');
        x_diff = chi_pred_new(:,i) - x_hat_pred;
        P_xz = P_xz + Wc(i) * (x_diff * z_diff');
    end

    % Kalman gain
    K_k = P_xz / S_k;

    % Update state and covariance
    z_k_current = z_k(k,:)';
    x_hat = x_hat_pred + K_k * (z_k_current - z_hat);
    P = P_pred - K_k * S_k * K_k';

    % Store results
    x_hat_store(k,:) = x_hat';
    P_store(:,:,k) = P;
    z_hat_store(k,:) = z_hat';
    z_post_store(k,:) = h_k(x_hat, zeros(m,1), const, body)';
end

%% - Plots

% Plot true vs predicted chief ECI position
figure;
subplot(3,1,1);
plot(t, x_c_gt(:,1), 'b-', 'DisplayName', 'True X');
hold on;
plot(t, x_hat_store(:,1), 'r--', 'DisplayName', 'Predicted X');
%title('Chief ECI Position - X');
xlabel('Time (s)');
ylabel('Position (m)');
legend;
grid on;

subplot(3,1,2);
plot(t, x_c_gt(:,2), 'b-', 'DisplayName', 'True Y');
hold on;
plot(t, x_hat_store(:,2), 'r--', 'DisplayName', 'Predicted Y');
%title('Chief ECI Position - Y');
xlabel('Time (s)');
ylabel('Position (m)');
legend;
grid on;

subplot(3,1,3);
plot(t, x_c_gt(:,3), 'b-', 'DisplayName', 'True Z');
hold on;
plot(t, x_hat_store(:,3), 'r--', 'DisplayName', 'Predicted Z');
%title('Chief ECI Position - Z');
xlabel('Time (s)');
ylabel('Position (m)');
legend;
grid on;


% Plot true vs predicted chief ECI velocity
figure;
subplot(3,1,1);
plot(t, x_c_gt(:,4),  'b-',  'DisplayName','True V_x');
hold on;
plot(t, x_hat_store(:,4), 'r--','DisplayName','Predicted V_x');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend;
grid on;

subplot(3,1,2);
plot(t, x_c_gt(:,5),  'b-',  'DisplayName','True V_y');
hold on;
plot(t, x_hat_store(:,5), 'r--','DisplayName','Predicted V_y');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend;
grid on;

subplot(3,1,3);
plot(t, x_c_gt(:,6),  'b-',  'DisplayName','True V_z');
hold on;
plot(t, x_hat_store(:,6), 'r--','DisplayName','Predicted V_z');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend;
grid on;


% Plot true vs predicted relative QNS OE with errors
% Plot 1: delta_a
figure;
subplot(2,1,1);
plot(t, x_d_gt(:,1), 'b-', 'DisplayName', 'True \delta a');
hold on;
plot(t, x_hat_store(:,7), 'r--', 'DisplayName', 'Predicted \delta a');
%title('Relative QNS OE - \delta a');
xlabel('Time (s)');
ylabel('\delta a');
legend;
grid on;

subplot(2,1,2);
plot(t, x_d_gt(:,1) - x_hat_store(:,7), 'k-', 'DisplayName', 'Error \delta a');
%title('Error in \delta a');
xlabel('Time (s)');
ylabel('Error');
legend;
grid on;

% Plot 2: delta_lambda
figure;
subplot(2,1,1);
plot(t, x_d_gt(:,2), 'b-', 'DisplayName', 'True \delta \lambda');
hold on;
plot(t, x_hat_store(:,8), 'r--', 'DisplayName', 'Predicted \delta \lambda');
title('Relative QNS OE - \delta \lambda');
xlabel('Time (s)');
ylabel('\delta \lambda (rad)');
legend;
grid on;

subplot(2,1,2);
plot(t, x_d_gt(:,2) - x_hat_store(:,8), 'k-', 'DisplayName', 'Error \delta \lambda');
%title('Error in \delta \lambda');
xlabel('Time (s)');
ylabel('Error (rad)');
legend;
grid on;

% Plot 3: delta_e_x
figure;
subplot(2,1,1);
plot(t, x_d_gt(:,3), 'b-', 'DisplayName', 'True \delta e_x');
hold on;
plot(t, x_hat_store(:,9), 'r--', 'DisplayName', 'Predicted \delta e_x');
%title('Relative QNS OE - \delta e_x');
xlabel('Time (s)');
ylabel('\delta e_x');
legend;
grid on;

subplot(2,1,2);
plot(t, x_d_gt(:,3) - x_hat_store(:,9), 'k-', 'DisplayName', 'Error \delta e_x');
%title('Error in \delta e_x');
xlabel('Time (s)');
ylabel('Error');
legend;
grid on;

% Plot 4: delta_e_y
figure;
subplot(2,1,1);
plot(t, x_d_gt(:,4), 'b-', 'DisplayName', 'True \delta e_y');
hold on;
plot(t, x_hat_store(:,10), 'r--', 'DisplayName', 'Predicted \delta e_y');
%title('Relative QNS OE - \delta e_y');
xlabel('Time (s)');
ylabel('\delta e_y');
legend;
grid on;

subplot(2,1,2);
plot(t, x_d_gt(:,4) - x_hat_store(:,10), 'k-', 'DisplayName', 'Error \delta e_y');
%title('Error in \delta e_y');
xlabel('Time (s)');
ylabel('Error');
legend;
grid on;

% Plot 5: delta_i_x
figure;
subplot(2,1,1);
plot(t, x_d_gt(:,5), 'b-', 'DisplayName', 'True \delta i_x');
hold on;
plot(t, x_hat_store(:,11), 'r--', 'DisplayName', 'Predicted \delta i_x');
title('Relative QNS OE - \delta i_x');
xlabel('Time (s)');
ylabel('\delta i_x (rad)');
legend;
grid on;

subplot(2,1,2);
plot(t, x_d_gt(:,5) - x_hat_store(:,11), 'k-', 'DisplayName', 'Error \delta i_x');
%title('Error in \delta i_x');
xlabel('Time (s)');
ylabel('Error (rad)');
legend;
grid on;

% Plot 6: delta_i_y
figure;
subplot(2,1,1);
plot(t, x_d_gt(:,6), 'b-', 'DisplayName', 'True \delta i_y');
hold on;
plot(t, x_hat_store(:,12), 'r--', 'DisplayName', 'Predicted \delta i_y');
%title('Relative QNS OE - \delta i_y');
xlabel('Time (s)');
ylabel('\delta i_y (rad)');
legend;
grid on;

subplot(2,1,2);
plot(t, x_d_gt(:,6) - x_hat_store(:,12), 'k-', 'DisplayName', 'Error \delta i_y');
%title('Error in \delta i_y');
xlabel('Time (s)');
ylabel('Error (rad)');
legend;
grid on;

%% - Plots V2

% Plot true vs predicted chief ECI position
figure;
subplot(3,1,1);
plot(t, x_c_gt(:,1), 'b-', 'DisplayName', 'True X');
hold on;
plot(t, x_hat_store(:,1), 'r--', 'DisplayName', 'Predicted X');
xlabel('Time (s)');
ylabel('Position (m)');
legend;
grid on;

subplot(3,1,2);
plot(t, x_c_gt(:,2), 'b-', 'DisplayName', 'True Y');
hold on;
plot(t, x_hat_store(:,2), 'r--', 'DisplayName', 'Predicted Y');
xlabel('Time (s)');
ylabel('Position (m)');
legend;
grid on;

subplot(3,1,3);
plot(t, x_c_gt(:,3), 'b-', 'DisplayName', 'True Z');
hold on;
plot(t, x_hat_store(:,3), 'r--', 'DisplayName', 'Predicted Z');
xlabel('Time (s)');
ylabel('Position (m)');
legend;
grid on;

% Plot chief ECI position error with 3-sigma bounds
figure;
subplot(3,1,1);
plot(t, x_c_gt(:,1) - x_hat_store(:,1), 'k-', 'DisplayName', 'Error X');
hold on;
plot(t, 3*sqrt(squeeze(P_store(1,1,:))), 'g--', 'DisplayName', '+3\sigma');
plot(t, -3*sqrt(squeeze(P_store(1,1,:))), 'g--', 'DisplayName', '-3\sigma');
xlabel('Time (s)');
ylabel('Error (m)');
legend;
grid on;

subplot(3,1,2);
plot(t, x_c_gt(:,2) - x_hat_store(:,2), 'k-', 'DisplayName', 'Error Y');
hold on;
plot(t, 3*sqrt(squeeze(P_store(2,2,:))), 'g--', 'DisplayName', '+3\sigma');
plot(t, -3*sqrt(squeeze(P_store(2,2,:))), 'g--', 'DisplayName', '-3\sigma');
xlabel('Time (s)');
ylabel('Error (m)');
legend;
grid on;

subplot(3,1,3);
plot(t, x_c_gt(:,3) - x_hat_store(:,3), 'k-', 'DisplayName', 'Error Z');
hold on;
plot(t, 3*sqrt(squeeze(P_store(3,3,:))), 'g--', 'DisplayName', '+3\sigma');
plot(t, -3*sqrt(squeeze(P_store(3,3,:))), 'g--', 'DisplayName', '-3\sigma');
xlabel('Time (s)');
ylabel('Error (m)');
legend;
grid on;

% Plot true vs predicted chief ECI velocity
figure;
subplot(3,1,1);
plot(t, x_c_gt(:,4), 'b-', 'DisplayName', 'True V_x');
hold on;
plot(t, x_hat_store(:,4), 'r--', 'DisplayName', 'Predicted V_x');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend;
grid on;

subplot(3,1,2);
plot(t, x_c_gt(:,5), 'b-', 'DisplayName', 'True V_y');
hold on;
plot(t, x_hat_store(:,5), 'r--', 'DisplayName', 'Predicted V_y');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend;
grid on;

subplot(3,1,3);
plot(t, x_c_gt(:,6), 'b-', 'DisplayName', 'True V_z');
hold on;
plot(t, x_hat_store(:,6), 'r--', 'DisplayName', 'Predicted V_z');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend;
grid on;

% Plot chief ECI velocity error with 3-sigma bounds
figure;
subplot(3,1,1);
plot(t, x_c_gt(:,4) - x_hat_store(:,4), 'k-', 'DisplayName', 'Error V_x');
hold on;
plot(t, 3*sqrt(squeeze(P_store(4,4,:))), 'g--', 'DisplayName', '+3\sigma');
plot(t, -3*sqrt(squeeze(P_store(4,4,:))), 'g--', 'DisplayName', '-3\sigma');
xlabel('Time (s)');
ylabel('Error (m/s)');
legend;
grid on;

subplot(3,1,2);
plot(t, x_c_gt(:,5) - x_hat_store(:,5), 'k-', 'DisplayName', 'Error V_y');
hold on;
plot(t, 3*sqrt(squeeze(P_store(5,5,:))), 'g--', 'DisplayName', '+3\sigma');
plot(t, -3*sqrt(squeeze(P_store(5,5,:))), 'g--', 'DisplayName', '-3\sigma');
xlabel('Time (s)');
ylabel('Error (m/s)');
legend;
grid on;

subplot(3,1,3);
plot(t, x_c_gt(:,6) - x_hat_store(:,6), 'k-', 'DisplayName', 'Error V_z');
hold on;
plot(t, 3*sqrt(squeeze(P_store(6,6,:))), 'g--', 'DisplayName', '+3\sigma');
plot(t, -3*sqrt(squeeze(P_store(6,6,:))), 'g--', 'DisplayName', '-3\sigma');
xlabel('Time (s)');
ylabel('Error (m/s)');
legend;
grid on;

% Plot true vs predicted relative QNS OE with errors
% Plot 1: delta_a
figure;
subplot(2,1,1);
plot(t, x_d_gt(:,1), 'b-', 'DisplayName', 'True \delta a');
hold on;
plot(t, x_hat_store(:,7), 'r--', 'DisplayName', 'Predicted \delta a');
xlabel('Time (s)');
ylabel('\delta a');
legend;
grid on;

subplot(2,1,2);
plot(t, x_d_gt(:,1) - x_hat_store(:,7), 'k-', 'DisplayName', 'Error \delta a');
hold on;
plot(t, 3*sqrt(squeeze(P_store(7,7,:))), 'g--', 'DisplayName', '+3\sigma');
plot(t, -3*sqrt(squeeze(P_store(7,7,:))), 'g--', 'DisplayName', '-3\sigma');
xlabel('Time (s)');
ylabel('Error');
legend;
grid on;

% Plot 2: delta_lambda
figure;
subplot(2,1,1);
plot(t, x_d_gt(:,2), 'b-', 'DisplayName', 'True \delta \lambda');
hold on;
plot(t, x_hat_store(:,8), 'r--', 'DisplayName', 'Predicted \delta \lambda');
xlabel('Time (s)');
ylabel('\delta \lambda (rad)');
legend;
grid on;

subplot(2,1,2);
plot(t, x_d_gt(:,2) - x_hat_store(:,8), 'k-', 'DisplayName', 'Error \delta \lambda');
hold on;
plot(t, 3*sqrt(squeeze(P_store(8,8,:))), 'g--', 'DisplayName', '+3\sigma');
plot(t, -3*sqrt(squeeze(P_store(8,8,:))), 'g--', 'DisplayName', '-3\sigma');
xlabel('Time (s)');
ylabel('Error (rad)');
legend;
grid on;

% Plot 3: delta_e_x
figure;
subplot(2,1,1);
plot(t, x_d_gt(:,3), 'b-', 'DisplayName', 'True \delta e_x');
hold on;
plot(t, x_hat_store(:,9), 'r--', 'DisplayName', 'Predicted \delta e_x');
xlabel('Time (s)');
ylabel('\delta e_x');
legend;
grid on;

subplot(2,1,2);
plot(t, x_d_gt(:,3) - x_hat_store(:,9), 'k-', 'DisplayName', 'Error \delta e_x');
hold on;
plot(t, 3*sqrt(squeeze(P_store(9,9,:))), 'g--', 'DisplayName', '+3\sigma');
plot(t, -3*sqrt(squeeze(P_store(9,9,:))), 'g--', 'DisplayName', '-3\sigma');
xlabel('Time (s)');
ylabel('Error');
legend;
grid on;

% Plot 4: delta_e_y
figure;
subplot(2,1,1);
plot(t, x_d_gt(:,4), 'b-', 'DisplayName', 'True \delta e_y');
hold on;
plot(t, x_hat_store(:,10), 'r--', 'DisplayName', 'Predicted \delta e_y');
xlabel('Time (s)');
ylabel('\delta e_y');
legend;
grid on;

subplot(2,1,2);
plot(t, x_d_gt(:,4) - x_hat_store(:,10), 'k-', 'DisplayName', 'Error \delta e_y');
hold on;
plot(t, 3*sqrt(squeeze(P_store(10,10,:))), 'g--', 'DisplayName', '+3\sigma');
plot(t, -3*sqrt(squeeze(P_store(10,10,:))), 'g--', 'DisplayName', '-3\sigma');
xlabel('Time (s)');
ylabel('Error');
legend;
grid on;

% Plot 5: delta_i_x
figure;
subplot(2,1,1);
plot(t, x_d_gt(:,5), 'b-', 'DisplayName', 'True \delta i_x');
hold on;
plot(t, x_hat_store(:,11), 'r--', 'DisplayName', 'Predicted \delta i_x');
xlabel('Time (s)');
ylabel('\delta i_x (rad)');
legend;
grid on;

subplot(2,1,2);
plot(t, x_d_gt(:,5) - x_hat_store(:,11), 'k-', 'DisplayName', 'Error \delta i_x');
hold on;
plot(t, 3*sqrt(squeeze(P_store(11,11,:))), 'g--', 'DisplayName', '+3\sigma');
plot(t, -3*sqrt(squeeze(P_store(11,11,:))), 'g--', 'DisplayName', '-3\sigma');
xlabel('Time (s)');
ylabel('Error (rad)');
legend;
grid on;

% Plot 6: delta_i_y
figure;
subplot(2,1,1);
plot(t, x_d_gt(:,6), 'b-', 'DisplayName', 'True \delta i_y');
hold on;
plot(t, x_hat_store(:,12), 'r--', 'DisplayName', 'Predicted \delta i_y');
xlabel('Time (s)');
ylabel('\delta i_y (rad)');
legend;
grid on;

subplot(2,1,2);
plot(t, x_d_gt(:,6) - x_hat_store(:,12), 'k-', 'DisplayName', 'Error \delta i_y');
hold on;
plot(t, 3*sqrt(squeeze(P_store(12,12,:))), 'g--', 'DisplayName', '+3\sigma');
plot(t, -3*sqrt(squeeze(P_store(12,12,:))), 'g--', 'DisplayName', '-3\sigma');
xlabel('Time (s)');
ylabel('Error (rad)');
legend;
grid on;


%% - Additional Plots

% True statistics (mean and std of state error in last orbit)
last_orbit_idx = t >= (14*T); % Last orbit (14T to 15T)
state_errors = x_gt - x_hat_store;
state_errors_last = state_errors(last_orbit_idx, :);
mean_errors = mean(state_errors_last);
std_errors = std(state_errors_last);

% Table of statistics
figure;
uitable('Data', [mean_errors; std_errors], ...
        'ColumnName', {'r_x', 'r_y', 'r_z', 'v_x', 'v_y', 'v_z', ...
                       '\delta a', '\delta \lambda', '\delta e_x', '\delta e_y', '\delta i_x', '\delta i_y'}, ...
        'RowName', {'Mean', 'Std'}, ...
        'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.8]);
title('True Statistics of State Errors (Last Orbit)');

% Pre-fit and post-fit residuals
pre_fit_residuals = z_k - z_hat_store;
post_fit_residuals = z_k - z_post_store;

% Plot residuals for r_c (position)
figure;
subplot(3,1,1);
plot(t, pre_fit_residuals(:,1), 'b-', 'DisplayName', 'Pre-fit Residual r_x');
hold on;
plot(t, post_fit_residuals(:,1), 'r--', 'DisplayName', 'Post-fit Residual r_x');
plot(t, sigma_pos*ones(size(t)), 'k:', 'DisplayName', '+1\sigma_{pos}');
plot(t, -sigma_pos*ones(size(t)), 'k-.', 'DisplayName', '-1\sigma_{pos}');
xlabel('Time (s)');
ylabel('Residual r_x (m)');
legend('Location', 'best');
grid on;

subplot(3,1,2);
plot(t, pre_fit_residuals(:,2), 'b-', 'DisplayName', 'Pre-fit Residual r_y');
hold on;
plot(t, post_fit_residuals(:,2), 'r--', 'DisplayName', 'Post-fit Residual r_y');
plot(t, sigma_pos*ones(size(t)), 'k:', 'DisplayName', '+1\sigma_{pos}');
plot(t, -sigma_pos*ones(size(t)), 'k-.', 'DisplayName', '-1\sigma_{pos}');
xlabel('Time (s)');
ylabel('Residual r_y (m)');
legend('Location', 'best');
grid on;

subplot(3,1,3);
plot(t, pre_fit_residuals(:,3), 'b-', 'DisplayName', 'Pre-fit Residual r_z');
hold on;
plot(t, post_fit_residuals(:,3), 'r--', 'DisplayName', 'Post-fit Residual r_z');
plot(t, sigma_pos*ones(size(t)), 'k:', 'DisplayName', '+1\sigma_{pos}');
plot(t, -sigma_pos*ones(size(t)), 'k-.', 'DisplayName', '-1\sigma_{pos}');
xlabel('Time (s)');
ylabel('Residual r_z (m)');
legend('Location', 'best');
grid on;

% Plot residuals for v_c (velocity)
figure;
subplot(3,1,1);
plot(t, pre_fit_residuals(:,4), 'b-', 'DisplayName', 'Pre-fit Residual v_x');
hold on;
plot(t, post_fit_residuals(:,4), 'r--', 'DisplayName', 'Post-fit Residual v_x');
plot(t, sigma_vel*ones(size(t)), 'k:', 'DisplayName', '+1\sigma_{vel}');
plot(t, -sigma_vel*ones(size(t)), 'k-.', 'DisplayName', '-1\sigma_{vel}');
xlabel('Time (s)');
ylabel('Residual v_x (m/s)');
legend('Location', 'best');
grid on;

subplot(3,1,2);
plot(t, pre_fit_residuals(:,5), 'b-', 'DisplayName', 'Pre-fit Residual v_y');
hold on;
plot(t, post_fit_residuals(:,5), 'r--', 'DisplayName', 'Post-fit Residual v_y');
plot(t, sigma_vel*ones(size(t)), 'k:', 'DisplayName', '+1\sigma_{vel}');
plot(t, -sigma_vel*ones(size(t)), 'k-.', 'DisplayName', '-1\sigma_{vel}');
xlabel('Time (s)');
ylabel('Residual v_y (m/s)');
legend('Location', 'best');
grid on;

subplot(3,1,3);
plot(t, pre_fit_residuals(:,6), 'b-', 'DisplayName', 'Pre-fit Residual v_z');
hold on;
plot(t, post_fit_residuals(:,6), 'r--', 'DisplayName', 'Post-fit Residual v_z');
plot(t, sigma_vel*ones(size(t)), 'k:', 'DisplayName', '+1\sigma_{vel}');
plot(t, -sigma_vel*ones(size(t)), 'k-.', 'DisplayName', '-1\sigma_{vel}');
xlabel('Time (s)');
ylabel('Residual v_z (m/s)');
legend('Location', 'best');
grid on;

% Plot residuals for range, azimuth, elevation
figure;
subplot(3,1,1);
plot(t, pre_fit_residuals(:,7), 'b-', 'DisplayName', 'Pre-fit Residual Range');
hold on;
plot(t, post_fit_residuals(:,7), 'r--', 'DisplayName', 'Post-fit Residual Range');
plot(t, sigma_range*ones(size(t)), 'k:', 'DisplayName', '+1\sigma_{range}');
plot(t, -sigma_range*ones(size(t)), 'k-.', 'DisplayName', '-1\sigma_{range}');
xlabel('Time (s)');
ylabel('Residual Range (m)');
legend('Location', 'best');
grid on;

subplot(3,1,2);
plot(t, pre_fit_residuals(:,8), 'b-', 'DisplayName', 'Pre-fit Residual Azimuth');
hold on;
plot(t, post_fit_residuals(:,8), 'r--', 'DisplayName', 'Post-fit Residual Azimuth');
plot(t, sigma_angles*ones(size(t)), 'k:', 'DisplayName', '+1\sigma_{angles}');
plot(t, -sigma_angles*ones(size(t)), 'k-.', 'DisplayName', '-1\sigma_{angles}');
xlabel('Time (s)');
ylabel('Residual Azimuth (rad)');
legend('Location', 'best');
grid on;

subplot(3,1,3);
plot(t, pre_fit_residuals(:,9), 'b-', 'DisplayName', 'Pre-fit Residual Elevation');
hold on;
plot(t, post_fit_residuals(:,9), 'r--', 'DisplayName', 'Post-fit Residual Elevation');
plot(t, sigma_angles*ones(size(t)), 'k:', 'DisplayName', '+1\sigma_{angles}');
plot(t, -sigma_angles*ones(size(t)), 'k-.', 'DisplayName', '-1\sigma_{angles}');
xlabel('Time (s)');
ylabel('Residual Elevation (rad)');
legend('Location', 'best');
grid on;


%% - Functions

function x_k = f_k(x_k1, u_k, w_k, delta_t, const, body)
    % Constants
    mu = const.(body).mu;
    J2 = const.(body).J2;
    R = const.(body).R;

    % Chief dynamics (ECI)
    x_c_k1 = x_k1(1:6);
    u_k_c = u_k(1:3);
    w_k_c = w_k(1:6);
    r_c_k1_n = norm(x_k1(1:3));
    x_c_k1_x = x_k1(1);
    x_c_k1_y = x_k1(2);
    x_c_k1_z = x_k1(3);
    help_vec = [x_c_k1_x*(1-5*x_c_k1_z^2/r_c_k1_n^2);
                x_c_k1_y*(1-5*x_c_k1_z^2/r_c_k1_n^2);
                x_c_k1_z*(3-5*x_c_k1_z^2/r_c_k1_n^2)];
    a_J2 = -(3*mu*J2*R^2/(2*r_c_k1_n^5))*help_vec;
    if delta_t == 0
        a_control = zeros(3,1);
    else
        a_control = u_k_c / delta_t;
    end
    two_body_acc = -(mu/r_c_k1_n^3)*x_k1(1:3);
    total_acc = two_body_acc + a_J2 + a_control;
    
    x_c_k1_dot = zeros(6,1);
    x_c_k1_dot(1:3) = x_k1(4:6);
    x_c_k1_dot(4:6) = total_acc;
    
    f_k_c = x_c_k1 + x_c_k1_dot*delta_t + w_k_c;

    % Deputy dynamics (rQNSOE)
    x_d_k1 = x_k1(7:12);
    u_k_d = u_k(4:6);
    w_k_d = w_k(7:12);
    chief_OE_k1 = utils.ECI2OE(x_k1(1:6), const, body);
    Phi_J2_qns_k1 = getQNS_J2_STM(chief_OE_k1, delta_t, const, body);
    B_k1 = getBcMatrix(chief_OE_k1, const, body);
    f_k_d = Phi_J2_qns_k1*x_d_k1 + delta_t*B_k1*u_k_d + w_k_d;

    x_k = [f_k_c; f_k_d];
end


function z_k = h_k(x_k, v_k, const, body)
    x_c_k = x_k(1:6);
    r_c_k = x_k(1:3);
    v_c_k = x_k(4:6);
    rQNSOE_d_k = x_k(7:12);
    OE_c_k = utils.ECI2OE(x_c_k,const,body);
    OE_d_k = rQNSOE2OE(OE_c_k,rQNSOE_d_k);
    x_d_k = utils.OE2ECI(OE_d_k,const,body);
    r_d_k = x_d_k(1:3);
    % v_d_k = x_d_k(4:6);
    delta_r = r_d_k - r_c_k;
    delta_r_n = norm(delta_r);
    range = sqrt(delta_r' * delta_r);
    delta_r_x = delta_r(1);
    delta_r_y = delta_r(2);
    delta_r_z = delta_r(3);
    azimuth = atan2(delta_r_y, delta_r_x);
    elevation = asin(delta_r_z/delta_r_n);
    z_k = [r_c_k; v_c_k; range; azimuth; elevation] + v_k;
end


%%

function Phi_J2_qns = getQNS_J2_STM(OE, tau, const, body)
    % Extract orbital elements
    a = OE(1); e = OE(2); i = OE(3);
    W = OE(4); w = OE(5);

    % Constants
    mu = const.(body).mu;
    J2 = const.(body).J2;
    R = const.(body).R;

    % Derived quantities
    eta = sqrt(1 - e^2);
    kappa = (3/4) * J2 * R^2 * sqrt(mu) / (a^(7/2) * eta^4);

    % Auxiliary terms
    E = 1 + eta;
    F = 4 + 3 * eta;
    G = 1 / eta^2;
    P = 3 * cos(i)^2 - 1;
    Q = 5 * cos(i)^2 - 1;
    R = cos(i);
    S = sin(2*i);
    T = sin(i)^2;

    % Perifocal eccentricity vector components
    exi = e * cos(w);
    eyi = e * sin(w);

    % Secular rates
    w_dot = kappa * Q;
    
    % Final perigee argument
    omega_f = w + w_dot * tau;
    exf = e * cos(omega_f);
    eyf = e * sin(omega_f);

    % Trig terms
    cos_wt = cos(w_dot * tau);
    sin_wt = sin(w_dot * tau);

    % Construct STM
    Phi_J2_qns = zeros(6,6);
    Phi_J2_qns(1,1) = 1;
    Phi_J2_qns(2,1) = -((3/2)*sqrt(mu/a^3) + (7/2)*kappa*E*P) * tau;
    Phi_J2_qns(2,2) = 1;
    Phi_J2_qns(2,3) = kappa * exi * F * G * P * tau;
    Phi_J2_qns(2,4) = kappa * eyi * F * G * P * tau;
    Phi_J2_qns(2,5) = -kappa * F * S * tau;

    Phi_J2_qns(3,1) = (7/2) * kappa * eyf * Q * tau;
    Phi_J2_qns(3,3) = cos_wt - 4 * kappa * exi * eyf * G * Q * tau;
    Phi_J2_qns(3,4) = -sin_wt - 4 * kappa * eyi * eyf * G * Q * tau;
    Phi_J2_qns(3,5) = 5 * kappa * eyf * S * tau;

    Phi_J2_qns(4,1) = -(7/2) * kappa * exf * Q * tau;
    Phi_J2_qns(4,3) = sin_wt + 4 * kappa * exi * exf * G * Q * tau;
    Phi_J2_qns(4,4) = cos_wt + 4 * kappa * eyi * exf * G * Q * tau;
    Phi_J2_qns(4,5) = -5 * kappa * exf * S * tau;

    Phi_J2_qns(5,5) = 1;
    Phi_J2_qns(6,1) = (7/2) * kappa * S * tau;
    Phi_J2_qns(6,3) = -4 * kappa * exi * G * S * tau;
    Phi_J2_qns(6,4) = -4 * kappa * eyi * G * S * tau;
    Phi_J2_qns(6,5) = 2 * kappa * T * tau;
    Phi_J2_qns(6,6) = 1;
end


function Bc = getBcMatrix(chief_OE, const, body)
     % Extract orbital elements
    a = chief_OE(1); e = chief_OE(2); i = chief_OE(3);
    W = chief_OE(4); w = chief_OE(5); f = chief_OE(6);

    % Constants
    mu = const.(body).mu;
    J2 = const.(body).J2;
    R = const.(body).R;

    % Derived quantities
    eta = sqrt(1 - e^2);
    nc = sqrt(mu / a^3);
    gamma = (3/4) * J2 * R^2 * sqrt(mu);
    kappa = gamma / (a^(7/2) * eta^4);

    % Trig values
    sin_wf = sin(w + f);
    cos_wf = cos(w + f);
    tan_i = tan(i);
    cos_f = cos(f);
    sin_f = sin(f);

    % Eccentricity vector components
    ex = e * cos(w);
    ey = e * sin(w);

    % Common factors
    denom = 1 + e * cos_f;
    ec2pf = (2 + e * cos_f); % shorthand

    % Compute B matrix
    Bc = zeros(6,3);
    Bc(1,1) = 2 * e * sin_f / eta;
    Bc(1,2) = 2 * (1 + e * cos_f) / eta;
    Bc(2,1) = -2 * eta^2 / denom;

    Bc(3,1) = eta * sin_wf;
    Bc(3,2) = eta * (ec2pf * cos_wf + ex) / denom;
    Bc(3,3) = eta * ey * sin_wf / (tan_i * denom);

    Bc(4,1) = -eta * cos_wf;
    Bc(4,2) = eta * (ec2pf * sin_wf + ey) / denom;
    Bc(4,3) = -eta * ex * sin_wf / (tan_i * denom);

    Bc(5,3) = eta * cos_wf / denom;

    Bc(6,3) = eta * sin_wf / denom;

    % Normalize
    Bc = Bc / (a * nc);
end


function deputy_OE = rQNSOE2OE(chief_OE, deputy_rQNSOE)
    % Unpack chief OEs
    a_c = chief_OE(1);
    e_c = chief_OE(2);
    i_c = chief_OE(3);
    W_c = chief_OE(4);
    w_c = chief_OE(5);
    f_c = chief_OE(6);
    E_c = 2 * atan( sqrt((1 - e_c)/(1 + e_c)) * tan(f_c/2) );  % Eccentric anomaly
    M_c = E_c - e_c * sin(E_c);                               % Mean anomaly
    %lambda_c = W_c + w_c + M_c;

    % Unpack deputy rQNS OEs
    delta_a = deputy_rQNSOE(1);
    delta_lambda = deputy_rQNSOE(2);
    delta_e_x = deputy_rQNSOE(3);
    delta_e_y = deputy_rQNSOE(4);
    delta_i_x = deputy_rQNSOE(5);
    delta_i_y = deputy_rQNSOE(6);

    % Deputy orbital elements
    a_d = a_c*(delta_a + 1);
    W_d = delta_i_y*sin(i_c) + W_c;
    i_d = delta_i_x + i_c;
    alpha = delta_e_y + e_c*sin(w_c);
    beta = delta_e_x + e_c*cos(w_c);
    w_d = atan2(alpha, beta);
    e_d = sqrt(alpha^2 + beta^2);
    M_d = delta_lambda - (W_d - W_c)*cos(i_c) + (M_c + w_c) - w_d;
    E_d = newton_raphson(M_d,e_d);
    f_d = 2 * atan2(sqrt(1 + e_d)*sin(E_d/2), sqrt(1 - e_d) * cos(E_d / 2));
    deputy_OE = [a_d, e_d, i_d, W_d, w_d, f_d];
end


function E = newton_raphson(M, e, epsilon)
    if nargin < 3
        epsilon = 1e-10;
    end

    E = M;
    max_iter = 1e5;

    for i = 1:max_iter
        f_E = E - e * sin(E) - M;
        f_prime_E = 1 - e * cos(E);
        increment = f_E / f_prime_E;
        E = E - increment;

        if abs(increment) <= epsilon
            break;
        end
    end
end



function rQNSOE = OE2rQNSOE(oe_chief, oe_deputy)
    % Extract chief elements
    a_0 = oe_chief(1);
    e_0 = oe_chief(2);
    i_0 = oe_chief(3);
    W_0 = oe_chief(4);
    w_0 = oe_chief(5);
    f_0 = oe_chief(6);

    % Extract deputy elements
    a_1 = oe_deputy(1);
    e_1 = oe_deputy(2);
    i_1 = oe_deputy(3);
    W_1 = oe_deputy(4);
    w_1 = oe_deputy(5);
    f_1 = oe_deputy(6);

    % True anomaly → Eccentric anomaly
    E_0 = 2 * atan( sqrt((1 - e_0)/(1 + e_0)) * tan(f_0 / 2) );
    E_1 = 2 * atan( sqrt((1 - e_1)/(1 + e_1)) * tan(f_1 / 2) );

    % Eccentric anomaly → Mean anomaly
    M_0 = E_0 - e_0 * sin(E_0);
    M_1 = E_1 - e_1 * sin(E_1);

    % Compute QNS-ROE
    delta_a  = (a_1 - a_0) / a_0;
    delta_lambda = (M_1 + w_1) - (M_0 + w_0) + (W_1 - W_0) * cos(i_0);
    delta_ex = e_1 * cos(w_1) - e_0 * cos(w_0);
    delta_ey = e_1 * sin(w_1) - e_0 * sin(w_0);
    delta_ix = i_1 - i_0;
    delta_iy = (W_1 - W_0) * sin(i_0);

    rQNSOE = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy];
end

function rQNSOE = ECI2rQNSOE(chief_ECI, deputy_ECI, const, body)
    N = size(chief_ECI, 1);
    rQNSOE = zeros(N, 6);

    for k = 1:N
        % Extract ECI states
        chief_state = chief_ECI(k, :)';
        deputy_state = deputy_ECI(k, :)';

        % Convert to Keplerian OE
        oe_chief = utils.ECI2OE(chief_state, const, body);
        oe_deputy = utils.ECI2OE(deputy_state, const, body);

        % Compute relative QNS
        rQNSOE(k, :) = OE2rQNSOE(oe_chief, oe_deputy);
    end
end