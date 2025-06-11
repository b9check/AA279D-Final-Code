%% - Reset

clc; clear; close all;

%% Initialization 

% Constants
mu = 3.986004418e14;  
R = 6378137;          
J2 = 1.08262668e-3; 

% Chief initial orbital elements
a_0 = 6771000 ;         % semi-major axis [m]
e_0 = 0.0005;           % eccentricity [-]
i_0 = deg2rad(51.64);   % inclination [rad]
W_0 = deg2rad(257);     % RAAN [rad]
w_0 = deg2rad(45);      % argument of perigee [rad]
f_0 = deg2rad(30);      % true anomaly [rad]

% Chief initial state in OE & ECI
initial_state_0_OE = [a_0, e_0, i_0, W_0, w_0, f_0];
initial_state_0_ECI = utils.OE2ECI(initial_state_0_OE, mu);
r0_init = initial_state_0_ECI(1:3);
v0_init = initial_state_0_ECI(4:6);

% Deputy initial relative quasi-non-singular OE
% initial_state_1_rQNS_OE = [0; deg2rad(-0.05); 0; 0; 0; 0];
%initial_state_1_rQNS_OE = [0; deg2rad(-0.5); 0; 0; 0; 0];
% delta_a = 0;
% delta_lambda = deg2rad(-0.0700);
% %delta_lambda = deg2rad(-10.0);
% delta_e_x = 0.0007;
% delta_e_y = 0.0007;
% delta_i_x = deg2rad(0.4989);
% delta_i_y = deg2rad(0.7850);


delta_a = 0;
delta_lambda = deg2rad(-3.3773);
delta_e_x = 0.0007;
delta_e_y = 0.0007;
delta_i_x = deg2rad(0.4989);
delta_i_y = deg2rad(0.7850);


initial_state_1_rQNS_OE = [delta_a; delta_lambda; delta_e_x; delta_e_y; delta_i_x; delta_i_y];

% Deputy initial state in OE & ECI
initial_state_1_OE = utils.rQNSOE2OE(initial_state_0_OE,initial_state_1_rQNS_OE);
initial_state_1_ECI = utils.OE2ECI(initial_state_1_OE, mu);
r1_init = initial_state_1_ECI(1:3);
v1_init = initial_state_1_ECI(4:6);


% Compute the distance between the two spacecraft
distance = norm(r1_init - r0_init);
% Display the distance
fprintf('Distance between the two spacecraft: %.2f meters\n', distance);


% Desired rQNSOE
epsilon = 1e-7;
desired_QNSROE = [epsilon; epsilon; epsilon; epsilon; epsilon; epsilon];
tau_diff = 1e-3; % Small time step for A_c numerical differentiation (seconds)
k = 1e3; % Scaling factor for P
N_param = 14; % Thrust concentration exponent 

%% - Generate Ground Truth State

% Orbital period
T = 2 * pi * sqrt(a_0^3 / mu); % Orbital period [s]

% Simulation parameters
tstart = 0;
tint = 10; % Match UKF time step
tend = 20 * T;
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Run the perturbed FODE propagator for chief in ECI
perturbated = true;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t, state, mu, R, J2, perturbated);
[t, state_ECI_0_p] = ode113(odefun, (tstart:tint:tend)', initial_state_0_ECI, options);

% Convert chief ECI to orbital elements and precompute A_c, B_c
N = length(t);
chief_OE_traj = zeros(N, 6);
A_all = zeros(6, 6, N);
B_all = zeros(6, 3, N);
for idx = 1:N
    state_ECI = state_ECI_0_p(idx, :)';
    chief_OE_traj(idx, :) = utils.ECI2OE(state_ECI, mu);
    Phi_J2_qns = control_helpers.getQNS_J2_STM(chief_OE_traj(idx, :), tau_diff, mu, J2, R);
    A_all(:, :, idx) = (Phi_J2_qns - eye(6)) / tau_diff;
    B_all(:, :, idx) = control_helpers.getBcMatrix(chief_OE_traj(idx, :), mu, J2, R);
end

% Propagate deputy QNS ROE with Lyapunov control after nb_orbit
nb_orbit = 2;
k_control_start = find(t >= nb_orbit * T, 1, 'first'); 
delta_alpha = initial_state_1_rQNS_OE;
delta_alpha_traj = zeros(N, 6);
u_traj = zeros(3, N); % Store control inputs
V_traj = zeros(N, 1); % Store Lyapunov function
delta_alpha_traj(1, :) = delta_alpha';
for idx = 1:N-1
    A_c = A_all(:, :, idx);
    B_c = B_all(:, :, idx);
    OE_curr = chief_OE_traj(idx, :);
    err = delta_alpha - desired_QNSROE;
    V_traj(idx) = 0.5 * (err' * err);
    if idx >= k_control_start
        P = control_helpers.getPMatrix(OE_curr, delta_alpha, desired_QNSROE, k, N_param);
        u_d = -pinv(B_c) * (A_c * delta_alpha + P * err);
    else
        u_d = zeros(3,1); % No control until nb_orbit
    end
    u_traj(:, idx) = u_d;
    delta_alpha_dot = A_c * delta_alpha + B_c * u_d;
    delta_alpha = delta_alpha + delta_alpha_dot * tint;
    delta_alpha_traj(idx+1, :) = delta_alpha';
end
% Compute final V and u
V_traj(N) = 0.5 * (delta_alpha - desired_QNSROE)' * (delta_alpha - desired_QNSROE);
if N >= k_control_start
    P = control_helpers.getPMatrix(chief_OE_traj(N, :), delta_alpha, desired_QNSROE, k, N_param);
    u_traj(:, N) = -pinv(B_all(:, :, N)) * (A_all(:, :, N) * delta_alpha + P * (delta_alpha - desired_QNSROE));
else
    u_traj(:, N) = zeros(3,1);
end

% Consolidated state
x_c_gt = state_ECI_0_p;
x_d_gt = delta_alpha_traj;
x_d_gt(:,2) = unwrap(x_d_gt(:,2));
x_gt = [x_c_gt, x_d_gt];

% Verification plots
% QNS ROE Evolution
figure('Name', 'Ground Truth QNS ROE Evolution', 'Position', [100, 100, 1200, 800]);
subplot(3, 2, 1); plot(t, delta_alpha_traj(:, 1), 'b-', 'LineWidth', 2); hold on;
plot(t, desired_QNSROE(1)*ones(N,1), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)'); ylabel('\delta a'); legend('Actual', 'Desired'); grid on; set(gca, 'FontSize', 12);
subplot(3, 2, 2); plot(t, delta_alpha_traj(:, 2), 'b-', 'LineWidth', 2); hold on;
plot(t, desired_QNSROE(2)*ones(N,1), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)'); ylabel('\delta \lambda (rad)'); legend('Actual', 'Desired'); grid on; set(gca, 'FontSize', 12);
subplot(3, 2, 3); plot(t, delta_alpha_traj(:, 3), 'b-', 'LineWidth', 2); hold on;
plot(t, desired_QNSROE(3)*ones(N,1), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)'); ylabel('\delta e_x'); legend('Actual', 'Desired'); grid on; set(gca, 'FontSize', 12);
subplot(3, 2, 4); plot(t, delta_alpha_traj(:, 4), 'b-', 'LineWidth', 2); hold on;
plot(t, desired_QNSROE(4)*ones(N,1), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)'); ylabel('\delta e_y'); legend('Actual', 'Desired'); grid on; set(gca, 'FontSize', 12);
subplot(3, 2, 5); plot(t, delta_alpha_traj(:, 5), 'b-', 'LineWidth', 2); hold on;
plot(t, desired_QNSROE(5)*ones(N,1), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)'); ylabel('\delta i_x (rad)'); legend('Actual', 'Desired'); grid on; set(gca, 'FontSize', 12);
subplot(3, 2, 6); plot(t, delta_alpha_traj(:, 6), 'b-', 'LineWidth', 2); hold on;
plot(t, desired_QNSROE(6)*ones(N,1), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)'); ylabel('\delta i_y (rad)'); legend('Actual', 'Desired'); grid on; set(gca, 'FontSize', 12);

% Lyapunov Function
figure('Name', 'Ground Truth Lyapunov Function', 'Position', [100, 100, 600, 400]);
plot(t, V_traj, 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('V'); grid on; set(gca, 'FontSize', 12);

% Control Inputs
figure('Name', 'Ground Truth Control Inputs', 'Position', [100, 100, 600, 600]);
subplot(3, 1, 1); plot(t, u_traj(1, :), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta V_r (m/s)'); grid on; set(gca, 'FontSize', 12);
subplot(3, 1, 2); plot(t, u_traj(2, :), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta V_t (m/s)'); grid on; set(gca, 'FontSize', 12);
subplot(3, 1, 3); plot(t, u_traj(3, :), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta V_n (m/s)'); grid on; set(gca, 'FontSize', 12);

% Cumulative Delta-V
deltaV_mag = sqrt(sum(u_traj.^2, 1));
cumulative_deltaV = cumsum(deltaV_mag);
fprintf('Ground Truth Final cumulative deltaV: %.4f m/s\n', cumulative_deltaV(end));
figure('Name', 'Ground Truth Cumulative Delta-V', 'Position', [100, 100, 600, 400]);
plot(t, cumulative_deltaV, 'r-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Cumulative \DeltaV (m/s)'); grid on; set(gca, 'FontSize', 12);

% Relative Distance
r_rel_traj = zeros(N,1);
for idx = 1:N
    chief_ECI = state_ECI_0_p(idx, :)';
    OE_c = chief_OE_traj(idx, :);
    delta_rQNS = delta_alpha_traj(idx, :)';
    deputy_OE = utils.rQNSOE2OE(OE_c, delta_rQNS);
    deputy_ECI = utils.OE2ECI(deputy_OE, mu);
    r_rel_traj(idx) = norm(deputy_ECI(1:3) - chief_ECI(1:3));
end
figure('Name', 'Ground Truth Relative Distance', 'Position', [100, 100, 600, 400]);
plot(t, r_rel_traj, 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Separation (m)'); grid on; set(gca, 'FontSize', 12);

%% - Generate Measurements

z_k = zeros(N, 9);
sigma_pos = 0.001; % m
sigma_vel = 0.001; % m/s
sigma_range = 0.001; % m
sigma_angles = 4.848e-3; % rad
R_k = diag([sigma_pos^2 * ones(1,3), sigma_vel^2 * ones(1,3), sigma_range^2, sigma_angles^2 * ones(1,2)]);
sqrt_R_k = sqrtm(R_k);

for k = 1:N
    x_k_gt = x_gt(k,:)';
    z_k_true = h_k(x_k_gt, zeros(9,1), mu);
    v_k = sqrt_R_k * randn(9,1);
    z_k(k,:) = z_k_true + v_k;
end

%% - UKF

% UKF Parameters
n = 12; % State dimension
m = 9;  % Measurement dimension
alpha =5*1e-1;
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
x_hat = [initial_state_0_ECI + [sigma_pos * randn(3,1); sigma_vel * randn(3,1)]; ...
         initial_state_1_rQNS_OE + [sigma_rQNS_dimless * randn(1,1); sigma_rQNS_angle * randn(1,1); ...
                                    sigma_rQNS_dimless * randn(2,1); sigma_rQNS_angle * randn(2,1)]];
P_c = 1e0 * eye(6); % Increased for robustness
P_d = 1e-1 * eye(6);
P = [P_c, zeros(6); zeros(6), P_d];

% Storage
x_hat_store = zeros(N, n);
P_store = zeros(n, n, N);
z_hat_store = zeros(N, m);
z_post_store = zeros(N, m);
u_ukf_traj = zeros(3, N);
V_ukf_traj = zeros(N, 1);
cumulative_deltaV_ukf = zeros(N, 1);

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
        if k >= k_control_start
            u_k1 = get_control_Lyapunov(chi(:,i), desired_QNSROE, k, N_param, tau_diff, mu, J2, R);
        else
            u_k1 = zeros(6,1);
        end
        chi_pred(:,i) = f_k(chi(:,i), u_k1, zeros(n,1), delta_t, mu, J2, R);
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

    % Compute control and Lyapunov function
    V_k = 0.5 * (x_hat_pred(7:12) - desired_QNSROE)' * (x_hat_pred(7:12) - desired_QNSROE);
    if k >= k_control_start
        u_k = get_control_Lyapunov(x_hat_pred, desired_QNSROE, k, N_param, tau_diff, mu, J2, R);
    else
        u_k = zeros(6,1);
    end
    u_ukf_traj(:, k) = u_k(4:6);
    V_ukf_traj(k) = V_k;
    deltaV_k = norm(u_k(4:6));
    if k > 1
        cumulative_deltaV_ukf(k) = cumulative_deltaV_ukf(k-1) + deltaV_k;
    else
        cumulative_deltaV_ukf(k) = deltaV_k;
    end

    % Update step
    sqrt_P_pred = chol((n + lambda) * P_pred, 'lower');
    chi_pred_new = zeros(n, 2*n+1);
    chi_pred_new(:,1) = x_hat_pred;
    for i = 1:n
        chi_pred_new(:,i+1) = x_hat_pred + sqrt_P_pred(:,i);
        chi_pred_new(:,i+n+1) = x_hat_pred - sqrt_P_pred(:,i);
    end
    Z_pred = zeros(m, 2*n+1);
    for i = 1:(2*n+1)
        Z_pred(:,i) = h_k(chi_pred_new(:,i), zeros(m,1), mu);
    end
    z_hat = zeros(m,1);
    for i = 1:(2*n+1)
        z_hat = z_hat + Wm(i) * Z_pred(:,i);
    end
    z_hat_store(k,:) = z_hat';

    S_k = R_k;
    P_xz = zeros(n,m);
    for i = 1:(2*n+1)
        z_diff = Z_pred(:,i) - z_hat;
        S_k = S_k + Wc(i) * (z_diff * z_diff');
        x_diff = chi_pred_new(:,i) - x_hat_pred;
        P_xz = P_xz + Wc(i) * (x_diff * z_diff');
    end
    K_k = P_xz / S_k;
    z_k_current = z_k(k,:)';
    x_hat = x_hat_pred + K_k * (z_k_current - z_hat);
    P = P_pred - K_k * S_k * K_k';

    % Store results
    x_hat_store(k,:) = x_hat';
    P_store(:,:,k) = P;
    z_post_store(k,:) = h_k(x_hat, zeros(m,1), mu)';
end


%% - Plots V2

% Chief ECI Position
figure('Name', 'Chief ECI Position', 'Position', [100, 100, 600, 600]);
subplot(3,1,1);
plot(t, x_c_gt(:,1), 'b-', 'DisplayName', 'True X', 'LineWidth', 2);
hold on;
plot(t, x_hat_store(:,1), 'r--', 'DisplayName', 'Predicted X', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('X (m)'); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);
subplot(3,1,2);
plot(t, x_c_gt(:,2), 'b-', 'DisplayName', 'True Y', 'LineWidth', 2);
hold on;
plot(t, x_hat_store(:,2), 'r--', 'DisplayName', 'Predicted Y', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Y (m)'); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);
subplot(3,1,3);
plot(t, x_c_gt(:,3), 'b-', 'DisplayName', 'True Z', 'LineWidth', 2);
hold on;
plot(t, x_hat_store(:,3), 'r--', 'DisplayName', 'Predicted Z', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Z (m)'); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);

% Chief ECI Position Error with 3-Sigma Bounds
figure('Name', 'Chief ECI Position Error', 'Position', [700, 100, 600, 600]);
subplot(3,1,1);
plot(t, x_c_gt(:,1) - x_hat_store(:,1), 'k-', 'DisplayName', 'Error X', 'LineWidth', 2);
hold on;
plot(t, 3*sqrt(squeeze(P_store(1,1,:))), 'g--', 'DisplayName', '+3\sigma', 'LineWidth', 1.5);
plot(t, -3*sqrt(squeeze(P_store(1,1,:))), 'g--', 'DisplayName', '-3\sigma', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error X (m)'); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);
subplot(3,1,2);
plot(t, x_c_gt(:,2) - x_hat_store(:,2), 'k-', 'DisplayName', 'Error Y', 'LineWidth', 2);
hold on;
plot(t, 3*sqrt(squeeze(P_store(2,2,:))), 'g--', 'DisplayName', '+3\sigma', 'LineWidth', 1.5);
plot(t, -3*sqrt(squeeze(P_store(2,2,:))), 'g--', 'DisplayName', '-3\sigma', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error Y (m)'); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);
subplot(3,1,3);
plot(t, x_c_gt(:,3) - x_hat_store(:,3), 'k-', 'DisplayName', 'Error Z', 'LineWidth', 2);
hold on;
plot(t, 3*sqrt(squeeze(P_store(3,3,:))), 'g--', 'DisplayName', '+3\sigma', 'LineWidth', 1.5);
plot(t, -3*sqrt(squeeze(P_store(3,3,:))), 'g--', 'DisplayName', '-3\sigma', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error Z (m)'); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);

% Chief ECI Velocity
figure('Name', 'Chief ECI Velocity', 'Position', [100, 700, 600, 600]);
subplot(3,1,1);
plot(t, x_c_gt(:,4), 'b-', 'DisplayName', 'True V_x', 'LineWidth', 2);
hold on;
plot(t, x_hat_store(:,4), 'r--', 'DisplayName', 'Predicted V_x', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('V_x (m/s)'); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);
subplot(3,1,2);
plot(t, x_c_gt(:,5), 'b-', 'DisplayName', 'True V_y', 'LineWidth', 2);
hold on;
plot(t, x_hat_store(:,5), 'r--', 'DisplayName', 'Predicted V_y', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('V_y (m/s)'); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);
subplot(3,1,3);
plot(t, x_c_gt(:,6), 'b-', 'DisplayName', 'True V_z', 'LineWidth', 2);
hold on;
plot(t, x_hat_store(:,6), 'r--', 'DisplayName', 'Predicted V_z', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('V_z (m/s)'); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);

% Chief ECI Velocity Error with 3-Sigma Bounds
figure('Name', 'Chief ECI Velocity Error', 'Position', [700, 700, 600, 600]);
subplot(3,1,1);
plot(t, x_c_gt(:,4) - x_hat_store(:,4), 'k-', 'DisplayName', 'Error V_x', 'LineWidth', 2);
hold on;
plot(t, 3*sqrt(squeeze(P_store(4,4,:))), 'g--', 'DisplayName', '+3\sigma', 'LineWidth', 1.5);
plot(t, -3*sqrt(squeeze(P_store(4,4,:))), 'g--', 'DisplayName', '-3\sigma', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error V_x (m/s)'); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);
subplot(3,1,2);
plot(t, x_c_gt(:,5) - x_hat_store(:,5), 'k-', 'DisplayName', 'Error V_y', 'LineWidth', 2);
hold on;
plot(t, 3*sqrt(squeeze(P_store(5,5,:))), 'g--', 'DisplayName', '+3\sigma', 'LineWidth', 1.5);
plot(t, -3*sqrt(squeeze(P_store(5,5,:))), 'g--', 'DisplayName', '-3\sigma', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error V_y (m/s)'); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);
subplot(3,1,3);
plot(t, x_c_gt(:,6) - x_hat_store(:,6), 'k-', 'DisplayName', 'Error V_z', 'LineWidth', 2);
hold on;
plot(t, 3*sqrt(squeeze(P_store(6,6,:))), 'g--', 'DisplayName', '+3\sigma', 'LineWidth', 1.5);
plot(t, -3*sqrt(squeeze(P_store(6,6,:))), 'g--', 'DisplayName', '-3\sigma', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error V_z (m/s)'); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);

% Deputy QNS ROE with Desired Reference
figure('Name', 'Deputy QNS ROE', 'Position', [100, 1300, 1200, 800]);
labels = {'\delta a', '\delta \lambda (rad)', '\delta e_x', '\delta e_y', '\delta i_x (rad)', '\delta i_y (rad)'};
for i = 1:6
    subplot(3,2,i);
    plot(t, x_d_gt(:,i), 'b-', 'DisplayName', 'True', 'LineWidth', 2); hold on;
    plot(t, x_hat_store(:,6+i), 'r--', 'DisplayName', 'Predicted', 'LineWidth', 1.5);
    plot(t, desired_QNSROE(i)*ones(N,1), 'k:', 'DisplayName', 'Desired', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel(labels{i}); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);
end

% Deputy QNS ROE Error with 3-Sigma Bounds
figure('Name', 'Deputy QNS ROE Error', 'Position', [700, 1300, 1200, 800]);
for i = 1:6
    subplot(3,2,i);
    plot(t, x_d_gt(:,i) - x_hat_store(:,6+i), 'k-', 'DisplayName', sprintf('Error %s', labels{i}), 'LineWidth', 2);
    hold on;
    plot(t, 3*sqrt(squeeze(P_store(6+i,6+i,:))), 'g--', 'DisplayName', '+3\sigma', 'LineWidth', 1.5);
    plot(t, -3*sqrt(squeeze(P_store(6+i,6+i,:))), 'g--', 'DisplayName', '-3\sigma', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel(sprintf('Error %s', labels{i})); legend('Location', 'best'); grid on; set(gca, 'FontSize', 12);
end

% Lyapunov Function for Estimated QNS ROE
figure('Name', 'UKF Lyapunov Function', 'Position', [100, 1900, 600, 400]);
plot(t, V_ukf_traj, 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('V'); grid on; set(gca, 'FontSize', 12);

% UKF Control Inputs
figure('Name', 'UKF Control Inputs', 'Position', [700, 1900, 600, 600]);
subplot(3,1,1); plot(t, u_ukf_traj(1,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta V_r (m/s)'); grid on; set(gca, 'FontSize', 12);
subplot(3,1,2); plot(t, u_ukf_traj(2,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta V_t (m/s)'); grid on; set(gca, 'FontSize', 12);
subplot(3,1,3); plot(t, u_ukf_traj(3,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta V_n (m/s)'); grid on; set(gca, 'FontSize', 12);

% UKF Cumulative Delta-V
deltaV_ukf_mag = sqrt(sum(u_ukf_traj.^2, 1));
cumulative_deltaV_ukf = cumsum(deltaV_ukf_mag);
figure('Name', 'UKF Cumulative Delta-V', 'Position', [100, 2500, 600, 400]);
plot(t, cumulative_deltaV_ukf, 'r-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Cumulative \DeltaV (m/s)'); grid on; set(gca, 'FontSize', 12);
fprintf('UKF Final cumulative deltaV: %.4f m/s\n', cumulative_deltaV_ukf(end));

%% - Additional Plots

% True statistics (mean and std of state error in last orbit)
last_orbit_idx = t >= (19*T); % Last orbit (14T to 15T)
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


%% - Control Statistics

% Calculate differences between ground truth and UKF control inputs
control_error = u_traj - u_ukf_traj; % Ground truth - UKF (RTN components)

% Calculate cumulative Delta V for ground truth and UKF (already computed)
% deltaV_mag = sqrt(sum(u_traj.^2, 1));
% cumulative_deltaV = cumsum(deltaV_mag);
% deltaV_ukf_mag = sqrt(sum(u_ukf_traj.^2, 1));
% cumulative_deltaV_ukf = cumsum(deltaV_ukf_mag);

% Calculate cumulative Delta V difference
cumulative_deltaV_diff = cumulative_deltaV - cumulative_deltaV_ukf;

% Calculate Lyapunov function difference
lyapunov_diff = V_traj - V_ukf_traj; % Ground truth - UKF

% Plot Delta V difference in RTN over time
figure('Name', 'Control Delta V Difference (RTN)', 'Position', [100, 3100, 600, 600]);
subplot(3,1,1);
plot(t, control_error(1,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta V_r Error (m/s)'); title('\Delta V_r (Ground Truth - UKF)'); grid on; set(gca, 'FontSize', 12);
subplot(3,1,2);
plot(t, control_error(2,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta V_t Error (m/s)'); title('\Delta V_t (Ground Truth - UKF)'); grid on; set(gca, 'FontSize', 12);
subplot(3,1,3);
plot(t, control_error(3,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta V_n Error (m/s)'); title('\Delta V_n (Ground Truth - UKF)'); grid on; set(gca, 'FontSize', 12);

% Plot cumulative Delta V difference over time
figure('Name', 'Cumulative Delta V Difference', 'Position', [700, 3100, 600, 400]);
plot(t, cumulative_deltaV_diff, 'r-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Cumulative \DeltaV Difference (m/s)'); title('Cumulative \DeltaV (Ground Truth - UKF)'); grid on; set(gca, 'FontSize', 12);
fprintf('Final cumulative Delta V difference (Ground Truth - UKF): %.4f m/s\n', cumulative_deltaV_diff(end));

% Plot Lyapunov function difference over time
figure('Name', 'Lyapunov Function Difference', 'Position', [100, 3600, 600, 400]);
plot(t, lyapunov_diff, 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Lyapunov Function Difference'); title('Lyapunov Function (Ground Truth - UKF)'); grid on; set(gca, 'FontSize', 12);

% Diagnostic: Check maximum control input differences
max_control_error = max(abs(control_error), [], 2);
fprintf('Maximum control input errors (RTN):\n');
fprintf('  \Delta V_r: %.4e m/s\n', max_control_error(1));
fprintf('  \Delta V_t: %.4e m/s\n', max_control_error(2));
fprintf('  \Delta V_n: %.4e m/s\n', max_control_error(3));


%% - Functions
function x_k = f_k(x_k1, u_k, w_k, delta_t, mu, J2, R)
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
    chief_OE_k1 = utils.ECI2OE(x_k1(1:6), mu);
    Phi_J2_qns_k1 = control_helpers.getQNS_J2_STM(chief_OE_k1, delta_t, mu, J2, R);
    B_k1 = control_helpers.getBcMatrix(chief_OE_k1, mu, J2, R);    
    f_k_d = Phi_J2_qns_k1*x_d_k1 + delta_t*B_k1*u_k_d + w_k_d;
    x_k = [f_k_c; f_k_d];
end


function z_k = h_k(x_k, v_k, mu)
    x_c_k = x_k(1:6);
    r_c_k = x_k(1:3);
    v_c_k = x_k(4:6);
    rQNSOE_d_k = x_k(7:12);
    OE_c_k = utils.ECI2OE(x_c_k,mu);
    OE_d_k = utils.rQNSOE2OE(OE_c_k,rQNSOE_d_k);
    x_d_k = utils.OE2ECI(OE_d_k,mu);
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


function u_k = get_control_Lyapunov(x_k, desired_QNSROE, k, N_param, tau_diff, mu, J2, R)
    % Compute Lyapunov-based control input for deputy satellite
    % Inputs:
    %   x_k: State vector [r_c; v_c; delta_a; delta_lambda; delta_e_x; delta_e_y; delta_i_x; delta_i_y]
    %   desired_QNSROE: Desired QNS ROE [6x1]
    %   const: Constants structure
    %   body: Celestial body ('earth')
    %   k: Lyapunov gain scaling factor
    %   N_param: Thrust concentration exponent
    %   tau_diff: Small time step for A_c approximation (s)
    % Outputs:
    %   u_k: Control input [u_c; u_d] = [0; 0; 0; DeltaV_r; DeltaV_t; DeltaV_n]

    % Extract deputy QNS ROE
    delta_alpha = x_k(7:12);

    % Compute chief orbital elements
    OE_curr = utils.ECI2OE(x_k(1:6), mu);

    % Compute A_c and B_c
    Phi_J2_qns = control_helpers.getQNS_J2_STM(OE_curr, tau_diff, mu, J2, R);
    A_c = (Phi_J2_qns - eye(6)) / tau_diff;
    B_c = control_helpers.getBcMatrix(OE_curr, mu, J2, R);

    % Compute P matrix
    P = control_helpers.getPMatrix(OE_curr, delta_alpha, desired_QNSROE, k, N_param);

    % Compute tracking error
    err = delta_alpha - desired_QNSROE;

    % Compute deputy control input
    u_d = -pinv(B_c) * (A_c * delta_alpha + P * err);

    % Construct full control input (no control for chief)
    u_k = [zeros(3,1); u_d];
end