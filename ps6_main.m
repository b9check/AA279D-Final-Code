%% - Reset

clc; clear; close all;

%% - Initialization

% Orbited body
body = 'earth';
const = utils.getConstants({body});

% Chief initial orbital elements
a_0 = 6771000 ;         % semi-major axis [m]
e_0 = 0.0005;           % eccentricity [-]
i_0 = deg2rad(51.64);   % inclination [rad]
W_0 = deg2rad(257);     % RAAN [rad]
w_0 = deg2rad(45);       % argument of perigee [rad]
f_0 = deg2rad(30);      % true anomaly [rad]

% Chief initial state in OE & ECI
initial_state_0_OE = [a_0, e_0, i_0, W_0, w_0, f_0];
initial_state_0_ECI = utils.OE2ECI(initial_state_0_OE, const, body);
r0_init = initial_state_0_ECI(1:3);
v0_init = initial_state_0_ECI(4:6);

% Deputy initial relative quasi-non-singular OE
% delta_a = 0/a_0;
% delta_lambda = 100/a_0;
% delta_e_x = 50/a_0;
% delta_e_y = 100/a_0;
% delta_i_x = 30/a_0;
% delta_i_y = 200/a_0;
% 
% delta_a = 5000/a_0;
% delta_lambda = 150/a_0;
% delta_e_x = 300/a_0;
% delta_e_y = 200/a_0;
% delta_i_x = 70/a_0;
% delta_i_y = 150/a_0;
% 
delta_a = 0;
delta_lambda = deg2rad(-3.3773);
delta_e_x = 0.0007;
delta_e_y = 0.0007;
delta_i_x = deg2rad(0.4989);
delta_i_y = deg2rad(0.7850);

% delta_a = 0;
% delta_lambda = deg2rad(-0.0700);
% delta_e_x = 0.0007;
% delta_e_y = 0.0007;
% delta_i_x = deg2rad(0.4989);
% delta_i_y = deg2rad(0.7850);


initial_state_1_rQNS_OE = [delta_a, delta_lambda, delta_e_x, delta_e_y, delta_i_x, delta_i_y];

% Deputy initial state in OE & ECI
initial_state_1_OE = rQNSOE2OE(initial_state_0_OE,initial_state_1_rQNS_OE);
initial_state_1_ECI = utils.OE2ECI(initial_state_1_OE, const, body);
r1_init = initial_state_1_ECI(1:3);
v1_init = initial_state_1_ECI(4:6);

% Compute the distance between the two spacecraft
distance = norm(r1_init - r0_init);

% Display the distance
fprintf('Distance between the two spacecraft: %.2f meters\n', distance);

% Desired rQNSOE
%desired_ROE = [0; deg2rad(-0.05); 0; 0; 0; 0];
%desired_ROE = [0; 0; 0; 0; 0; 0];
%desired_QNSROE = OE2QNSOE(desired_ROE');
%desired_QNSROE   = OE2QNSOE(desired_ROE')';
%desired_QNSROE = [0; 0; 0; 0; 0; 0];
epsilon = 1e-7;
desired_QNSROE = [epsilon; epsilon; epsilon; epsilon; epsilon; epsilon];


%% - Groundtruth Chief & Control Matrices

% Orbital period
mu = const.(body).mu;
T = 2 * pi * sqrt(a_0^3 / mu);

% Simulation parameters
tstart = 0;
tint = 10;
tend = 20*T;
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Run the perturbed FODE propagator for chief in ECI
perturbated = true;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t,state,const,body,perturbated);                               
[t, state_ECI_0_p] = ode113(odefun, (tstart:tint:tend)', initial_state_0_ECI, options);


% Preallocate storage for A and B matrices at each time step
N = length(t);
A_all = zeros(6, 6, N);
B_all = zeros(6, 3, N);

% Small time step for A_c numerical differentiation
tau_diff = 1e-3; % seconds

for idx = 1:N
    % Get current chief state in ECI and convert to OE
    state_ECI = state_ECI_0_p(idx, :)';
    OE_curr = utils.ECI2OE(state_ECI, const, body);

    % Compute STM over small tau_diff and approximate A_c
    Phi_J2_qns = getQNS_J2_STM(OE_curr, tau_diff, const, body);
    A_all(:, :, idx) = (Phi_J2_qns - eye(6)) / tau_diff;

    % Compute B_c
    B_all(:, :, idx) = getBcMatrix(OE_curr, const, body);
end


%% -  Noiseless Lyapunov

% Convert ECI trajectory to orbital elements
N = length(t);
chief_OE_traj = zeros(N, 6);
for idx = 1:N
    state_ECI = state_ECI_0_p(idx, :)';
    chief_OE_traj(idx, :) = utils.ECI2OE(state_ECI, const, body);
end

% Verify A_c(2,1) against expected delta_lambda drift
% n_c = sqrt(mu / a_0^3);
% fprintf('Expected A_c(2,1) ≈ %e, Actual A_c(2,1) = %e\n', -1.5*n_c, A_all(2,1,1));

% Control parameters
k = 1e3; % Scaling factor for P
N_param = 14; % Thrust concentration exponent (avoid conflict with N)
dt = tint; % Time step for Euler integration (seconds)

% Initialize deputy ROE and storage
%delta_alpha = initial_state_1_rQNS_OE; % Current ROE
delta_alpha     = initial_state_1_rQNS_OE';
delta_alpha_traj = zeros(N, 6); % Store ROE trajectory
u_traj = zeros(3, N); % Store control inputs
V_traj = zeros(N, 1); % Store Lyapunov function
delta_alpha_traj(1, :) = delta_alpha'; % Initial state


% Step-by-step control loop
for idx = 1:N-1
    % Current time and matrices
    t_curr = t(idx);
    A_c = A_all(:, :, idx);
    B_c = B_all(:, :, idx);
    OE_curr = chief_OE_traj(idx, :);

    % Compute P
    P = getPMatrix(OE_curr, delta_alpha, desired_QNSROE, k, N_param);

    % Compute Lyapunov function
    err = delta_alpha - desired_QNSROE;
    V_traj(idx)= 0.5 * (err' * err);
    % V_traj(idx) = 0.5 * (delta_alpha - desired_QNSROE)' * (delta_alpha - desired_QNSROE);

    % Compute control input (fixed desired ROE, no RG)
    %u = -pinv(B_c) * (A_c * delta_alpha + P * (delta_alpha - desired_QNSROE));
    u = -pinv(B_c) * ( A_c*delta_alpha + P*err );
    u_traj(:, idx) = u;

    % Update ROE using Euler integration
    %delta_alpha_dot = A_c * delta_alpha + B_c * u;
    %delta_alpha = delta_alpha + delta_alpha_dot * dt;
    %delta_alpha_traj(idx+1, :) = delta_alpha';
    delta_alpha_dot = A_c*delta_alpha + B_c*u; % 6×1
    delta_alpha = delta_alpha + delta_alpha_dot*dt;
    delta_alpha_traj(idx+1,:) = delta_alpha';         % row storage
end

% Compute final V and u
V_traj(N) = 0.5 * (delta_alpha - desired_QNSROE)' * (delta_alpha - desired_QNSROE);
u_traj(:, N) = -pinv(B_all(:, :, N)) * (A_all(:, :, N) * delta_alpha + ...
    getPMatrix(chief_OE_traj(N, :), delta_alpha, desired_QNSROE, k, N_param) * (delta_alpha - desired_QNSROE));

% Plot results
% ROE trajectories with desired ROE
figure('Name', 'QNS ROE Evolution', 'Position', [100, 100, 1200, 800]);
subplot(3, 2, 1); plot(t, delta_alpha_traj(:, 1), 'b-', 'LineWidth', 2); hold on;
plot(t, desired_QNSROE(1)*ones(N,1), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)'); ylabel('\delta a'); %title('Relative Semi-Major Axis');
legend('Actual', 'Desired'); grid on; set(gca, 'FontSize', 12);
subplot(3, 2, 2); plot(t, delta_alpha_traj(:, 2), 'b-', 'LineWidth', 2); hold on;
plot(t, desired_QNSROE(2)*ones(N,1), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)'); ylabel('\delta \lambda (rad)'); %title('Relative Mean Longitude');
legend('Actual', 'Desired'); grid on; set(gca, 'FontSize', 12);
subplot(3, 2, 3); plot(t, delta_alpha_traj(:, 3), 'b-', 'LineWidth', 2); hold on;
plot(t, desired_QNSROE(3)*ones(N,1), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)'); ylabel('\delta e_x'); %title('Eccentricity Vector X');
legend('Actual', 'Desired'); grid on; set(gca, 'FontSize', 12);
subplot(3, 2, 4); plot(t, delta_alpha_traj(:, 4), 'b-', 'LineWidth', 2); hold on;
plot(t, desired_QNSROE(4)*ones(N,1), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)'); ylabel('\delta e_y'); %title('Eccentricity Vector Y');
legend('Actual', 'Desired'); grid on; set(gca, 'FontSize', 12);
subplot(3, 2, 5); plot(t, delta_alpha_traj(:, 5), 'b-', 'LineWidth', 2); hold on;
plot(t, desired_QNSROE(5)*ones(N,1), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)'); ylabel('\delta i_x'); %title('Inclination Vector X');
legend('Actual', 'Desired'); grid on; set(gca, 'FontSize', 12);
subplot(3, 2, 6); plot(t, delta_alpha_traj(:, 6), 'b-', 'LineWidth', 2); hold on;
plot(t, desired_QNSROE(6)*ones(N,1), 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)'); ylabel('\delta i_y'); %title('Inclination Vector Y');
legend('Actual', 'Desired'); grid on; set(gca, 'FontSize', 12);

% Lyapunov function
figure('Name', 'Lyapunov Function', 'Position', [100, 100, 600, 400]);
plot(t, V_traj, 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('V'); %title('Lyapunov Function Over Time');
grid on; set(gca, 'FontSize', 12);

% Control inputs
figure('Name', 'Control Inputs', 'Position', [100, 100, 600, 600]);
subplot(3, 1, 1); plot(t, u_traj(1, :), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta V_r (m/s)'); %title('Radial Control Input');
grid on; set(gca, 'FontSize', 12);
subplot(3, 1, 2); plot(t, u_traj(2, :), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta V_t (m/s)'); %title('Tangential Control Input');
grid on; set(gca, 'FontSize', 12);
subplot(3, 1, 3); plot(t, u_traj(3, :), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\Delta V_n (m/s)'); %title('Normal Control Input');
grid on; set(gca, 'FontSize', 12);

deltaV_mag = sqrt(sum(u_traj.^2, 1));
cumulative_deltaV = cumsum(deltaV_mag);
fprintf('Final cumulative deltaV: %.4f m/s\n', cumulative_deltaV(end));
figure('Name', 'Cumulative Delta-V', 'Position', [100, 100, 600, 400]);
plot(t, cumulative_deltaV, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Cumulative \DeltaV (m/s)');
%title('Total Cumulative \DeltaV over Time');
grid on; set(gca, 'FontSize', 12);

% Relative distance 
r_rel_traj = zeros(N,1);
for idx = 1:N
    chief_ECI   = state_ECI_0_p(idx, :)';
    OE_c        = chief_OE_traj(idx, :);
    delta_rQNS  = delta_alpha_traj(idx, :);
    deputy_OE   = rQNSOE2OE(OE_c, delta_rQNS);
    deputy_ECI  = utils.OE2ECI(deputy_OE, const, body);
    r_rel_traj(idx) = norm(deputy_ECI(1:3) - chief_ECI(1:3));
end

figure('Name','Relative Distance','Position',[100,100,600,400]);
plot(t, r_rel_traj, 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Separation (m)');
%title('Chief–Deputy Relative Distance');
grid on; set(gca,'FontSize',12);

% Save results
save('control_results.mat', 't', 'delta_alpha_traj', 'V_traj', ...
     'u_traj', 'A_all', 'B_all', 'chief_OE_traj', 'r_rel_traj');


%% - Noisy Sensors and Actuators

% Additive white Gaussian noise on measurement and actuation

% Noise parameters
sigma_sensor = [1e-5; 1e-5; 1e-6; 1e-6; 1e-7; 1e-7]; % std dev for each ROE component
sigma_act     = [1e-5; 1e-5; 1e-5];                 % std dev for each ΔV channel (m/s)

% Preallocate storage
delta_alpha_noisy_traj = zeros(N,6);
u_noisy_traj           = zeros(3,N);
V_noisy_traj           = zeros(N,1);

% Initialize
delta_alpha_true       = initial_state_1_rQNS_OE';  % true ROE state
delta_alpha_noisy_traj(1,:) = delta_alpha_true';
V_noisy_traj(1)        = 0.5 * norm(delta_alpha_true - desired_QNSROE)^2;

% Reseed for reproducibility
rng(0);

for idx = 1:N-1
    % 1) True dynamics matrices at this step
    A_c = A_all(:,:,idx);
    B_c = B_all(:,:,idx);
    
    % 2) Sensor measurement with noise
    meas_delta = delta_alpha_true + sigma_sensor .* randn(6,1);

    % 3) Compute control based on noisy measurement
    err_meas = meas_delta - desired_QNSROE;
    P       = getPMatrix(chief_OE_traj(idx,:), meas_delta, desired_QNSROE, k, N_param);
    u_nom   = -pinv(B_c)*(A_c*meas_delta + P*err_meas);

    % 4) Apply actuator noise
    u_applied = u_nom + sigma_act .* randn(3,1);
    u_noisy_traj(:,idx) = u_applied;

    % 5) Propagate true ROE with noisy actuation
    delta_dot      = A_c*delta_alpha_true + B_c*u_applied;
    delta_alpha_true = delta_alpha_true + delta_dot*dt;
    delta_alpha_noisy_traj(idx+1,:) = delta_alpha_true';

    % 6) Noisy‐based Lyapunov evaluation (true state vs. desired)
    V_noisy_traj(idx+1) = 0.5 * norm(delta_alpha_true - desired_QNSROE)^2;
end

% Final actuator at N
u_noisy_traj(:,N) = -pinv(B_all(:,:,N))*(A_all(:,:,N)*delta_alpha_true + ...
    getPMatrix(chief_OE_traj(N,:), delta_alpha_true, desired_QNSROE, k, N_param)*(delta_alpha_true - desired_QNSROE));

% ROE Evolution: Noise‐Free vs. Noisy
figure('Name','ROE Evolution: Noise‐Free vs. Noisy','Position',[100,100,1200,800]);
for comp = 1:6
    subplot(3,2,comp);
    % noise-free
    plot(t, delta_alpha_traj(:,comp), 'b-', 'LineWidth', 2); hold on;
    % noisy
    plot(t, delta_alpha_noisy_traj(:,comp), '--', 'LineWidth', 1.5);
    % desired
    plot(t, desired_QNSROE(comp)*ones(N,1), 'r--', 'LineWidth', 1.5);
    hold off;

    % switch comp
    %     case 1, title('Relative Semi‐Major Axis (\delta a)'),            ylabel('\delta a');
    %     case 2, title('Relative Mean Longitude (\delta \lambda)'),     ylabel('\delta \lambda (rad)');
    %     case 3, title('Eccentricity Vector X (\delta e_x)'),           ylabel('\delta e_x');
    %     case 4, title('Eccentricity Vector Y (\delta e_y)'),           ylabel('\delta e_y');
    %     case 5, title('Inclination Vector X (\delta i_x)'),            ylabel('\delta i_x');
    %     case 6, title('Inclination Vector Y (\delta i_y)'),            ylabel('\delta i_y');
    % end
    switch comp
        case 1,            ylabel('\delta a');
        case 2,   ylabel('\delta \lambda (rad)');
        case 3,  ylabel('\delta e_x');
        case 4,   ylabel('\delta e_y');
        case 5,   ylabel('\delta i_x');
        case 6,  ylabel('\delta i_y');
    end

    xlabel('Time (s)');
    legend('Noise‐Free','Noisy','Desired','Location','Best');
    grid on; set(gca,'FontSize',12);
end

% Control Inputs: Noise‐Free vs. Noisy
figure('Name','Control Inputs: Noise‐Free vs. Noisy','Position',[100,100,600,600]);
channels = {'Radial','Tangential','Normal'};
labels   = {'\Delta V_r (m/s)','\Delta V_t (m/s)','\Delta V_n (m/s)'};
for kidx = 1:3
    subplot(3,1,kidx);
    plot(t, u_traj(kidx,:),    'b-',  'LineWidth', 2); hold on;
    plot(t, u_noisy_traj(kidx,:),'--',  'LineWidth', 1.5); hold off;
    xlabel('Time (s)'); ylabel(labels{kidx});
    %title([channels{kidx} ' Control Input']);
    legend('Noise‐Free','Noisy','Location','Best');
    grid on; set(gca,'FontSize',12);
end


%Plot comparison of Lyapunov functions
figure('Name','Lyapunov: Noise-Free vs. Noisy','Position',[100,100,800,400]);
plot(t, V_traj, 'LineWidth',2); hold on;
plot(t, V_noisy_traj, '--', 'LineWidth',1.5); hold off;
xlabel('Time (s)'); ylabel('V'); % title('Lyapunov Function Comparison');
legend('Noise-Free','With Sensor \& Actuator Noise','Location','Best');
grid on; set(gca,'FontSize',12);

% Plot separation error growth
figure('Name','Separation Error','Position',[100,100,800,400]);
% Noise-free separation
plot(t, r_rel_traj, 'LineWidth',2); hold on;
% Noisy separation: recompute deputy ECI and separation
r_rel_noisy = zeros(N,1);
for idx=1:N
    OE_c     = chief_OE_traj(idx,:);
    delta_rQ = delta_alpha_noisy_traj(idx,:);
    deputy_OE = rQNSOE2OE(OE_c, delta_rQ);
    dep_ECI   = utils.OE2ECI(deputy_OE, const, body);
    chief_ECI = state_ECI_0_p(idx,:)';
    r_rel_noisy(idx) = norm(dep_ECI(1:3)-chief_ECI(1:3));
end
plot(t, r_rel_noisy, '--','LineWidth',1.5); hold off;
xlabel('Time (s)'); ylabel('Separation (m)');
%title('Chief–Deputy Separation: Noise-Free vs. Noisy');
legend('Noise-Free','With Noise','Location','Best');
grid on; set(gca,'FontSize',12);

% Cumulative delta-V: Noise-Free vs. Noisy
deltaV_mag_free = sqrt(sum(u_traj.^2, 1));       % norm of control input at each step
cumulative_deltaV_free = cumsum(deltaV_mag_free);

deltaV_mag_noisy = sqrt(sum(u_noisy_traj.^2, 1)); % noisy case
cumulative_deltaV_noisy = cumsum(deltaV_mag_noisy);

fprintf('Final cumulative deltaV (noise-free): %.4f m/s\n', cumulative_deltaV_free(end));
fprintf('Final cumulative deltaV (noisy):      %.4f m/s\n', cumulative_deltaV_noisy(end));

figure('Name', 'Cumulative Delta-V Comparison', 'Position', [100, 100, 600, 400]);
plot(t, cumulative_deltaV_free, 'b-', 'LineWidth', 2); hold on;
plot(t, cumulative_deltaV_noisy, 'r--', 'LineWidth', 1.5); hold off;
xlabel('Time (s)');
ylabel('Cumulative \DeltaV (m/s)');
legend('Noise-Free','Noisy','Location','Best');
grid on; set(gca, 'FontSize', 12);



%% - Functions

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


function P_bar = getPMatrix(chief_OE, current_QNSROE, desired_QNSROE, k, N)
    % Extract chief orbital elements
    a = chief_OE(1); % Semi-major axis (m)
    e = chief_OE(2); % Eccentricity (-)
    i = chief_OE(3); % Inclination (rad)
    W = chief_OE(4); % RAAN (rad)
    w = chief_OE(5); % Argument of perigee (rad)
    f = chief_OE(6); % True anomaly (rad)

    % Compute mean anomaly (M) from true anomaly (f)
    E = 2 * atan2(sqrt(1 - e) * sin(f/2), sqrt(1 + e) * cos(f/2)); % Eccentric anomaly
    M = E - e * sin(E); % Mean anomaly

    % Mean argument of latitude
    phi = w + M;

    % Compute tracking errors
    delta_rQNSOE = current_QNSROE - desired_QNSROE;
    delta_rQNSOE_ex = delta_rQNSOE(3);
    delta_rQNSOE_ey = delta_rQNSOE(4);
    delta_rQNSOE_ix = delta_rQNSOE(5);
    delta_rQNSOE_iy = delta_rQNSOE(6);

    % Compute fuel-optimal angles
    phi_ip = atan2(delta_rQNSOE_ey, delta_rQNSOE_ex); % In-plane optimal angle
    phi_oop = atan2(delta_rQNSOE_iy, delta_rQNSOE_ix); % Out-of-plane optimal angle
    phi_lambda = w; % Radial thrust at perigee (approximate)

    % Angular offsets
    J = phi - phi_ip;
    H = phi - phi_oop;
    K = phi - phi_lambda;

    % Ensure N is even and positive
    if mod(N, 2) ~= 0 || N < 2
        error('N must be an even positive integer >= 2');
    end

    % Construct P_bar
    P_bar = (1/k) * diag([cos(J)^N, cos(K)^N, cos(J)^N, cos(J)^N, cos(H)^N, cos(H)^N]);
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


function QNSOE = OE2QNSOE(OE_array)
    N = size(OE_array, 1);
    QNSOE = zeros(N, 6);

    for k = 1:N
        a = OE_array(k,1);
        e = OE_array(k,2);
        i = OE_array(k,3);
        W = OE_array(k,4);
        w = OE_array(k,5);
        M = OE_array(k,6);

        lambda = W + w + M;
        e_x = e * cos(w);
        e_y = e * sin(w);
        i_x = sin(i) * cos(W);
        i_y = sin(i) * sin(W);

        QNSOE(k,:) = [a, lambda, e_x, e_y, i_x, i_y];
    end
end