%% Define Rendezvous 
alpha_c = [6771e3; 5.0e-4; deg2rad(51.64); deg2rad(257); deg2rad(45); deg2rad(30)];
delta_alpha_initial = [0; 0.001; deg2rad(0.5); deg2rad(1); deg2rad(1); deg2rad(-5)];
delta_alpha_desired = [0; 0; 0; 0; 0; deg2rad(-0.05)];
alpha_d_0 = alpha_c + delta_alpha_initial;
alpha_d_f = alpha_c + delta_alpha_desired;


QNS_chief = utils.OE2QNSOE(alpha_c);
QNSROE_0 = utils.ROE2QNSROE(alpha_c, delta_alpha_initial);
QNSROE_F = utils.ROE2QNSROE(alpha_c, delta_alpha_desired);


%% Graph Reachable Sets for Max Delta V
mu      = 3.986004418e14;
N_nu    = 200;   % orbit samples
N_phi   = 200;   % direction samples
N_th    = 200;    % polar samples for inclination plane
dv_max = 7.36;

% Compute hulls and sample points
[kpolys, allPts] = computeReachableHulls(alpha_c, mu, N_nu, N_phi, N_th, dv_max);
QNSROE_des = alpha_c(1)*(QNSROE_F - QNSROE_0);

figure(1);
subplot(1,3,1)
plane = 1;
pts   = allPts{plane};       % 2×M array of sampled (Δx,Δy)
idx   = kpolys{plane};       % hull indices into pts
plot(pts(1,:), pts(2,:), '.', 'MarkerSize',2); hold on
plot(pts(1,idx), pts(2,idx), 'r-', 'LineWidth',2);
xlabel('a\Delta\deltaa (m)');
ylabel('a\Delta\delta\lambda (m)');

subplot(1,3,2)
plane = 2;
pts   = allPts{plane};       % 2×M array of sampled (Δx,Δy)
idx   = kpolys{plane};       % hull indices into pts
plot(pts(1,:), pts(2,:), '.', 'MarkerSize',2); hold on
plot(pts(1,idx), pts(2,idx), 'r-', 'LineWidth',2);
xlabel('a\Delta\deltae_x (m)');
ylabel('a\Delta\deltae_y (m)');

subplot(1,3,3)
plane = 3;
pts   = allPts{plane};       % 2×M array of sampled (Δx,Δy)
idx   = kpolys{plane};       % hull indices into pts
plot(pts(1,:), pts(2,:), '.', 'MarkerSize',2); hold on
plot(pts(1,idx), pts(2,idx), 'r-', 'LineWidth',2);
xlabel('a\Delta\deltai_x (m)');
ylabel('a\Delta\deltai_y (m)');
legend('samples', 'hull')


figure(2);
subplot(1,3,1)
plane = 1;
pts   = allPts{plane};       % 2×M array of sampled (Δx,Δy)
idx   = kpolys{plane};       % hull indices into pts
plot(pts(1,:), pts(2,:), '.', 'MarkerSize',2); hold on
plot(pts(1,idx), pts(2,idx), 'r-', 'LineWidth',2);
scatter(QNSROE_des(1), QNSROE_des(2), 'g', 'filled')
xlabel('a\Delta\deltaa (m)');
ylabel('a\Delta\delta\lambda (m)');

subplot(1,3,2)
plane = 2;
pts   = allPts{plane};       % 2×M array of sampled (Δx,Δy)
idx   = kpolys{plane};       % hull indices into pts
plot(pts(1,:), pts(2,:), '.', 'MarkerSize',2); hold on
plot(pts(1,idx), pts(2,idx), 'r-', 'LineWidth',2);
scatter(QNSROE_des(3), QNSROE_des(4), 'g', 'filled')
xlabel('a\Delta\deltae_x (m)');
ylabel('a\Delta\deltae_y (m)');

subplot(1,3,3)
plane = 3;
pts   = allPts{plane};       % 2×M array of sampled (Δx,Δy)
idx   = kpolys{plane};       % hull indices into pts
plot(pts(1,:), pts(2,:), '.', 'MarkerSize',2); hold on
plot(pts(1,idx), pts(2,idx), 'r-', 'LineWidth',2);
scatter(QNSROE_des(5), QNSROE_des(6), 'g', 'filled')
xlabel('a\Delta\deltai_x (m)');
ylabel('a\Delta\deltai_y (m)');
legend('samples', 'hull', 'required')



%% METHOD 1: Calculate optimal nu / min delta v using planes
deltaQNS = QNSROE_F - QNSROE_0;
planes = {[1 2], [3 4], [5 6]}; 
mu   = 3.986004418e14;
N    = 200;  % number of ν samples in dv_min_calc
for p = 1:3
  rows = planes{p};
  dpl       = deltaQNS(rows);                 
  [dv_p, nu_p] = dv_min_calc_planes(rows, dpl, alpha_c, N, mu);
  dv_min(p) = dv_p;
  nu_opt(p) = nu_p;
end
[dv_tot, idx] = max(dv_min);  
nu_star       = nu_opt(idx);    % the anomaly where that “worst” plane is met
fprintf('Method 1 (Planes): Minimum ∆v = %.3f m/s at ν = %.1f°\n', dv_tot, rad2deg(nu_star));


%% METHOD 2: Calculate optimal nu / min delta v using all 6 at once
[v_best, nu_best] = dv_min_calc_global(alpha_c, N, deltaQNS, mu);
fprintf('Method 2 (True optimal): Minimum ∆v = %.3f m/s at ν = %.1f°\n', norm(v_best), rad2deg(nu_best));
disp('RTN vector:'), disp(v_best)


%% Impulse control
% Burn schedule
J2 = 0.001082;
Re = 6378e3;
tau   = 12*3600;
t_burn = [0.025:0.025:0.975]*tau;
t0     = 0;
t_f    = tau;
M      = numel(t_burn);
nu0 = alpha_c(6);
a = alpha_c(1);
e = alpha_c(2);


% Build A
A = zeros(6,3*M);
for j = 1:M
  dt_b    = t_burn(j)-t0;
  dt_to_f = t_f - t_burn(j);

  Phi0_burn   = PHI_J2(alpha_c,mu,J2,Re,dt_b);
  Phiburn_f   = PHI_J2(alpha_c,mu,J2,Re,dt_to_f);
  nu_j   = propagate_nu(nu0,e,a,dt_b,mu, J2, Re, alpha_c(3));  
  Gamma_j = computeGamma(alpha_c,mu,nu_j);

  A(:,3*(j-1)+1:3*j) = Phiburn_f * Gamma_j;
end

% Pseudostate and solve
Phi0_f = PHI_J2(alpha_c,mu,J2,Re,t_f - t0);
b      = QNSROE_F - Phi0_f*QNSROE_0;

dv_stack = pinv(A)*b;            % minimal ‖dv‖₂ solution
t_star   = max(abs(dv_stack));   % largest single‐burn magnitude
disp('Burns:')
for j = 1:M
  dvj = dv_stack(3*(j-1)+1:3*j);
  fprintf(' burn %d → [%.3f  %.3f  %.3f]\n', j, dvj);
end

dvR = dv_stack(1:3:end);    % radial
dvT = dv_stack(2:3:end);    % tangential
dvN = dv_stack(3:3:end);    % normal

burn_mags = sqrt(dvR.^2 + dvT.^2 + dvN.^2);
total_dv  = sum(burn_mags);
fprintf('Maximum impulse = %.3f m/s\n', t_star);
fprintf('Total ∆v = %.3f m/s\n', total_dv);

figure(3)
subplot(3,1,1)
scatter(t_burn/3600, dvR, 'filled')   % time in hours
ylabel('Δv_R (m/s)')
xlabel('Time (h)')
grid on
title('Radial Impulses')

subplot(3,1,2)
scatter(t_burn/3600, dvT, 'filled')
ylabel('Δv_T (m/s)')
xlabel('Time (h)')
grid on
title('Tangential Impulses')

subplot(3,1,3)
scatter(t_burn/3600, dvN, 'filled')
ylabel('Δv_N (m/s)')
xlabel('Time (h)')
grid on
title('Normal Impulses')



%% Ground‐Truth Simulation, with and without control error

% 1% relative noise on each burn
err_sigma = 0.05;

% preallocate containers
T_nominal             = [];
T_perturbed           = [];
QNSROE_nominal        = [];
QNSROE_perturbed      = [];

for caseIdx = 1:2
    % pick exact or noisy burns
    if caseIdx == 1
        dv_use = dv_stack;  
    else
        dv_use = dv_stack .* (1 + err_sigma * randn(size(dv_stack)));
    end

    % initialize state
    dep_ECI = utils.OE2ECI(alpha_d_0, mu);
    chi_ECI = utils.OE2ECI(alpha_c,   mu);
    dep_RTN   = utils.ECI2RTN(chi_ECI', dep_ECI');
    state           = [chi_ECI; dep_RTN];
    t_prev          = 0;

    % storage for this run
    T_all = [];
    X_all = [];

    % loop burns
    for j = 1:M
        % coast to burn j
        tspan = [t_prev, t_burn(j)];
        [T, X] = ode113(@(t,x) propagators.NonlinearPropagator.nonlinear_state_dot(t,x, mu, Re, J2), tspan, state);
        T_all = [T_all; T(1:end-1)];
        X_all = [X_all; X(1:end-1,:)];
        
        % apply burn j
        state = X(end,:)';
        dvj   = dv_use(3*(j-1)+1 : 3*j);
        state(10:12) = state(10:12) + dvj;
        t_prev = t_burn(j);
    end

    % final coast
    [T, X] = ode113(@(t,x) propagators.NonlinearPropagator.nonlinear_state_dot(t,x, mu, Re, J2), [t_prev, tau], state);
    T_all   = [T_all;   T];
    X_all   = [X_all;   X];

    % reconstruct QNSROE history
    N = size(T_all,1);
    Qhist = zeros(N,6);
    for k = 1:N
        r_c   = X_all(k,1:3)';   v_c   = X_all(k,4:6)';
        r_rel = X_all(k,7:9)';   v_rel = X_all(k,10:12)';
        [r_d, v_d] = utils.RTN2ECI(r_c, v_c, r_rel, v_rel);

        chief_OE = utils.ECI2OE([r_c; v_c], mu);
        deputy_OE = utils.ECI2OE([r_d; v_d], mu);

        dSROE           = deputy_OE - chief_OE;
        Qhist(k,:)      = utils.ROE2QNSROE(chief_OE, dSROE)';
    end

    % store results
    if caseIdx == 1
        T_nominal        = T_all;
        QNSROE_nominal   = Qhist;
    else
        T_perturbed      = T_all;
        QNSROE_perturbed = Qhist;
    end
end

dl_nom = unwrap(QNSROE_nominal(:,2));
dl_err = unwrap(QNSROE_perturbed(:,2));


figure(4);
% δa
subplot(3,1,1)
plot(T_nominal/3600,    QNSROE_nominal(:,1),   'b-','LineWidth',1.5); hold on
plot(T_perturbed/3600,  QNSROE_perturbed(:,1), 'r--','LineWidth',1.5);
ylabel('\delta a');    grid on
legend('Nominal','With Error','Location','best')

% δλ
subplot(3,1,2)
plot(T_nominal/3600,    dl_nom,   'b-','LineWidth',1.5); hold on
plot(T_perturbed/3600,  dl_err, 'r--','LineWidth',1.5);
ylabel('\delta \lambda');    grid on

% e_x
subplot(3,1,3)
plot(T_nominal/3600,    QNSROE_nominal(:,3),   'b-','LineWidth',1.5); hold on
plot(T_perturbed/3600,  QNSROE_perturbed(:,3), 'r--','LineWidth',1.5);
ylabel('e_x');    grid on

figure(5)
% e_y
subplot(3,1,1)
plot(T_nominal/3600,    QNSROE_nominal(:,4),   'b-','LineWidth',1.5); hold on
plot(T_perturbed/3600,  QNSROE_perturbed(:,4), 'r--','LineWidth',1.5);
ylabel('e_y');    grid on

% δi_x
subplot(3,1,2)
plot(T_nominal/3600,    QNSROE_nominal(:,5),   'b-','LineWidth',1.5); hold on
plot(T_perturbed/3600,  QNSROE_perturbed(:,5), 'r--','LineWidth',1.5);
ylabel('\delta i_x');    grid on

% δi_y
subplot(3,1,3)
plot(T_nominal/3600,    QNSROE_nominal(:,6),   'b-','LineWidth',1.5); hold on
plot(T_perturbed/3600,  QNSROE_perturbed(:,6), 'r--','LineWidth',1.5);
ylabel('\delta i_y');
xlabel('Time [h]');    grid on





function nu = propagate_nu(nu0,e,a,dt,mu,J2,Re,inc)
  % Unperturbed Kepler propagation to get keplerian nu_kep
  E0    = 2*atan( sqrt((1-e)/(1+e)) * tan(nu0/2) );
  M0    = E0 - e*sin(E0);
  n     = sqrt(mu/a^3);
  M     = mod(M0 + n*dt, 2*pi);
  E     = M;
  for it=1:20
    f  = E - e*sin(E) - M;
    fp = 1 - e*cos(E);
    dE = f/fp;
    E  = E - dE;
    if abs(dE)<1e-12, break; end
  end
  nu_kep = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );

  % J2 secular perigee drift rate
  eta   = sqrt(1-e^2);
  kappa = (3*J2*Re^2*sqrt(mu))/(4*a^(7/2)*eta^4);
  Q     = 5*cos(inc)^2 - 1;
  omega_dot = kappa * Q;

  % Combine Keplerian nu with perigee drift
  nu = mod(nu_kep + omega_dot*dt, 2*pi);
end





function [v_best, nu_best] = dv_min_calc_global(alpha_c, N, deltaQNS, mu)
a = alpha_c(1);
e = alpha_c(2);
omega = alpha_c(5);

nus    = linspace(0,2*pi,N);
dv_best = inf;
nu_best = 0;
for k = 1:N
  Gk    = computeGamma(alpha_c,mu,nus(k));   % 6×3
  vk    = pinv(Gk)*deltaQNS;                   % 3×1 least-norm solution
  dv_k  = norm(vk);
  if dv_k < dv_best
    dv_best = dv_k;
    nu_best = nus(k);
    v_best  = vk;
  end
end
end



function [dv_min, nu_opt] = dv_min_calc_planes(plane_rows, d, alpha_c, N, mu)
a = alpha_c(1);
e = alpha_c(2);
omega = alpha_c(5);

nus    = linspace(0,2*pi,N);
dv_min = inf;
nu_opt = 0;
for k=1:N
    G = computeGamma(alpha_c,mu,nus(k)); % 6×3

    % Decide which columns to use
    if isequal(plane_rows, [5 6])
        % Inclination plane: only normal impulse (3rd col)
        P = G(plane_rows,3);   % 2×1
        % d = P * v_N
        % Solve for v_N (scalar)
        v = P \ d; % Least-squares solution (can use pinv(P)*d if you want)
        test_norm = norm(v); % This is a scalar
    else
        % In-plane: radial/tangential only
        P = G(plane_rows,1:2); % 2×2
        v = pinv(P)*d;         % 2×1
        test_norm = norm(v);   % Vector norm
    end

    if test_norm < dv_min
        dv_min = test_norm;
        nu_opt = nus(k);
    end
end
end


function [kpolys, allPts] = computeReachableHulls(alpha_c, mu, N_nu, N_phi, N_th, dv_max)
% Returns convex‐hull indices and sampled points in each QNSROE plane.

a     = alpha_c(1);
e     = alpha_c(2);
omega = alpha_c(5);

nus    = linspace(0,2*pi, N_nu);
phis   = linspace(0,2*pi, N_phi);
thetas = linspace(0,pi,   N_th);

planes = {[1 2], [3 4], [5 6]};
kpolys = cell(1,3);
allPts  = cell(1,3);

for p = 1:3
  rows = planes{p};
  pts  = zeros(2,0);

  for k = 1:N_nu
    G = computeGamma(alpha_c,mu,nus(k));  % 6×3 control‐input matrix

    if p < 3
      % In‐plane burns: only radial/tangential
      P2 = G(rows,1:2);
      for m = 1:N_phi
        u2  = dv_max*[cos(phis(m)); sin(phis(m))];  % 2×1
        pts(:,end+1) = P2*u2;
      end

    else
      % Inclination plane: only normal burns
      P2 = G(rows,3);    % 2×1
      for ti = 1:N_th
        uN = dv_max*cos(thetas(ti));      
        pts(:,end+1) = P2*uN;
      end
    end
  end

  % store points & convex hull
  allPts{p}  = a*pts;
  kpolys{p}  = convhull(pts(1,:), pts(2,:));
  
end
end


function Gamma = computeGamma(alpha_c, mu, nu_k)
a     = alpha_c(1);
e     = alpha_c(2);
omega = alpha_c(5);

n   = sqrt(mu / a^3); % mean motion
eta = sqrt(1 - e^2); % η = √(1−e²)
theta_k = nu_k + omega; 

cc = 1 + e*cos(nu_k);

% build Gamma
Gamma = (1/(n*a)) * [
    (2/eta)*e*sin(nu_k), (2/eta)*cc, 0;
   -(2*eta^2)/cc, 0, 0;
    eta*sin(nu_k), eta*(e + ((2+e*cos(nu_k)).*cos(nu_k))./cc), 0;
   -(eta/e)*cos(nu_k), (eta/e)*(((2+e*cos(nu_k)).*sin(nu_k))./cc), 0;
    0, 0, eta*(cos(theta_k)./cc);
    0, 0, eta*sin(theta_k)./cc];
end


function Phi = PHI_J2(alpha_c, mu, J2, Re, tau)
  a = alpha_c(1);
  e = alpha_c(2);
  inc = alpha_c(3);
  omega = alpha_c(5);

  % mean motion
  n = sqrt(mu/(a^3));

  % eccentricity-dependent terms
  eta   = sqrt(1 - e^2);
  kappa = (3*J2*Re^2*sqrt(mu))/(4*a^(7/2)*eta^4);
  E     = 1 + eta;
  F     = 4 + 3*eta;
  G     = 1/eta^2;

  % inclination-dependent terms
  P = 3*cos(inc)^2 - 1;
  Q = 5*cos(inc)^2 - 1;
  S = sin(2*inc);
  T = sin(inc)^2;

  % exi, exf, eyi, eyf
  omega_dot = kappa*Q;
  omega_f = omega + omega_dot*tau;
  exi = e*cos(omega);
  eyi = e*sin(omega);
  exf = e*cos(omega_f);
  eyf = e*sin(omega_f);


  % assemble the STM
  Phi = [
    1,                                  0,                  0,                    0,                   0,                0;
   -7*kappa*eta*P*tau-1.5*n*tau,        1,     7*kappa*exi*P*tau/eta,     7*kappa*eyi*P*tau/eta,    -7*kappa*eta*S*tau,    0;
    3.5*kappa*eyf*Q*tau,                 0, cos(omega_dot*tau)-4*kappa*exi*eyf*G*Q*tau, -sin(omega_dot*tau)-4*kappa*eyi*eyf*G*Q*tau,  5*kappa*eyf*S*tau, 0;
   -3.5*kappa*exf*Q*tau,                 0, sin(omega_dot*tau)+4*kappa*exi*exf*G*Q*tau,  cos(omega_dot*tau)+4*kappa*eyi*exf*G*Q*tau, -5*kappa*exf*S*tau, 0;
    0,                                  0,                  0,                    0,                   1,                0;
    3.5*kappa*S*tau,                    0, -4*kappa*exi*G*S*tau,   -4*kappa*eyi*G*S*tau,  2*kappa*T*tau,      1
  ];
end