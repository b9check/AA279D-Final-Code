%% Problem Set 2

%% 1b
% Constants
mu = 3.986004418e14;  
Re = 6378137;          
J2 = 1.08262668e-3;   

% Chief orbital elements
chief_OE = [6771e3; 0.0005; deg2rad(51.64); deg2rad(257); deg2rad(0); deg2rad(30)];
deputy_OE = [6771e3; 0.0015; deg2rad(52.14); deg2rad(257.5); deg2rad(0.5); deg2rad(25)];

% Convert OE to ECI and RTN state vectors
chief_ECI  = utils.OE2ECI(chief_OE, mu);
deputy_ECI = utils.OE2ECI(deputy_OE, mu);
deputy_RTN = utils.ECI2RTN(chief_ECI', deputy_ECI');

% Initialize combined state: [chief_ECI; deputy_RTN]
state0 = zeros(12, 1);
state0(1:6, :) = chief_ECI;
state0(7:12, :) = deputy_RTN;

% Propagate nonlinear dynamics
T = 2*pi*sqrt(chief_OE(1)^3 / mu); % Orbital period
tspan = 0:1:T*10;
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
odefun = @(t, state) propagators.NonlinearPropagator.nonlinear_state_dot(t, state, mu, Re, J2);
[t_out1, x_out1] = ode113(odefun, tspan, state0, options);

% PLOT
R1 = x_out1(:, 7);
T1 = x_out1(:, 8);
N1 = x_out1(:, 9);
vR1 = x_out1(:, 10);
vT1 = x_out1(:, 11);
vN1 = x_out1(:, 12);

% FIGURE 1 — RTN Position Components
figure(1)
subplot(2,2,1)
plot(R1, T1)
xlabel('R (m)')
ylabel('T (m)')
grid on

subplot(2,2,2)
plot(R1, N1)
xlabel('R (m)')
ylabel('N (m)')
grid on

subplot(2,2,3)
plot(N1, T1)
xlabel('N (m)')
ylabel('T (m)')
grid on

subplot(2,2,4)
plot3(R1, T1, N1)
xlabel('R (m)')
ylabel('T (m)')
zlabel('N (m)')
grid on

% FIGURE 2 — RTN Velocity Components
figure(2)
subplot(2,2,1)
plot(vR1, vT1)
xlabel('R_{dot} (m/s)')
ylabel('T_{dot} (m/s)')
grid on

subplot(2,2,2)
plot(vR1, vN1)
xlabel('R_{dot} (m/s)')
ylabel('N_{dot} (m/s)')
grid on

subplot(2,2,3)
plot(vN1, vT1)
xlabel('N_{dot} (m/s)')
ylabel('T_{dot} (m/s)')
grid on

subplot(2,2,4)
plot3(vR1, vT1, vN1)
xlabel('R_{dot} (m/s)')
ylabel('T_{dot} (m/s)')
zlabel('N_{dot} (m/s)')
grid on

% FIGURE 3 — RTN Position Components vs Orbital Periods
figure(3)
subplot(3,1,1)
plot(t_out1/T, R1)
ylabel('R (m)')
grid on

subplot(3,1,2)
plot(t_out1/T, T1)
ylabel('T (m)')
grid on

subplot(3,1,3)
plot(t_out1/T, N1)
xlabel('Orbital Periods')
ylabel('N (m)')
grid on


% FIGURE 4 — RTN Velocity Components vs Time
figure(4)
subplot(3,1,1)
plot(t_out1/T, vR1)
ylabel('v_R (m/s)')
grid on

subplot(3,1,2)
plot(t_out1/T, vT1)
ylabel('v_T (m/s)')
grid on

subplot(3,1,3)
plot(t_out1/T, vN1)
xlabel('Orbital Periods')
ylabel('v_N (m/s)')
grid on



%% 1c
% Propagate chief
T = 2*pi*sqrt(chief_OE(1)^3 / mu); % Period
tspan = 0:1:10*T;  
state0 = chief_ECI;

% Call ode113 to propagate chief without J2 perturbations
perturbated = false;
odefun = @(t, state) propagators.FodePropagator.getStatedot(t, state, mu, Re, J2, perturbated);
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
[t_out2, x_out2] = ode113(odefun, tspan, state0, opts);

% Propagate deputy
state0 = utils.OE2ECI(deputy_OE, mu);
[t_out3, x_out3] = ode113(odefun, tspan, state0, opts);
deputy_ECI = utils.OE2ECI(deputy_OE, mu);

RTNs = zeros(length(tspan), 6);
for i = 1:length(tspan)
    RTN = utils.ECI2RTN(x_out2(i,:), x_out3(i,:));
    RTNs(i, :) = RTN;
end


% PLOT
R2 = RTNs(:, 1);
T2 = RTNs(:, 2);
N2 = RTNs(:, 3);
vR2 = RTNs(:, 4);
vT2 = RTNs(:, 5);
vN2 = RTNs(:, 6);

% FIGURE 5 — RTN Position Components
figure(5)
subplot(2,2,1)
plot(R2, T2, 'r')
xlabel('R (m)')
ylabel('T (m)')
grid on

subplot(2,2,2)
plot(R2, N2, 'r')
xlabel('R (m)')
ylabel('N (m)')
grid on

subplot(2,2,3)
plot(N2, T2, 'r')
xlabel('N (m)')
ylabel('T (m)')
grid on

subplot(2,2,4)
plot3(R2, T2, N2, 'r')
xlabel('R (m)')
ylabel('T (m)')
zlabel('N (m)')
grid on

% FIGURE 6 — RTN Velocity Components
figure(6)
subplot(2,2,1)
plot(vR2, vT2, 'r')
xlabel('R_{dot} (m/s)')
ylabel('T_{dot} (m/s)')
grid on

subplot(2,2,2)
plot(vR2, vN2, 'r')
xlabel('R_{dot} (m/s)')
ylabel('N_{dot} (m/s)')
grid on

subplot(2,2,3)
plot(vN2, vT2, 'r')
xlabel('N_{dot} (m/s)')
ylabel('T_{dot} (m/s)')
grid on

subplot(2,2,4)
plot3(vR2, vT2, vN2, 'r')
xlabel('R_{dot} (m/s)')
ylabel('T_{dot} (m/s)')
zlabel('N_{dot} (m/s)')
grid on


% FIGURE 7 — RTN Position Components vs Orbital Periods
figure(7)
subplot(3,1,1)
plot(t_out2/T, R2)
ylabel('R (m)')
grid on

subplot(3,1,2)
plot(t_out2/T, T2)
ylabel('T (m)')
grid on

subplot(3,1,3)
plot(t_out2/T, N2)
xlabel('Orbital Periods')
ylabel('N (m)')
grid on


% FIGURE 8 — RTN Velocity Components vs Time
figure(8)
subplot(3,1,1)
plot(t_out2/T, vR2)
ylabel('v_R (m/s)')
grid on

subplot(3,1,2)
plot(t_out2/T, vT2)
ylabel('v_T (m/s)')
grid on

subplot(3,1,3)
plot(t_out2/T, vN2)
xlabel('Orbital Periods')
ylabel('v_N (m/s)')
grid on


%% 1d Comparison

% FIGURE 9 — RTN Position Components
figure(9)
subplot(2,2,1)
plot(R1, T1, 'b', R2, T2, 'r--')
xlabel('R (m)')
ylabel('T (m)')
legend('R1 vs T1', 'R2 vs T2')
grid on

subplot(2,2,2)
plot(R1, N1, 'b', R2, N2, 'r--')
xlabel('R (m)')
ylabel('N (m)')
legend('R1 vs N1', 'R2 vs N2')
grid on

subplot(2,2,3)
plot(N1, T1, 'b', N2, T2, 'r--')
xlabel('N (m)')
ylabel('T (m)')
legend('N1 vs T1', 'N2 vs T2')
grid on

subplot(2,2,4)
plot3(R1, T1, N1, 'b', R2, T2, N2, 'r--')
xlabel('R (m)')
ylabel('T (m)')
zlabel('N (m)')
legend('RTN1', 'RTN2')
grid on

% FIGURE 10 — RTN Velocity Components
figure(10)
subplot(2,2,1)
plot(vR1, vT1, 'b', vR2, vT2, 'r--')
xlabel('R_{dot} (m/s)')
ylabel('T_{dot} (m/s)')
legend('vR1 vs vT1', 'vR2 vs vT2')
grid on

subplot(2,2,2)
plot(vR1, vN1, 'b', vR2, vN2, 'r--')
xlabel('R_{dot} (m/s)')
ylabel('N_{dot} (m/s)')
legend('vR1 vs vN1', 'vR2 vs vN2')
grid on

subplot(2,2,3)
plot(vN1, vT1, 'b', vN2, vT2, 'r--')
xlabel('N_{dot} (m/s)')
ylabel('T_{dot} (m/s)')
legend('vN1 vs vT1', 'vN2 vs vT2')
grid on

subplot(2,2,4)
plot3(vR1, vT1, vN1, 'b', vR2, vT2, vN2, 'r--')
xlabel('R_{dot} (m/s)')
ylabel('T_{dot} (m/s)')
zlabel('N_{dot} (m/s)')
legend('vRTN1', 'vRTN2')
grid on


% FIGURE 11 — RTN Position Components vs Orbital Periods
figure(11)
subplot(3,1,1)
plot(t_out1/T, R1, 'b', t_out2/T, R2, 'r--')
ylabel('R (m)')
legend('R1', 'R2')
grid on

subplot(3,1,2)
plot(t_out1/T, T1, 'b', t_out2/T, T2, 'r--')
ylabel('T (m)')
legend('T1', 'T2')
grid on

subplot(3,1,3)
plot(t_out1/T, N1, 'b', t_out2/T, N2, 'r--')
xlabel('Orbital Periods')
ylabel('N (m)')
legend('N1', 'N2')
grid on

% FIGURE 12 — RTN Velocity Components vs Orbital Periods
figure(12)
subplot(3,1,1)
plot(t_out1/T, vR1, 'b', t_out2/T, vR2, 'r--')
ylabel('v_R (m/s)')
legend('vR1', 'vR2')
grid on

subplot(3,1,2)
plot(t_out1/T, vT1, 'b', t_out2/T, vT2, 'r--')
ylabel('v_T (m/s)')
legend('vT1', 'vT2')
grid on

subplot(3,1,3)
plot(t_out1/T, vN1, 'b', t_out2/T, vN2, 'r--')
xlabel('Orbital Periods')
ylabel('v_N (m/s)')
legend('vN1', 'vN2')
grid on


%% 1d — New set of ICs with ∆a ≠ 0 (drift-inducing)

% New deputy with different semi-major axis
deputy_OE_drift = [6781e3; 0.0015; deg2rad(52.14); deg2rad(257.5); deg2rad(0.5); deg2rad(25)];

% Convert OE to ECI and RTN
deputy_ECI_drift = utils.OE2ECI(deputy_OE_drift, mu);
deputy_RTN_drift = utils.ECI2RTN(chief_ECI', deputy_ECI_drift');

% Initial state [chief; relative RTN deputy]
state0_drift = zeros(12,1);
state0_drift(1:6) = chief_ECI;
state0_drift(7:12) = deputy_RTN_drift;

% Reuse tspan and options
[t_drift_rel, x_drift_rel] = ode113(@(t, state) propagators.NonlinearPropagator.nonlinear_state_dot(t, state, mu, Re, J2*0), tspan, state0_drift, options);

% Extract
R1d = x_drift_rel(:, 7);
T1d = x_drift_rel(:, 8);
N1d = x_drift_rel(:, 9);
vR1d = x_drift_rel(:, 10);
vT1d = x_drift_rel(:, 11);
vN1d = x_drift_rel(:, 12);

% Absolute motion method
state0 = utils.OE2ECI(deputy_OE_drift, mu);
[t_drift_abs, x_drift_abs] = ode113(@(t, state) propagators.FodePropagator.getStatedot(t, state, mu, Re, J2*0, false), tspan, state0, opts);

% RTN relative transformation
RTNs_drift = zeros(length(tspan),6);
for i = 1:length(tspan)
    RTNs_drift(i,:) = utils.ECI2RTN(x_out2(i,:), x_drift_abs(i,:));
end

R2d = RTNs_drift(:, 1); T2d = RTNs_drift(:, 2); N2d = RTNs_drift(:, 3);
vR2d = RTNs_drift(:, 4); vT2d = RTNs_drift(:, 5); vN2d = RTNs_drift(:, 6);

%% Comparison — Drift Case
% FIGURE 13 — RTN Position
figure(13)
subplot(2,2,1); plot(R1d,T1d,'b', R2d,T2d,'r--'); xlabel('R (m)'); ylabel('T (m)'); legend('RTN','FODE'); grid on
subplot(2,2,2); plot(R1d,N1d,'b', R2d,N2d,'r--'); xlabel('R (m)'); ylabel('N (m)'); grid on
subplot(2,2,3); plot(N1d,T1d,'b', N2d,T2d,'r--'); xlabel('N (m)'); ylabel('T (m)'); grid on
subplot(2,2,4); plot3(R1d,T1d,N1d,'b', R2d,T2d,N2d,'r--'); xlabel('R (m)'); ylabel('T (m)'); zlabel('N (m)'); grid on

% FIGURE 14 — RTN Velocity
figure(14)
subplot(2,2,1); plot(vR1d,vT1d,'b', vR2d,vT2d,'r--'); xlabel('v_R (m/s)'); ylabel('v_T (m/s)'); legend('RTN','FODE'); grid on
subplot(2,2,2); plot(vR1d,vN1d,'b', vR2d,vN2d,'r--'); xlabel('v_R (m/s)'); ylabel('v_N (m/s)'); grid on
subplot(2,2,3); plot(vN1d,vT1d,'b', vN2d,vT2d,'r--'); xlabel('v_N (m/s)'); ylabel('v_T (m/s)'); grid on
subplot(2,2,4); plot3(vR1d,vT1d,vN1d,'b', vR2d,vT2d,vN2d,'r--'); xlabel('v_R (m/s)'); ylabel('v_T (m/s)'); zlabel('v_N (m/s)'); grid on

% FIGURE 15 — RTN Position vs Time
figure(15)
subplot(3,1,1); plot(t_drift_rel/T, R1d,'b', t_drift_rel/T, R2d,'r--'); ylabel('R (m)'); legend('RTN','FODE'); grid on
subplot(3,1,2); plot(t_drift_rel/T, T1d,'b', t_drift_rel/T, T2d,'r--'); ylabel('T (m)'); grid on
subplot(3,1,3); plot(t_drift_rel/T, N1d,'b', t_drift_rel/T, N2d,'r--'); ylabel('N (m)'); xlabel('Orbital Periods'); grid on

% FIGURE 16 — RTN Velocity vs Time
figure(16)
subplot(3,1,1); plot(t_drift_rel/T, vR1d,'b', t_drift_rel/T, vR2d,'r--'); ylabel('v_R (m/s)'); legend('RTN','FODE'); grid on
subplot(3,1,2); plot(t_drift_rel/T, vT1d,'b', t_drift_rel/T, vT2d,'r--'); ylabel('v_T (m/s)'); grid on
subplot(3,1,3); plot(t_drift_rel/T, vN1d,'b', t_drift_rel/T, vN2d,'r--'); ylabel('v_N (m/s)'); xlabel('Orbital Periods'); grid on


%% 1e: Impulse calculation
a = chief_OE(1);
v_end = norm(x_out2(end, 4:6));
da = a - deputy_OE_drift(1);
dvT = (mu*da) / (2*a^2*v_end);

disp(v_end)
disp(dvT)


%% 1f: Simulate
chief_initial = x_drift_rel(end, 1:6);
deputy_initial = x_drift_rel(end, 7:12) + [0, 0, 0, 0, dvT, 0];
state0 = [chief_initial'; deputy_initial']';

% Propagate nonlinear dynamics
T = 2*pi*sqrt(chief_OE(1)^3 / mu); % Orbital period
tspan = 0:1:T*10;
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
odefun = @(t, state) propagators.NonlinearPropagator.nonlinear_state_dot(t, state, mu, Re, J2);
[t_out3, x_out3] = ode113(odefun, tspan, state0, options);

RTN_full = [x_drift_rel(:, 7:12); x_out3(2:end, 7:12)];
t_full = [t_drift_rel; t_out3(2:end) + t_drift_rel(end)];
t_orb_full = t_full / T;

% Split position and velocity
R_full  = RTN_full(:, 1);
T_full  = RTN_full(:, 2);
N_full  = RTN_full(:, 3);
vR_full = RTN_full(:, 4);
vT_full = RTN_full(:, 5);
vN_full = RTN_full(:, 6);

% Figure 17
figure(17)
subplot(3,1,1)
plot(t_orb_full, R_full)
ylabel('R (m)')
grid on

subplot(3,1,2)
plot(t_orb_full, T_full)
ylabel('T (m)')
grid on

subplot(3,1,3)
plot(t_orb_full, N_full)
ylabel('N (m)')
xlabel('Orbital Periods')
grid on

% Figure 18
figure(18)
subplot(3,1,1)
plot(t_orb_full, vR_full)
ylabel('v_R (m/s)')
grid on

subplot(3,1,2)
plot(t_orb_full, vT_full)
ylabel('v_T (m/s)')
grid on

subplot(3,1,3)
plot(t_orb_full, vN_full)
ylabel('v_N (m/s)')
xlabel('Orbital Periods')
grid on

