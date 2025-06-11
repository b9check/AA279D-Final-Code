%% ps4_main_regressionMean.m
% Propagate osculating J2, compute QNS & rQNS, then overlay best‐fit slopes

close all; clear;

%% 1. Setup
mu = 3.986004418e14;
Re = 6378137;
J2 = 1.08262668e-3;

% ICs (a)
chief_OE  = [6771e3; 0.0005; deg2rad(51.64); deg2rad(257); 0; deg2rad(30)];
deputy_OE = [6771e3; 0.0015; deg2rad(52.14); deg2rad(257.5); deg2rad(0.5); deg2rad(25)];
% ICs (b)
rQNS_b    = [0,100,50,100,30,200]'/chief_OE(1);
deputy2_OE= utils.rQNSOE2OE(chief_OE, rQNS_b);

% part g
qnsroes0 = utils.OE2rQNSOE(chief_OE, deputy_OE);
rQNS_b = qnsroes0 - [0 0 0 0 qnsroes0(5) 0];
deputy2_OE = utils.rQNSOE2OE(chief_OE, rQNS_b);

T     = 2*pi*sqrt(chief_OE(1)^3/mu);
tspan = (0:1:15*T)';

ode_noJ2 = @(t,s) propagators.KeplerianPropagator.getStatedot(t,s,mu,Re,J2,false);
ode_J2   = @(t,s) propagators.KeplerianPropagator.getStatedot(t,s,mu,Re,J2,true);
opts     = odeset('RelTol',1e-12,'AbsTol',1e-14);

%% 2. Propagate all four osculating cases
[tc1, xc1] = ode113(ode_noJ2, tspan, chief_OE,  opts);
[td1, xd1] = ode113(ode_noJ2, tspan, deputy_OE, opts);
[tc2, xc2] = ode113(ode_J2,   tspan, chief_OE,  opts);
[td2, xd2] = ode113(ode_J2,   tspan, deputy_OE, opts);
[tc3, xc3] = ode113(ode_noJ2, tspan, chief_OE,  opts);
[td3, xd3] = ode113(ode_noJ2, tspan, deputy2_OE,opts);
[tc4, xc4] = ode113(ode_J2,   tspan, chief_OE,  opts);
[td4, xd4] = ode113(ode_J2,   tspan, deputy2_OE,opts);
cases = { {xc1,xd1}, {xc2,xd2}, {xc3,xd3}, {xc4,xd4} };

%% 3. Build QNSOE & rQNSOE for each case
N=length(tspan);
QNS  = zeros(N,6,4);
rQNS = zeros(N,6,4);
for c=1:4
    xc = cases{c}{1};
    xd = cases{c}{2};
    for k=1:N
        oe_c = utils.ECI2OE(utils.OE2ECI(xc(k,:)',mu),mu);
        oe_d = utils.ECI2OE(utils.OE2ECI(xd(k,:)',mu),mu);
        QNS(k,:,c)  = utils.OE2QNSOE( oe_d );
        rQNS(k,:,c) = utils.OE2rQNSOE( oe_c, oe_d );
    end
    QNS(:,2,c)  = unwrap(QNS(:,2,c));
    rQNS(:,2,c) = unwrap(rQNS(:,2,c));
end

%% 4. Plot QNSOE and rQNSOE with regression mean

labels = {'\delta a','\delta \lambda','\delta e_x','\delta e_y','\delta i_x','\delta i_y'};
rlabels = {'a \delta a','a \delta \lambda','a \delta e_x','a \delta e_y','a \delta i_x','a \delta i_y'};
for c = 1:4
    t_orb = tspan / T;

    % QNSOE
    fig = 2*(c-1) + 1;
    figure(fig); clf; set(gcf,'Color','w');
    for j = 1:6
        y = QNS(:,j,c);
        p = polyfit(t_orb, y, 1);
        yfit = p(1)*t_orb + p(2);
        subplot(6,1,j)
        plot(t_orb, y, 'b', t_orb, yfit, 'r--','LineWidth',1.2);
        if j == 1
            ylabel([labels{j} ' (m)']);
        elseif j == 2 || j == 5 || j == 6
            ylabel([labels{j} ' (rad)']);
        else
            ylabel(labels{j});
        end
        grid on
        if j < 6
            set(gca, 'XTickLabel', []);
        else
            xlabel('Orbital Periods');
            legend('Osculating', 'Mean', 'Location', 'best');
        end
    end

    % rQNSOE (converted to meters)
    fig = fig + 1;
    figure(fig); clf; set(gcf,'Color','w');
    a0 = chief_OE(1);
    for j = 1:6
        y = a0 * rQNS(:,j,c);
        p = polyfit(t_orb, y, 1);
        yfit = p(1)*t_orb + p(2);
        subplot(6,1,j)
        plot(t_orb, y, 'b', t_orb, yfit, 'r--','LineWidth',1.2);
        ylabel([rlabels{j} ' (m)']); grid on
        if j < 6
            set(gca, 'XTickLabel', []);
        else
            xlabel('Orbital Periods');
            legend('Osculating', 'Mean', 'Location', 'best');
        end
    end
end

%% 5. RTN trajectories for case2 (IC a, with J2)
RTN_unp = zeros(N,3);
RTN_J2  = zeros(N,3);
for k=1:N
    c1 = utils.ECI2RTN( utils.OE2ECI(xc1(k,:)',mu), utils.OE2ECI(xd1(k,:)',mu) );
    RTN_unp(k,:) = c1(1:3)';
    c2 = utils.ECI2RTN( utils.OE2ECI(xc2(k,:)',mu), utils.OE2ECI(xd2(k,:)',mu) );
    RTN_J2(k,:)  = c2(1:3)';
end
figure(9); clf; set(gcf,'Color','w');

% 3D plot (T vs R vs N)
subplot(221)
plot3(RTN_unp(:,2),RTN_unp(:,1),RTN_unp(:,3),'r--',RTN_J2(:,2),RTN_J2(:,1),RTN_J2(:,3),'b-');
grid on; legend('NoJ2','J2');
xlabel('T (m)'); ylabel('R (m)'); zlabel('N (m)');

% T vs R
subplot(222)
plot(RTN_unp(:,2),RTN_unp(:,1),'r--',RTN_J2(:,2),RTN_J2(:,1),'b-');
axis equal; grid on;
xlabel('T (m)'); ylabel('R (m)');

% N vs R
subplot(223)
plot(RTN_unp(:,3),RTN_unp(:,1),'r--',RTN_J2(:,3),RTN_J2(:,1),'b-');
axis equal; grid on;
xlabel('N (m)'); ylabel('R (m)');

% T vs N
subplot(224)
plot(RTN_unp(:,2),RTN_unp(:,3),'r--',RTN_J2(:,2),RTN_J2(:,3),'b-');
axis equal; grid on;
xlabel('T (m)'); ylabel('N (m)');


%% 6. Vector plots for case 1 and case 2
a0 = chief_OE(1);

%% Figure 10 — Case 1 (no J2)
figure(10); clf; set(gcf,'Color','w');

% a) delta e vector
ex = a0 * rQNS(:,3,1);
ey = a0 * rQNS(:,4,1);
subplot(3,1,1);
plot(ex, ey, 'b.'); hold on;
p = polyfit(ex, ey, 1);
xx = linspace(min(ex), max(ex), 2);
plot(xx, p(1)*xx + p(2), 'r--', 'LineWidth', 1.5);
axis equal; grid on;
xlabel('a \delta e_x [m]'); ylabel('a \delta e_y [m]');
legend('Osc','Mean','Location','best');

% b) delta i vector
ix = a0 * rQNS(:,5,1);
iy = a0 * rQNS(:,6,1);
subplot(3,1,2);
plot(ix, iy, 'b.'); hold on;
p = polyfit(ix, iy, 1);
xx = linspace(min(ix), max(ix), 2);
plot(xx, p(1)*xx + p(2), 'r--', 'LineWidth', 1.5);
axis equal; grid on;
xlabel('a \delta i_x [m]'); ylabel('a \delta i_y [m]');
legend('Osc','Mean','Location','best');

% c) delta lambda vs delta a
x_lambda = a0 * rQNS(:,2,1);
y_da = a0 * rQNS(:,1,1);
subplot(3,1,3);
plot(x_lambda, y_da, 'b.'); hold on;
p = polyfit(x_lambda, y_da, 1);
xx = linspace(min(x_lambda), max(x_lambda), 2);
plot(xx, p(1)*xx + p(2), 'r--', 'LineWidth', 1.5);
axis equal; grid on;
xlabel('a \delta \lambda [m]'); ylabel('a \delta a [m]');
legend('Osc','Mean','Location','best');

%% Figure 11 — Case 2 (with J2)
figure(11); clf; set(gcf,'Color','w');

% a) delta e vector
ex = a0 * rQNS(:,3,2);
ey = a0 * rQNS(:,4,2);
subplot(3,1,1);
plot(ex, ey, 'b.'); hold on;
p = polyfit(ex, ey, 1);
xx = linspace(min(ex), max(ex), 2);
plot(xx, p(1)*xx + p(2), 'r--', 'LineWidth', 1.5);
axis equal; grid on;
xlabel('a \delta e_x [m]'); ylabel('a \delta e_y [m]');
legend('Osc','Mean','Location','best');

% b) delta i vector
ix = a0 * rQNS(:,5,2);
iy = a0 * rQNS(:,6,2);
subplot(3,1,2);
plot(ix, iy, 'b.'); hold on;
p = polyfit(ix, iy, 1);
xx = linspace(min(ix), max(ix), 2);
plot(xx, p(1)*xx + p(2), 'r--', 'LineWidth', 1.5);
axis equal; grid on;
xlabel('a \delta i_x [m]'); ylabel('a \delta i_y [m]');
legend('Osc','Mean','Location','best');

% c) delta lambda vs delta a
x_lambda = a0 * rQNS(:,2,2);
y_da = a0 * rQNS(:,1,2);
subplot(3,1,3);
plot(x_lambda, y_da, 'b.'); hold on;
p = polyfit(x_lambda, y_da, 1);
xx = linspace(min(x_lambda), max(x_lambda), 2);
plot(xx, p(1)*xx + p(2), 'r--', 'LineWidth', 1.5);
axis equal; grid on;
xlabel('a \delta \lambda [m]'); ylabel('a \delta a [m]');
legend('Osc','Mean','Location','best');



%% 7.
qnsroes0 = utils.OE2rQNSOE(chief_OE, deputy_OE);
qnsroes0_m = chief_OE(1)*utils.OE2rQNSOE(chief_OE, deputy_OE);
qnsroes1 = qnsroes0 - [0 0 0 0 qnsroes0(5) 0];
deputy_OE_new = utils.rQNSOE2OE(chief_OE, qnsroes1);
disp(deputy_OE_new(2:end))


%% 8. STM-only vector plots for two initial conditions
%-----------------------------------------------------------
% Use the closed-form J2 STM to propagate the initial relative QNS for:
%  - IC set 2 (original QNSROE)
%  - IC set 6 (with delta_i_x zeroed)
% Then generate two figures (11 and 12) with three STM-only subplots each:
%   a) a*delta_e_x vs a*delta_e_y
%   b) a*delta_i_x vs a*delta_i_y
%   c) a*delta_lambda vs a*delta_a

% Unpack chief and constants
a0       = chief_OE(1);
e0       = chief_OE(2);
i0       = chief_OE(3);
omega0   = chief_OE(5);
n        = sqrt(mu/a0^3);
eta      = sqrt(1 - e0^2);
p        = a0*(1 - e0^2);
kappa    = 1.5 * J2 * n * (Re/p)^2;
P        = 0.5*(1 - 5*cos(i0)^2);
Qmat     = 0.25*(5*cos(i0)^2 - 1);
S        = 0.25*sin(i0)*cos(i0);
Tmat     = Qmat;
F        = 1/eta;
G        = 1/eta^3;
dot_omega= kappa * Qmat;

% Define initial normalized QNS vectors
qns0_row = utils.OE2rQNSOE(chief_OE, deputy_OE);
qns1_row = qns0_row;      % copy for zeroed delta_i_x
qns1_row(5) = 0;
qns0 = qns0_row(:);
qns1 = qns1_row(:);

% Allocate
N      = length(tspan);
anal0  = zeros(N,6);
anal1  = zeros(N,6);

% STM propagation
for k = 1:N
    tau = tspan(k);
    wf  = omega0 + dot_omega * tau;
    exi = e0*cos(omega0);
    eyi = e0*sin(omega0);
    exf = e0*cos(wf);
    eyf = e0*sin(wf);
    Phi = eye(6);
    Phi(2,1) = -(1.5*n + 3.5*kappa*eta*P)*tau;
    Phi(2,3) =  kappa*exi*F*G*P *tau;
    Phi(2,4) =  kappa*eyi*F*G*P *tau;
    Phi(2,5) = -kappa*F*S     *tau;
    Phi(3,1) =  3.5*kappa*eyf*Qmat *tau;
    Phi(3,3) =  cos(dot_omega*tau)-4*kappa*exi*eyf*G*Qmat *tau;
    Phi(3,4) = -sin(dot_omega*tau)-4*kappa*eyi*eyf*G*Qmat *tau;
    Phi(3,5) =  5*kappa*eyf*S        *tau;
    Phi(4,1) = -3.5*kappa*exf*Qmat *tau;
    Phi(4,3) =  sin(dot_omega*tau)+4*kappa*exi*exf*G*Qmat *tau;
    Phi(4,4) =  cos(dot_omega*tau)+4*kappa*eyi*exf*G*Qmat *tau;
    Phi(4,5) = -5*kappa*exf*S        *tau;
    Phi(6,1) =  3.5*kappa*S         *tau;
    Phi(6,3) = -4*kappa*exi*G*S     *tau;
    Phi(6,4) = -4*kappa*eyi*G*S     *tau;
    Phi(6,5) =  2*kappa*Tmat        *tau;
    anal0(k,:) = (a0*(Phi*qns0))';
    anal1(k,:) = (a0*(Phi*qns1))';
end

% STM-only plots
% Figure 11: original ICs
figure(12); clf; set(gcf,'Color','w');
subplot(3,1,1);
plot(anal0(:,3), anal0(:,4), 'b-'); axis equal; grid on;
xlabel('a \delta e_x [m]'); ylabel('a \delta e_y [m]');
subplot(3,1,2);
plot(anal0(:,5), anal0(:,6), 'b-'); axis equal; grid on;
xlabel('a \delta i_x [m]'); ylabel('a \delta i_y [m]');
subplot(3,1,3);
plot(anal0(:,2), anal0(:,1), 'b-'); axis equal; grid on;
xlabel('a \delta \lambda [m]'); ylabel('a \delta a [m]');

% Figure 12: delta_i_x zeroed ICs
figure(13); clf; set(gcf,'Color','w');
subplot(3,1,1);
plot(anal1(:,3), anal1(:,4), 'b-'); axis equal; grid on;
xlabel('a \delta e_x [m]'); ylabel('a \delta e_y [m]');
subplot(3,1,2);
plot(anal1(:,5), anal1(:,6), 'b-'); axis equal; grid on;
xlabel('a \delta i_x [m]'); ylabel('a \delta i_y [m]');
subplot(3,1,3);
plot(anal1(:,2), anal1(:,1), 'b-'); axis equal; grid on;
xlabel('a \delta \lambda [m]'); ylabel('a \delta a [m]');
