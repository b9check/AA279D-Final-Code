%% 1.3
close all
clear all
% Constants
mu = 3.986004418e14;  
Re = 6378137;          
J2 = 1.08262668e-3;   

% ICs
chief_OE = [6771e3; 0.0005; deg2rad(51.64); deg2rad(257); deg2rad(0); deg2rad(30)];
deputy_OE = [6771e3; 0.0015; deg2rad(52.14); deg2rad(257.5); deg2rad(0.5); deg2rad(25)];
chief_mean0 = utils.getMeanInitOE(chief_OE, chief_OE(2), mu, Re, J2);
deputy_mean0 = utils.getMeanInitOE(deputy_OE, deputy_OE(2), mu, Re, J2);

% QNS IC's (b)
deputy_rQNSOE_2 = [0, 100, 50, 100, 30, 200]/chief_OE(1);
deputy_OE_2 = utils.rQNSOE2OE(chief_OE, deputy_rQNSOE_2);
deputy2_mean0 = utils.getMeanInitOE(deputy_OE_2, deputy_OE_2(2), mu, Re, J2);


% Propagate chief with Keplerian propagator, no J2
odefun = @(t, state) propagators.KeplerianPropagator.getStatedot(t, state, mu, Re, J2, false);
odefun_J2 = @(t, state) propagators.KeplerianPropagator.getStatedot(t, state, mu, Re, J2, true);
odefun_mean = @(t, state) propagators.KeplerianPropagator.getMeanStatedot(t, state, mu, Re, J2, false);
odefun_mean_J2 = @(t, state) propagators.KeplerianPropagator.getMeanStatedot(t, state, mu, Re, J2, true);
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

% Define time span
T = 2*pi*sqrt(chief_OE(1)^3 / mu); % Period
tspan = 0:1:15*T;  

% Simulate chief with all 4 odefuns
[t_c, x_c] = ode113(odefun, tspan, chief_OE, opts);
[t_c_m, x_c_m] = ode113(odefun_mean, tspan, chief_OE, opts);
[t_c_J2, x_c_J2] = ode113(odefun_J2, tspan, chief_OE, opts);
[t_c_mJ2, x_c_mJ2] = ode113(odefun_mean_J2, tspan, chief_mean0, opts);

% Simulate deputy with all 4 odefuns for IC's from (a)
[t_d1, x_d1] = ode113(odefun, tspan, deputy_OE, opts);
[t_d1_m, x_d1_m] = ode113(odefun_mean, tspan, deputy_OE, opts);
[t_d1_J2, x_d1_J2] = ode113(odefun_J2, tspan, deputy_OE, opts);
[t_d1_mJ2, x_d1_mJ2] = ode113(odefun_mean_J2, tspan, deputy_mean0, opts);

% Simulate deputy with all 4 odefuns for IC's from (b)
state0 = deputy_OE_2;
[t_d2, x_d2] = ode113(odefun, tspan, deputy_OE_2, opts);
[t_d2_m, x_d2_m] = ode113(odefun_mean, tspan, deputy_OE_2, opts);
[t_d2_J2, x_d2_J2] = ode113(odefun_J2, tspan, deputy_OE_2, opts);
[t_d2_mJ2, x_d2_mJ2] = ode113(odefun_mean_J2, tspan, deputy2_mean0, opts);

% Case 1: IC (a), no J2
QNSOEs1 = zeros(length(tspan), 6);
QNSROEs1 = zeros(length(tspan), 6);
meanQNSOEs1 = zeros(length(tspan), 6);
meanQNSROEs1 = zeros(length(tspan), 6);

% Case 2: IC (a), with J2
QNSOEs2 = zeros(length(tspan), 6);
QNSROEs2 = zeros(length(tspan), 6);
meanQNSOEs2 = zeros(length(tspan), 6);
meanQNSROEs2 = zeros(length(tspan), 6);

% Case 3: IC (b), no J2
QNSOEs3 = zeros(length(tspan), 6);
QNSROEs3 = zeros(length(tspan), 6);
meanQNSOEs3 = zeros(length(tspan), 6);
meanQNSROEs3 = zeros(length(tspan), 6);

% Case 4: IC (b), with J2
QNSOEs4 = zeros(length(tspan), 6);
QNSROEs4 = zeros(length(tspan), 6);
meanQNSOEs4 = zeros(length(tspan), 6);
meanQNSROEs4 = zeros(length(tspan), 6);


for i = 1:length(tspan)
    % Extract chief OE's
    chief_OE_i = x_c(i,:);
    chief_OE_J2_i = x_c_J2(i,:);
    chief_OE_im = x_c_m(i,:);
    chief_OE_J2_im = x_c_mJ2(i,:);
    
    % Extract deputy OE's for IC's (a)
    deputy_OE_i = x_d1(i,:);
    deputy_OE_J2_i = x_d1_J2(i,:);
    deputy_OE_im = x_d1_m(i,:);
    deputy_OE_J2_im = x_d1_mJ2(i,:);

    % Extract deputy OE's for IC's (b)
    deputy2_OE_i = x_d2(i,:);
    deputy2_OE_J2_i = x_d2_J2(i,:);
    deputy2_OE_im = x_d2_m(i,:);
    deputy2_OE_J2_im = x_d2_mJ2(i,:);
    
    % Case 1: IC (a), no J2
    QNSOEs1(i, :) = utils.OE2QNSOE(deputy_OE_i);
    QNSROEs1(i, :) = utils.OE2rQNSOE(chief_OE_i, deputy_OE_i);
    meanQNSOEs1(i,:) = utils.OE2QNSOE(deputy_OE_im);
    meanQNSROEs1(i,:) = utils.OE2rQNSOE(chief_OE_im, deputy_OE_im);

    % Case 2: IC (a), with J2
    QNSOEs2(i, :) = utils.OE2QNSOE(deputy_OE_J2_i);
    QNSROEs2(i, :) = utils.OE2rQNSOE(chief_OE_J2_i, deputy_OE_J2_i);
    meanQNSOEs2(i,:) = utils.OE2QNSOE(deputy_OE_J2_im);
    meanQNSROEs2(i,:) = utils.OE2rQNSOE(chief_OE_J2_im, deputy_OE_J2_im);

    % Case 3: IC (b), no J2
    QNSOEs3(i, :) = utils.OE2QNSOE(deputy2_OE_i);
    QNSROEs3(i, :) = utils.OE2rQNSOE(chief_OE_i, deputy2_OE_i);
    meanQNSOEs3(i,:) = utils.OE2QNSOE(deputy2_OE_im);
    meanQNSROEs3(i,:) = utils.OE2rQNSOE(chief_OE_im, deputy2_OE_im);

    % Case 4: IC (b), with J2
    QNSOEs4(i, :) = utils.OE2QNSOE(deputy2_OE_J2_i);
    QNSROEs4(i, :) = utils.OE2rQNSOE(chief_OE_J2_i, deputy2_OE_J2_i);
    meanQNSOEs4(i,:) = utils.OE2QNSOE(deputy2_OE_J2_im);
    meanQNSROEs4(i,:) = utils.OE2rQNSOE(chief_OE_J2_im, deputy2_OE_J2_im);
    
end

% Unwrap lambda
QNSOEs1(:,2)       = unwrap(QNSOEs1(:,2));
QNSROEs1(:,2)      = unwrap(QNSROEs1(:,2));
QNSOEs2(:,2)       = unwrap(QNSOEs2(:,2));
QNSROEs2(:,2)      = unwrap(QNSROEs2(:,2));
QNSOEs3(:,2)       = unwrap(QNSOEs3(:,2));
QNSROEs3(:,2)      = unwrap(QNSROEs3(:,2));
QNSOEs4(:,2)       = unwrap(QNSOEs4(:,2));
QNSROEs4(:,2)      = unwrap(QNSROEs4(:,2));
meanQNSOEs1(:,2)   = unwrap(meanQNSOEs1(:,2));
meanQNSROEs1(:,2)  = unwrap(meanQNSROEs1(:,2));
meanQNSOEs2(:,2)   = unwrap(meanQNSOEs2(:,2));
meanQNSROEs2(:,2)  = unwrap(meanQNSROEs2(:,2));
meanQNSOEs3(:,2)   = unwrap(meanQNSOEs3(:,2));
meanQNSROEs3(:,2)  = unwrap(meanQNSROEs3(:,2));
meanQNSOEs4(:,2)   = unwrap(meanQNSOEs4(:,2));
meanQNSROEs4(:,2)  = unwrap(meanQNSROEs4(:,2));


% ---- Plot QNSOEs of Deputy ----
t_plot = tspan;
figure(1)
labels = {'a', '\lambda', 'e_x', 'e_y', 'i_x', 'i_y'};
for j = 1:6
    subplot(6,1,j)
    plot(t_plot, QNSOEs1(:,j), 'b', t_plot, meanQNSOEs1(:,j), 'r--')
    ylabel(['\delta ', labels{j}])
    legend('Osculating','Mean')
    grid on
    if j < 6
        set(gca, 'XTickLabel', [])
    else
        xlabel('Time (s)')
    end
end
sgtitle("QNSOE's with IC (a), no J2")

figure(2)
for j = 1:6
    subplot(6,1,j)
    plot(t_plot, QNSROEs1(:,j), 'b', t_plot, meanQNSROEs1(:,j), 'r--')
    ylabel(['\delta ', labels{j}])
    legend('Osculating','Mean')
    grid on
    if j < 6
        set(gca, 'XTickLabel', [])
    else
        xlabel('Time (s)')
    end
end
sgtitle("QNSROE's with IC (a), no J2")


figure(3)
for j = 1:6
    subplot(6,1,j)
    plot(t_plot, QNSOEs2(:,j), 'b', t_plot, meanQNSOEs2(:,j), 'r--')
    ylabel(['\delta ', labels{j}])
    legend('Osculating','Mean')
    grid on
    if j < 6
        set(gca, 'XTickLabel', [])
    else
        xlabel('Time (s)')
    end
end
sgtitle("QNSOE's with IC (a), with J2")

figure(4)
for j = 1:6
    subplot(6,1,j)
    plot(t_plot, QNSROEs2(:,j), 'b', t_plot, meanQNSROEs2(:,j), 'r--')
    ylabel(['\delta ', labels{j}])
    legend('Osculating','Mean')
    grid on
    if j < 6
        set(gca, 'XTickLabel', [])
    else
        xlabel('Time (s)')
    end
end
sgtitle("QNSROE's with IC (a), with J2")

figure(5)
for j = 1:6
    subplot(6,1,j)
    plot(t_plot, QNSOEs3(:,j), 'b', t_plot, meanQNSOEs3(:,j), 'r--')
    ylabel(['\delta ', labels{j}])
    legend('Osculating','Mean')
    grid on
    if j < 6
        set(gca, 'XTickLabel', [])
    else
        xlabel('Time (s)')
    end
end
sgtitle("QNSOE's with IC (b), no J2")

figure(6)
for j = 1:6
    subplot(6,1,j)
    plot(t_plot, QNSROEs3(:,j), 'b', t_plot, meanQNSROEs3(:,j), 'r--')
    ylabel(['\delta ', labels{j}])
    legend('Osculating','Mean')
    grid on
    if j < 6
        set(gca, 'XTickLabel', [])
    else
        xlabel('Time (s)')
    end
end
sgtitle("QNSROE's with IC (b), no J2")

figure(7)
for j = 1:6
    subplot(6,1,j)
    plot(t_plot, QNSOEs4(:,j), 'b', t_plot, meanQNSOEs4(:,j), 'r--')
    ylabel(['\delta ', labels{j}])
    legend('Osculating','Mean')
    grid on
    if j < 6
        set(gca, 'XTickLabel', [])
    else
        xlabel('Time (s)')
    end
end
sgtitle("QNSOE's with IC (b), with J2")

figure(8)
for j = 1:6
    subplot(6,1,j)
    plot(t_plot, QNSROEs4(:,j), 'b', t_plot, meanQNSROEs4(:,j), 'r--')
    ylabel(['\delta ', labels{j}])
    legend('Osculating','Mean')
    grid on
    if j < 6
        set(gca, 'XTickLabel', [])
    else
        xlabel('Time (s)')
    end
end
sgtitle("QNSROE's with IC (b), with J2")



%% 1.4
% Compute RTN trajectories for IC (a) cases: no J2 vs J2
RTN_unp = zeros(length(tspan),3);
RTN_J2  = zeros(length(tspan),3);

for k = 1:length(tspan)
    % Chief ECI (row - column)
    chief_ECI_unp = utils.OE2ECI(x_c(k,:)', mu);
    chief_ECI_J2  = utils.OE2ECI(x_c_J2(k,:)', mu);

    % Deputy ECI
    dep_ECI_unp = utils.OE2ECI(x_d1(k,:)', mu);
    dep_ECI_J2  = utils.OE2ECI(x_d1_J2(k,:)', mu);

    % Compute full RTN state vector, then extract position
    fullRTN_unp = utils.ECI2RTN(chief_ECI_unp, dep_ECI_unp);
    RTN_unp(k,:) = fullRTN_unp(1:3)';    % R,T,N position

    fullRTN_J2 = utils.ECI2RTN(chief_ECI_J2, dep_ECI_J2);
    RTN_J2(k,:)  = fullRTN_J2(1:3)';
end


% Plot relative position in RTN (m)
figure(9); clf; set(gcf, 'Color', 'w');

% 3D RTN trajectory
subplot(2,2,1)
plot3( RTN_J2(:,2),  RTN_J2(:,1),  RTN_J2(:,3),  'b', ...
       RTN_unp(:,2), RTN_unp(:,1), RTN_unp(:,3), 'r--' );
xlabel('T [m]'); ylabel('R [m]'); zlabel('N [m]');
legend('J2','Unperturbed');

% T–R plane
subplot(2,2,2)
plot( RTN_J2(:,2),  RTN_J2(:,1),  'b', ...
      RTN_unp(:,2), RTN_unp(:,1), 'r--' );
axis equal; grid on
xlabel('T [m]'); ylabel('R [m]');

% N–R plane
subplot(2,2,3)
plot( RTN_J2(:,3),  RTN_J2(:,1),  'b', ...
      RTN_unp(:,3), RTN_unp(:,1), 'r--' );
axis equal; grid on
xlabel('N [m]'); ylabel('R [m]');

% T–N plane
subplot(2,2,4)
plot( RTN_J2(:,2),  RTN_J2(:,3),  'b', ...
      RTN_unp(:,2), RTN_unp(:,3), 'r--' );
axis equal; grid on
xlabel('T [m]'); ylabel('N [m]');



%% 1.5 — Plot perturbed vs propagated‐mean relative QNS‐ROEs (Case a, with J2)
a0 = chief_OE(1);
N  = numel(tspan);

% (1) Rebuild the osculating‐J2 QNSROE series
Qosc = zeros(N,6);
for k = 1:N
    ce_osc = utils.ECI2OE(utils.OE2ECI(x_c_J2(k,:)',mu),mu);
    de_osc = utils.ECI2OE(utils.OE2ECI(x_d1_J2(k,:)',mu),mu);
    Qosc(k,:) = utils.OE2rQNSOE(ce_osc, de_osc);
end
Qosc(:,2) = unwrap(Qosc(:,2));

% (2) Build the mean‐J2 QNSROE series from the numeric‐mean propagation
Qmean = zeros(N,6);
for k = 1:N
    ce_m = utils.ECI2OE(utils.OE2ECI(x_c_mJ2(k,:)',mu),mu);
    de_m = utils.ECI2OE(utils.OE2ECI(x_d1_mJ2(k,:)',mu),mu);
    Qmean(k,:) = utils.OE2rQNSOE(ce_m, de_m);
end
Qmean(:,2) = unwrap(Qmean(:,2));

% (3) Plot each pair in its own figure:

% Figure 10: δe_x vs δe_y
figure(10); clf; set(gcf,'Color','w');
plot(Qosc(:,3), Qosc(:,4), 'b-', Qmean(:,3), Qmean(:,4), 'r--','LineWidth',1.4);
axis equal; grid on;
xlabel('\delta e_x'); ylabel('\delta e_y');
title('Relative Eccentricity Vector (J_2)');
legend('Osculating','Mean','Location','best');

% Figure 11: δi_x vs δi_y
figure(11); clf; set(gcf,'Color','w');
plot(Qosc(:,5), Qosc(:,6), 'b-', Qmean(:,5), Qmean(:,6), 'r--','LineWidth',1.4);
axis equal; grid on;
xlabel('\delta i_x'); ylabel('\delta i_y');
title('Relative Inclination Vector (J_2)');
legend('Osculating','Mean','Location','best');

% Figure 12: δλ vs δa (both in meters)
delta_la = Qosc(:,2) * a0;
delta_a  = Qosc(:,1) * a0;
mean_la  = Qmean(:,2) * a0;
mean_a   = Qmean(:,1) * a0;

figure(12); clf; set(gcf,'Color','w');
plot(delta_la, delta_a, 'b-', mean_la, mean_a, 'r--','LineWidth',1.4);
axis equal; grid on;
xlabel('\delta \lambda [m]'); ylabel('\delta a [m]');
title('Relative Mean Longitude vs Semi-major Axis (J_2)');
legend('Osculating','Mean','Location','best');